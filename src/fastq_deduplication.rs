use fastq::{each_zipped, Parser, Record};
use hashbrown::{HashMap, HashSet};
use indicatif::ParallelProgressIterator;
use log::{debug, info, warn};
use rand::prelude::*;
use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use std::ops::Add;
use std::path::Path;
use std::prelude::rust_2021::*;
use std::{iter::Iterator, str::FromStr};
use tempfile::tempdir;
use twox_hash::xxh3::hash64_with_seed;

use crate::utils::{get_fastq_reader_file_handles, get_fastq_writer_file_handles, write_records};

#[derive(Debug, Serialize, Deserialize, PartialEq, Clone, Copy)]
pub struct FastqReadDeduplicationStats {
    read_pairs_total: u64,
    read_pairs_duplicated: u64,
    read_pairs_unique: u64,
}

impl Add for FastqReadDeduplicationStats {
    type Output = Self;

    fn add(self, other: Self) -> Self {
        Self {
            read_pairs_total: self.read_pairs_total + other.read_pairs_total,
            read_pairs_duplicated: self.read_pairs_duplicated + other.read_pairs_duplicated,
            read_pairs_unique: self.read_pairs_unique + other.read_pairs_unique,
        }
    }
}

impl FastqReadDeduplicationStats {
    pub fn new() -> Self {
        Self {
            read_pairs_total: 0,
            read_pairs_duplicated: 0,
            read_pairs_unique: 0,
        }
    }

    pub fn get_read_pairs_total(&self) -> u64 {
        self.read_pairs_total
    }

    pub fn get_read_pairs_duplicated(&self) -> u64 {
        self.read_pairs_duplicated
    }

    pub fn get_read_pairs_unique(&self) -> u64 {
        self.read_pairs_unique
    }

    pub fn get_read_pairs_duplicated_percentage(&self) -> f64 {
        self.read_pairs_duplicated as f64 / self.read_pairs_total as f64
    }

    pub fn get_read_pairs_unique_percentage(&self) -> f64 {
        self.read_pairs_unique as f64 / self.read_pairs_total as f64
    }
}

struct ShardDuplicates {
    fq1: String,
    fq2: String,
    shard_inner_duplicate_positions: Vec<usize>,
    shard_reads_seen: HashSet<u64>,
    shard_duplicate_read_hashes: HashSet<u64>,
}

pub struct FastqDeduplicator {
    paths: Vec<(String, String)>,
    output_paths: Vec<(String, String)>,
    duplicates: HashMap<(String, String), Vec<usize>>,
    shuffle_shard_order: bool,
    compress_output: bool,
}

impl FastqDeduplicator {
    pub fn new(
        paths: Vec<(String, String)>,
        output_paths: Option<Vec<(String, String)>>,
        shuffle_shard_order: bool,
    ) -> Self
where {
        // If no output paths are provided, create them from the input paths
        let output_paths: Vec<(String, String)> = match output_paths {
            Some(paths) => paths
                .into_iter()
                .map(|(p1, p2)| (p1.into(), p2.into()))
                .collect(),
            None => paths
                .iter()
                .map(|(r1, r2)| {
                    let o1 = r1.replace(".fastq.gz", "_deduplicated.fastq.gz");
                    let o2 = r2.replace(".fastq.gz", "_deduplicated.fastq.gz");
                    (o1, o2)
                })
                .collect(),
        };

        if paths.len() != output_paths.len() {
            panic!("The number of input paths must match the number of output paths");
        }

        // Check if the output paths are compressed
        let compress_output = match output_paths[0].0.split('.').last() {
            Some("gz") => true,
            _ => false,
        };

        Self {
            paths,
            output_paths,
            duplicates: HashMap::new(),
            shuffle_shard_order,
            compress_output,
        }
    }

    fn get_compression(&self) -> niffler::Format {
        if self.compress_output {
            niffler::Format::Gzip
        } else {
            niffler::Format::No
        }
    }

    fn get_compression_level(&self) -> Option<niffler::Level> {
        if self.compress_output {
            Some(niffler::Level::Five)
        } else {
            None
        }
    }

    fn get_hashed_sequences_by_shard(&self) -> Result<Vec<ShardDuplicates>, std::io::Error> {
        info!("Identifying duplicates by shard");

        // Iterate over the paths in parallel using rayon par_iter
        let shard_duplicates = self
            .paths
            .par_iter()
            .map(|(r1, r2)| {
                let mut file_handles = get_fastq_reader_file_handles(vec![r1, r2]).unwrap();
                let mut duplicate_read_positions = Vec::new();
                let mut reads_seen = HashSet::new();
                let mut read_count = 0;

                each_zipped(
                    file_handles.remove(0),
                    file_handles.remove(0),
                    |r1, r2| match (r1, r2) {
                        (Some(rec1), Some(rec2)) => {

                            // Get the hash of the sequence
                            let sequences = [rec1.seq(), rec2.seq()].concat();
                            let sequences_hashed = hash64_with_seed(&sequences, 42);

                            // Check if the hash is in the reads_seen set
                            if reads_seen.contains(&sequences_hashed) {
                                duplicate_read_positions.push(read_count);
                            } else {
                                reads_seen.insert(sequences_hashed);
                            }

                            read_count += 1;

                            (true, true)
                        }
                        _ => (false, false),
                    },
                )
                .expect("Error reading fq");

                ShardDuplicates {
                    fq1: r1.to_string(),
                    fq2: r2.to_string(),
                    shard_inner_duplicate_positions: duplicate_read_positions,
                    shard_reads_seen: reads_seen,
                    shard_duplicate_read_hashes: HashSet::new(),
                }
            })
            .collect::<Vec<ShardDuplicates>>();

        Ok(shard_duplicates)
    }

    fn unique_reads_identify(&mut self) -> Result<Vec<ShardDuplicates>, std::io::Error> {
        // Get the hashed sequences by shard
        let shard_duplicates = self.get_hashed_sequences_by_shard()?;

        // Shuffle the order of the shards
        let mut shard_duplicates = if self.shuffle_shard_order {
            let mut rng = thread_rng();
            let mut shard_duplicates = shard_duplicates;
            shard_duplicates.shuffle(&mut rng);
            shard_duplicates
        } else {
            shard_duplicates
        };

        // See if any read is duplicated in any other shard
        // Exclude the current shard from the comparison
        // Extend the reads_seen set with the current shard's reads
        let mut reads_seen = HashSet::new();

        for mut shard in shard_duplicates.iter_mut() {
            let shard_reads_seen = &shard.shard_reads_seen;

            shard.shard_duplicate_read_hashes = shard_reads_seen
                .intersection(&reads_seen)
                .cloned()
                .collect();

            reads_seen.extend(shard_reads_seen);
        }

        Ok(shard_duplicates)
    }

    pub fn write_unique_reads(&mut self) -> Result<FastqReadDeduplicationStats, std::io::Error> {
        let shard_duplicates = self.unique_reads_identify()?;
        
        info!("Writing unique reads");
        // Iterate over the paths in parallel using rayon par_iter
        let stats = self
            .paths
            .par_iter()
            .zip(self.output_paths.par_iter())
            .zip(shard_duplicates.into_par_iter())
            .map(|(((r1, r2), (o1, o2)), shard)| {
                let mut file_handles = get_fastq_reader_file_handles(vec![r1, r2]).unwrap();
                let mut output_file_handles = get_fastq_writer_file_handles(
                    vec![o1, o2],
                    self.get_compression(),
                    self.get_compression_level(),
                )
                .expect("Error opening output file");

                let mut dedup_stats = FastqReadDeduplicationStats::new();

                each_zipped(
                    file_handles.remove(0),
                    file_handles.remove(0),
                    |r1, r2| match (r1, r2) {
                        (Some(rec1), Some(rec2)) => {
                            dedup_stats.read_pairs_total += 1;
                            let read_position = dedup_stats.read_pairs_total;

                            // Get the hash of the sequence
                            let sequences = [rec1.seq(), rec2.seq()].concat();
                            let sequences_hashed = hash64_with_seed(&sequences, 42);

                            // Check if the hash is in the reads_seen set
                            if (!shard
                                .shard_duplicate_read_hashes
                                .contains(&sequences_hashed))
                                && (!shard
                                    .shard_inner_duplicate_positions
                                    .contains(&(read_position as usize)))
                            {
                                write_records(vec![rec1, rec2], &mut output_file_handles);
                                dedup_stats.read_pairs_unique += 1;
                            } else {
                                dedup_stats.read_pairs_duplicated += 1;
                            }

                            (true, true)
                        }
                        _ => (false, false),
                    },
                )
                .expect("Error reading fq");

                dedup_stats
            })
            .reduce(|| FastqReadDeduplicationStats::new(), |a, b| a + b);

        Ok(stats)
    }
}

// Test the deduplication class with a single pair of fastq files
#[cfg(test)]
mod tests {
    use super::*;
    use std::fs::File;
    use std::io::prelude::*;
    use tempfile::tempdir;

    use crate::fastq_deduplication::FastqDeduplicator;
    use crate::fastq_deduplication::FastqReadDeduplicationStats;

    #[test]
    fn test_deduplication() {
        let fq1 = "tests/fastq_deduplicate/duplicated_1.fastq.gz".to_string();
        let fq2 = "tests/fastq_deduplicate/duplicated_2.fastq.gz".to_string();
        let fq1_2 = "tests/fastq_deduplicate/duplicated2_1.fastq.gz".to_string();
        let fq2_2 = "tests/fastq_deduplicate/duplicated2_2.fastq.gz".to_string();

        // Create a temporary directory
        let temp_dir = tempdir().unwrap();

        // Create a temporary file
        let mut temp_file1 = tempfile::NamedTempFile::new_in(temp_dir.path())
            .unwrap()
            .path()
            .as_os_str()
            .to_str()
            .unwrap()
            .to_string();
        let mut temp_file2 = tempfile::NamedTempFile::new_in(temp_dir.path())
            .unwrap()
            .path()
            .as_os_str()
            .to_str()
            .unwrap()
            .to_string();

        let infiles = vec![(fq1, fq2), (fq1_2, fq2_2)];
        let outfiles = vec![(temp_file1, temp_file2)];

        // Initialise deduplication class
        let mut deduplication = FastqDeduplicator::new(infiles, None, false);
        let deduplication_stats = deduplication.write_unique_reads().unwrap();

        println!("deduplication:{:?}", deduplication_stats);

        // Check the number of unique reads
        assert_eq!(
            deduplication_stats.read_pairs_unique,
            982,
        );
    }
}
