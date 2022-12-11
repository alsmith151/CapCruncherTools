use bloom::{BloomFilter, ASMS};
use fastq::{each_zipped, Parser, Record};
use hashbrown::{HashMap, HashSet};
use indicatif::ParallelProgressIterator;
use itertools::Itertools;
use log::{debug, info, warn};
use rand::prelude::*;
use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use std::iter::Iterator;
use std::ops::Add;
use std::path::Path;
use std::prelude::rust_2021::*;

use crate::utils::{get_fastq_reader_file_handles, get_fastq_writer_file_handles, write_records};

#[derive(Debug, Serialize, Deserialize)]
pub struct DeduplicationStats {
    read_pairs_total: u64,
    read_pairs_duplicated: u64,
    read_pairs_unique: u64,
}

impl DeduplicationStats {
    pub fn new(read_pairs_total: u64, read_pairs_unique: u64) -> Self {
        Self {
            read_pairs_total,
            read_pairs_unique,
            read_pairs_duplicated: read_pairs_total - read_pairs_unique,
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
}

impl Add for DeduplicationStats {
    type Output = Self;
    fn add(self, rhs: Self) -> Self::Output {
        Self {
            read_pairs_total: self.read_pairs_total + rhs.read_pairs_total,
            read_pairs_duplicated: self.read_pairs_duplicated + rhs.read_pairs_duplicated,
            read_pairs_unique: self.read_pairs_unique + rhs.read_pairs_unique,
        }
    }
}

pub struct FastqDeduplicator {
    paths: Vec<(String, String)>,
    output_paths: Vec<(String, String)>,
    duplicates: HashMap<(String, String), Vec<usize>>,
    expected_num_items: u32,
    error_rate: f32,
    shuffle_shard_order: bool,
    compress_output: bool,
}

impl FastqDeduplicator {
    pub fn new(
        paths: Vec<(String, String)>,
        output_paths: Option<Vec<(String, String)>>,
        error_rate: Option<f32>,
        expected_num_items: Option<u32>,
        shuffle_shard_order: bool,
    ) -> Self {
        
        let output_paths = match output_paths {
            Some(paths) => paths,
            None => paths
                .iter()
                .map(|(r1, r2)| {
                    let o1 = r1.replace(".fastq.gz", "_deduplicated.fastq.gz");
                    let o2 = r2.replace(".fastq.gz", "_deduplicated.fastq.gz");
                    (o1, o2)
                })
                .collect(),
        };

        let error = match error_rate {
            Some(error) => error,
            None => 1e-4,
        };

        let expected_num_items = match expected_num_items {
            Some(expected_num_items) => expected_num_items,
            None => (paths.len() as u32 * 1e6 as u32) as u32,
        };

        let compress_output = match output_paths[0].0.split('.').last() {
            Some("gz") => true,
            _ => false,
        };

        Self {
            paths,
            output_paths,
            duplicates: HashMap::new(),
            expected_num_items: expected_num_items,
            error_rate: error,
            shuffle_shard_order,
            compress_output,
        }
    }

    pub fn identify_duplicates(&mut self) -> Result<(), std::io::Error> {
        
        info!("Identifying duplicates");

        // Shuffle the paths if shuffle_shard_order is true
        let mut paths = self.paths.clone();
        if self.shuffle_shard_order {
            let rng = &mut thread_rng();
            paths.shuffle(rng);
        }

        // Create a bloom filter with the expected number of items and the false positive rate
        let expected_num_items = self.expected_num_items;
        let false_positive_rate = self.error_rate;
        let mut items_seen = BloomFilter::with_rate(false_positive_rate, expected_num_items);

        // Create a HashMap to store the duplicate read positions
        let mut duplicates = HashMap::new();

        // Iterate over the paths and identify the duplicate read positions
        for (r1, r2) in paths.iter() {

            info!("Processing {} and {}", r1, r2);

            let mut file_handles = get_fastq_reader_file_handles(vec![r1, r2]).unwrap();
            let mut duplicate_read_positions = Vec::new();
            let mut read_count = 0;

            // Iterate over the read pairs and identify the duplicate read positions
            each_zipped(
                file_handles.remove(0),
                file_handles.remove(0),
                |r1, r2| match (r1, r2) {
                    (Some(rec1), Some(rec2)) => {

                        // Check if 100_000 reads have been processed
                        if read_count % 100_000 == 0 {
                            info!("Processed {} reads", read_count);
                        }

                        let sequences = [rec1.seq(), rec2.seq()].concat();
                        match items_seen.insert(&sequences) {
                            true => (),
                            false => {
                                duplicate_read_positions.push(read_count);
                            }
                        }

                        read_count += 1;
                        (true, true)
                    }
                    _ => (false, false),
                },
            )?;
            duplicates.insert((r1.to_owned(), r2.to_owned()), duplicate_read_positions);
        }

        self.duplicates = duplicates;
        Ok(())
    }

    pub fn remove_duplicate_sequences(&self) -> Result<DeduplicationStats, std::io::Error> {
        // Identify the output format
        let output_format = match self.compress_output {
            true => niffler::Format::Gzip,
            false => niffler::Format::No,
        };

        info!("Removing duplicate sequences");

        // Iterate over the paths and remove the duplicate read positions
        // and write the deduplicated reads to the output paths
        // Return the deduplication statistics
        let deduplication_stats = self
            .paths
            .par_iter()
            .zip(self.output_paths.to_owned())
            .map(|((r1, r2), (o1, o2))| {
                let mut file_handles = get_fastq_reader_file_handles(vec![r1, r2]).unwrap();
                let mut writer_handles =
                    get_fastq_writer_file_handles(vec![o1, o2], output_format, None)
                        .expect("Cannot open files for writing");

                let mut read_number: usize = 0;
                let mut duplicate_count: usize = 0;

                let fq_file_names_key = (r1.to_owned(), r2.to_owned());
                let duplicates = self.duplicates.get(&(fq_file_names_key)).unwrap();

                each_zipped(
                    file_handles.remove(0),
                    file_handles.remove(0),
                    |r1, r2| match (r1, r2) {
                        (Some(rec1), Some(rec2)) => {
                            if !duplicates.contains(&read_number) {
                                write_records(vec![rec1, rec2], &mut writer_handles);
                            } else {
                                duplicate_count += 1;
                            }

                            read_number += 1;
                            (true, true)
                        }
                        _ => (false, false),
                    },
                )
                .expect("Cannot read fastq files");

                info!(
                    "Removed {} duplicates from {} and {}",
                    duplicate_count, r1, r2
                );
                DeduplicationStats::new(read_number as u64, (read_number - duplicate_count) as u64)
            })
            .reduce(|| DeduplicationStats::new(0, 0), |acc, x| acc + x);

        Ok(deduplication_stats)
    }
}

// Test FastqDeduplicator
// Use files from tests/fastq_deduplicate
// Test identify and remove duplicates functions on copies of the same file
#[cfg(test)]
mod tests {
    use super::*;
    use tempfile;

    #[test]
    fn test_fastq_deduplicator() {
        let tmpdir = tempfile::tempdir().unwrap();

        let paths = vec![(
            "tests/fastq_deduplicate/duplicated_1.fastq.gz".to_owned(),
            "tests/fastq_deduplicate/duplicated_2.fastq.gz".to_owned(),
        )];

        let output_paths = vec![(
            tmpdir
                .path()
                .join("reads_1_dedup.fq.gz".to_owned())
                .to_str()
                .unwrap()
                .to_owned(),
            tmpdir
                .path()
                .join("reads_2_dedup.fq.gz".to_owned())
                .to_str()
                .unwrap()
                .to_owned(),
        )];

        let mut deduplicator = FastqDeduplicator::new(paths, Some(output_paths), Some(0.01), None, true);
        deduplicator.identify_duplicates().unwrap();
        let deduplication_stats = deduplicator.remove_duplicate_sequences().unwrap();

        println!("{:?}", deduplication_stats);

        assert_eq!(deduplication_stats.get_read_pairs_total(), 1520);
        assert_eq!(deduplication_stats.get_read_pairs_unique(), 982);
        assert_eq!(deduplication_stats.get_read_pairs_duplicated(), 538);
    }

    #[test]
    fn test_fastq_deduplicator_multiple_files() {
        let tmpdir = tempfile::tempdir().unwrap();

        let paths = vec![
            (
                "tests/fastq_deduplicate/duplicated_1.fastq.gz".to_owned(),
                "tests/fastq_deduplicate/duplicated_2.fastq.gz".to_owned(),
            ),
            (
                "tests/fastq_deduplicate/duplicated2_1.fastq.gz".to_owned(),
                "tests/fastq_deduplicate/duplicated2_2.fastq.gz".to_owned(),
            ),
        ];

        let output_paths = vec![
            (
                tmpdir
                    .path()
                    .join("reads_1_dedup.fq.gz".to_owned())
                    .to_str()
                    .unwrap()
                    .to_owned(),
                tmpdir
                    .path()
                    .join("reads_2_dedup.fq.gz".to_owned())
                    .to_str()
                    .unwrap()
                    .to_owned(),
            ),
            (
                tmpdir
                    .path()
                    .join("reads2_1_dedup.fq.gz".to_owned())
                    .to_str()
                    .unwrap()
                    .to_owned(),
                tmpdir
                    .path()
                    .join("reads2_2_dedup.fq.gz".to_owned())
                    .to_str()
                    .unwrap()
                    .to_owned(),
            ),
        ];

        let mut deduplicator = FastqDeduplicator::new(paths, Some(output_paths), Some(0.01), None, true);
        deduplicator.identify_duplicates().unwrap();
        let deduplication_stats = deduplicator.remove_duplicate_sequences().unwrap();

        println!("{:?}", deduplication_stats);

        assert_eq!(deduplication_stats.get_read_pairs_total(), 3040);
        assert_eq!(deduplication_stats.get_read_pairs_unique(), 982);
        assert_eq!(deduplication_stats.get_read_pairs_duplicated(), 2058);
    }
}
