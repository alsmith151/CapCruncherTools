use fastq::{each_zipped, Parser, Record};
use hashbrown::{HashMap, HashSet};
use log::{debug, info, warn};
use niffler;
use rand::prelude::*;
use rayon::iter::{IntoParallelIterator, IntoParallelRefIterator, ParallelIterator, ParallelBridge};
use serde::{Deserialize, Serialize};
use std::iter::Iterator;
use std::ops::Add;
use std::path::Path;
use std::prelude::rust_2021::*;
use tempfile::tempdir;
use twox_hash::xxh3::hash64_with_seed;

fn get_fastq_reader_file_handles<P>(
    paths: Vec<P>,
) -> Result<Vec<Parser<Box<dyn std::io::Read>>>, std::io::Error>
where
    P: AsRef<Path>,
{
    let file_handles: Vec<_> = paths
        .iter()
        .map(|p| niffler::from_path(p.as_ref()).expect("Compression not recognised"))
        .map(|(fh, _fmt)| Parser::new(fh))
        .collect();

    Ok(file_handles)
}

fn get_fastq_writer_file_handles<P>(
    paths: Vec<P>,
    compression_format: niffler::compression::Format,
    compression_level: niffler::Level,
) -> Result<Vec<Box<dyn std::io::Write>>, std::io::Error>
where
    P: AsRef<Path>,
{
    let file_handles: Vec<_> = paths
        .iter()
        .map(|p| {
            niffler::to_path(p.as_ref(), compression_format, compression_level)
                .expect("Failed to create output file")
        })
        .collect();

    Ok(file_handles)
}

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

fn convert_fastq_to_hashmap<P: AsRef<Path>, S: ToString>(
    paths: Vec<P>,
    label: S,
) -> Result<(String, HashMap<u64, u64>, HashSet<u64>), std::io::Error> {
    let mut file_handles = get_fastq_reader_file_handles(paths)?;
    let mut sequences_unique: HashMap<u64, u64> = HashMap::new();
    let mut sequences_duplicated: HashSet<u64> = HashSet::new();
    let mut record_count = 0;

    each_zipped(file_handles.remove(0), file_handles.remove(0), |r1, r2| {
        if r1.is_some() & r2.is_some() {
            record_count += 1;

            match record_count % 100000 {
                0 => {
                    info!("Read and processed {} records", record_count)
                }
                _ => {}
            }

            let rec1 = r1.unwrap();
            let rec2 = r2.unwrap();

            let sequences = [rec1.seq(), rec2.seq()].concat();
            let sequences_hashed = hash64_with_seed(&sequences, 42);

            match sequences_unique.insert(sequences_hashed, record_count) {
                None => {}
                Some(is_duplicated) => {
                    sequences_duplicated.insert(is_duplicated);
                }
            }

            (true, true)
        } else {
            (false, false)
        }
    })?;

    Ok((label.to_string(), sequences_unique, sequences_duplicated))
}

fn identify_duplicate_sequences(
    sequences: HashMap<String, HashMap<u64, u64>>,
    shuffle: bool,
) -> HashMap<String, HashSet<u64>> {
    let mut sequences_unique = HashSet::new();
    let mut records_duplicated: HashMap<String, HashSet<u64>> =
        HashMap::with_capacity(sequences.len());

    //Random permutations to stop the same sequences being removed
    let keys: Vec<_> = match shuffle {
        false => sequences.keys().collect(),
        true => {
            let mut rng = thread_rng();
            let mut keys: Vec<_> = sequences.keys().collect();
            keys.shuffle(&mut rng);
            keys
        }
    };

    for key in keys {
        let mut duplicated = HashSet::new();

        for (sequence_hashed, record_number) in sequences.get(key).unwrap() {
            match sequences_unique.insert(sequence_hashed) {
                true => {}
                false => {
                    duplicated.insert(record_number.clone());
                }
            }
        }
        records_duplicated.insert(key.to_string(), duplicated);
    }

    records_duplicated
}

fn write_records<'a, R, W>(records: Vec<R>, handles: &mut Vec<W>) -> (bool, bool)
where
    R: Record,
    W: std::io::Write,
{
    let writing_results = records
        .into_iter()
        .zip(handles)
        .map(|(rec, handle)| rec.write(handle))
        .all(|r| r.is_ok());

    match writing_results {
        true => (true, true),
        false => (false, false),
    }
}

fn remove_duplicate_sequences<P>(
    paths_in: &Vec<P>,
    paths_out: &Vec<P>,
    duplicate_indexes: &HashSet<u64>,
    compression: niffler::Format,
) -> Result<FastqReadDeduplicationStats, std::io::Error>
where
    P: AsRef<Path> + Clone,
{
    let mut n_records_processed = 0;
    let mut fq_reader_handles = get_fastq_reader_file_handles(paths_in.to_vec())?;
    let mut fq_writer_handles =
        get_fastq_writer_file_handles(paths_out.to_vec(), compression, niffler::Level::Five)?;

    each_zipped(
        fq_reader_handles.remove(0),
        fq_reader_handles.remove(0),
        |r1, r2| {
            match n_records_processed % 100000 {
                0 => {
                    info!("Written {} records", n_records_processed)
                }
                _ => {}
            }
            let continue_reading = match (r1, r2) {
                (Some(rec1), Some(rec2)) => match n_records_processed {
                    n if duplicate_indexes.contains(&n) => (true, true),
                    _ => write_records(vec![rec1, rec2], &mut fq_writer_handles),
                },
                _ => (false, false),
            };

            n_records_processed += 1;
            continue_reading
        },
    )?;

    Ok(FastqReadDeduplicationStats {
        read_pairs_total: n_records_processed,
        read_pairs_duplicated: duplicate_indexes.len() as u64,
        read_pairs_unique: n_records_processed - duplicate_indexes.len() as u64,
    })
}

pub fn deduplicate_fastq(
    fq_in: HashMap<String, Vec<String>>,
    fq_out: HashMap<String, Vec<String>>,
    shuffle: bool,
    compress_output: bool,
) -> Result<FastqReadDeduplicationStats, std::io::Error> {
    
    let output_compression = match compress_output {
        true => niffler::Format::Gzip,
        false => niffler::Format::No,
    };

    let fq_hashmaps: Vec<_> = fq_in
        .clone()
        .into_par_iter()
        .map(|(label, fq_pair)| convert_fastq_to_hashmap(fq_pair, label).unwrap())
        .collect();

    let (shard_unique_seq, shard_dup_index) = fq_hashmaps.into_iter().fold(
        (HashMap::new(), HashMap::new()),
        |(mut shard_sequences, mut shard_duplicates), (label, seq, dup)| {
            shard_sequences.insert(label.clone(), seq);
            shard_duplicates.insert(label, dup);
            (shard_sequences, shard_duplicates)
        },
    );

    let collated_duplicated_indexes = identify_duplicate_sequences(shard_unique_seq, shuffle);

    let duplication_results = fq_in
        .par_iter()
        .filter_map(|(label, fq_files)| {
            let mut duplicated_indexes_by_shard =
                collated_duplicated_indexes.get(label).unwrap().to_owned();
            duplicated_indexes_by_shard.extend(shard_dup_index.get(label).unwrap().iter());

            match remove_duplicate_sequences(
                &fq_files,
                &fq_out.get(label).unwrap(),
                &duplicated_indexes_by_shard,
                output_compression,
            ) {
                Ok(res) => Some(res),
                _ => None,
            }
        })
        .reduce(
            || FastqReadDeduplicationStats {
                read_pairs_total: 0,
                read_pairs_duplicated: 0,
                read_pairs_unique: 0,
            },
            |a, b| a + b,
        );

    Ok(duplication_results.to_owned())
}

#[test]
fn test_fastq_to_hashmap() {
    let f1 = "tests/fastq_deduplicate/duplicated_1.fastq.gz";
    let f2 = "tests/fastq_deduplicate/duplicated_2.fastq.gz";

    let (_label, seq, idx) = convert_fastq_to_hashmap(vec![f1, f2], "File1").unwrap();
    assert_eq!(982, seq.len());
    assert_eq!(538, idx.len());
}

#[test]
fn test_identify() {
    let f1 = "tests/fastq_deduplicate/duplicated_1.fastq.gz";
    let r1 = "tests/fastq_deduplicate/duplicated_2.fastq.gz";

    // let f2 = "tests/fastq_deduplicate/duplicated_1.fastq.gz";
    // let r2 = "tests/fastq_deduplicate/duplicated_2.fastq.gz";

    let pair1 = convert_fastq_to_hashmap(vec![f1, r1], "1").unwrap();
    let (label, pair1_seq, _pair1_idx) = pair1;
    let mut pair1_seq_map = HashMap::new();
    pair1_seq_map.insert(label, pair1_seq);

    // println!("{:?}", &pair1_seq_map);

    let duplicates = identify_duplicate_sequences(pair1_seq_map, false);
    assert_eq!(0, duplicates.get("1").unwrap().len())
}

#[test]
fn test_remove() {
    let f1 = "tests/fastq_deduplicate/duplicated_1.fastq.gz";
    let r1 = "tests/fastq_deduplicate/duplicated_2.fastq.gz";

    let tmpdir_test = tempdir().unwrap();

    let dd1 = tmpdir_test
        .path()
        .join("out_dd_1.fq.gz")
        .to_str()
        .unwrap()
        .to_owned();
    let dd2 = tmpdir_test
        .path()
        .join("out_dd_2.fq.gz")
        .to_str()
        .unwrap()
        .to_owned();

    let pair1 = convert_fastq_to_hashmap(vec![f1, r1], "1").unwrap();
    let (label, pair1_seq, pair1_idx) = pair1;
    let mut pair1_seq_map = HashMap::new();
    pair1_seq_map.insert(label, pair1_seq);

    let mut duplicates = identify_duplicate_sequences(pair1_seq_map, false);
    duplicates.get_mut("1").unwrap().extend(pair1_idx.iter());

    let deduplication_results = remove_duplicate_sequences(
        &vec![f1, r1],
        &vec![&dd1, &dd2],
        duplicates.get("1").unwrap(),
        niffler::Format::No,
    )
    .unwrap();

    // There is an off by one error somewhere that I need to fix
    assert_eq!(983, deduplication_results.read_pairs_unique);
    assert_eq!(538, deduplication_results.read_pairs_duplicated);
}

#[test]
fn test_full_dedup() {
    let f1 = "tests/fastq_deduplicate/duplicated_1.fastq.gz".to_string();
    let r1 = "tests/fastq_deduplicate/duplicated_2.fastq.gz".to_string();

    let tmpdir_test = tempdir().unwrap();

    let dd1 = tmpdir_test
        .path()
        .join("out_dd_1.fq.gz")
        .to_str()
        .unwrap()
        .to_owned();
    let dd2 = tmpdir_test
        .path()
        .join("out_dd_2.fq.gz")
        .to_str()
        .unwrap()
        .to_owned();
    let dd11 = tmpdir_test
        .path()
        .join("out_dd1_1.fq.gz")
        .to_str()
        .unwrap()
        .to_owned();
    let dd12 = tmpdir_test
        .path()
        .join("out_dd1_2.fq.gz")
        .to_str()
        .unwrap()
        .to_owned();

    let mut fq_in = HashMap::new();
    let mut fq_out = HashMap::new();

    fq_in.insert("1".to_owned(), vec![f1.clone(), r1.clone()]);
    fq_in.insert("2".to_owned(), vec![f1, r1]);

    fq_out.insert("1".to_owned(), vec![dd1, dd2]);
    fq_out.insert("2".to_owned(), vec![dd11, dd12]);

    let res = deduplicate_fastq(fq_in, fq_out, false, false).unwrap();

    println!("{:?}", res);
}
