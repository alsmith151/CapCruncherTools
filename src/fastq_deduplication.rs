use ahash::{AHasher};
use fastq::{each_zipped, Parser, Record};
use hashbrown::{HashMap, HashSet};
use niffler;
use rand::prelude::*;
use rayon::iter::{IntoParallelIterator, IntoParallelRefIterator, ParallelIterator};
use serde::{Deserialize, Serialize};
use std::hash::{Hasher};
use std::ops::Add;
use std::path::{Path};
use std::prelude::rust_2021::*;

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
) -> Result<(String, HashMap<u64, u64>), std::io::Error> {
    let mut file_handles = get_fastq_reader_file_handles(paths)?;
    let mut sequences_unique: HashMap<u64, u64> = HashMap::new();
    let mut record_count = 0;

    each_zipped(file_handles.remove(0), file_handles.remove(0), |r1, r2| {
        if r1.is_some() & r2.is_some() {
            record_count += 1;
            let rec1 = r1.unwrap();
            let rec2 = r2.unwrap();

            let sequences = [rec1.seq(), rec2.seq()].concat();
            let mut hasher = AHasher::new_with_keys(42, 69);
            hasher.write(&sequences);
            let sequences_hashed = hasher.finish();

            sequences_unique.insert(sequences_hashed, record_count);

            (true, true)
        } else {
            (false, false)
        }
    })?;

    Ok((label.to_string(), sequences_unique))
}

fn identify_duplicate_sequences(
    sequences: Vec<(String, HashMap<u64, u64>)>,
) -> HashMap<String, HashSet<u64>> {
    let mut sequences_unique = HashSet::new();
    let mut records_duplicated: HashMap<String, HashSet<u64>> =
        HashMap::with_capacity(sequences.len());

    //Random permutations to stop the same sequences being removed
    let mut sequences_to_dedup = sequences.clone();
    let mut rng = thread_rng();
    sequences_to_dedup.shuffle(&mut rng);

    for (label, sequences_hashed_map) in sequences_to_dedup {
        let mut duplicated = HashSet::new();

        for (sequence_hashed, record_number) in sequences_hashed_map {
            match sequences_unique.insert(sequence_hashed) {
                true => {}
                false => {
                    duplicated.insert(record_number);
                }
            }
        }
        records_duplicated.insert(label.to_string(), duplicated);
    }

    records_duplicated
}

fn write_records<'a, R, W>(records: Vec<R>, handles: &mut Vec<W>) -> (bool, bool)
where
    R: Record,
    W: std::io::Write,
{
    let writing_results = records
        .iter()
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
) -> Result<FastqReadDeduplicationStats, std::io::Error>
where
    P: AsRef<Path> + Clone,
{
    let mut n_records_processed = 0;
    let mut fq_reader_handles = get_fastq_reader_file_handles(paths_in.to_vec())?;
    let mut fq_writer_handles = get_fastq_writer_file_handles(
        paths_out.to_vec(),
        niffler::Format::Gzip,
        niffler::Level::Five,
    )?;

    each_zipped(
        fq_reader_handles.remove(0),
        fq_reader_handles.remove(0),
        |r1, r2| {
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
) -> Result<FastqReadDeduplicationStats, std::io::Error>

{
    let fq_hashmaps: Vec<_> = fq_in
        .clone()
        .into_par_iter()
        .map(|(label, fq_pair)| convert_fastq_to_hashmap(fq_pair, label).unwrap())
        .collect();

    let duplicated_indexes = identify_duplicate_sequences(fq_hashmaps);

    let duplication_results = fq_in
        .par_iter()
        .filter_map(|(label, fq_files)| {
            let duplicate_indexes_sample_specific = duplicated_indexes.get(label).unwrap();

            match remove_duplicate_sequences(
                &fq_files,
                &fq_out.get(label).unwrap(),
                duplicate_indexes_sample_specific,
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
    let f1 = "test_DUP_1.fq.gz";
    let f2 = "test_DUP_1.fq.gz";

    let filt = convert_fastq_to_hashmap(vec![f1, f2], "File1").unwrap();
    // println!("{:?}", &filt);
    // println!("{:?}", &filt.1.len());
}

#[test]
fn test_identify() {
    let f1 = "test_1.fastq.gz";
    let r1 = "test_1.fastq.gz";

    let f2 = "test_DUP_1.fq.gz";
    let r2 = "test_DUP_2.fq.gz";

    let pair1 = convert_fastq_to_hashmap(vec![f1, r1], "1").unwrap();
    let pair2 = convert_fastq_to_hashmap(vec![f2, r2], "2").unwrap();

    let duplicates = identify_duplicate_sequences(vec![pair1, pair2]);
    // println!("{:?}", duplicates);
}

#[test]
fn test_remove() {
    let f1 = "test_1.fastq.gz";
    let r1 = "test_1.fastq.gz";

    let f2 = "test_DUP_1.fq.gz";
    let r2 = "test_DUP_2.fq.gz";

    let dd1 = "out_dd_1.fq.gz";
    let dd2 = "out_dd_2.fq.gz";

    let pair1 = convert_fastq_to_hashmap(vec![f1, r1], "1").unwrap();
    let pair2 = convert_fastq_to_hashmap(vec![f2, r2], "2").unwrap();

    let duplicates = identify_duplicate_sequences(vec![pair1, pair2]);

    let removed =
        remove_duplicate_sequences(&vec![f1, r1], &vec![dd1, dd2], duplicates.get("1").unwrap());
    println!("{:?}", removed.unwrap());
}

#[test]
fn test_full_dedup() {

    let f1 = "test_1.fastq.gz".to_owned();
    let r1 = "test_1.fastq.gz".to_owned();

    let f2 = "test_DUP_1.fq.gz".to_owned();
    let r2 = "test_DUP_2.fq.gz".to_owned();

    let dd1 = "out_dd_1.fq.gz".to_owned();
    let dd2 = "out_dd_2.fq.gz".to_owned();

    let dd11 = "out_dd2_1.fq.gz".to_owned();
    let dd12 = "out_dd2_2.fq.gz".to_owned();


    let mut fq_in = HashMap::new();
    let mut fq_out = HashMap::new();
    
    fq_in.insert("1".to_owned(), vec![f1, r1]);
    fq_in.insert("2".to_owned(), vec![f2, r2]);

    fq_out.insert("1".to_owned(), vec![dd1, dd2]);
    fq_out.insert("2".to_owned(), vec![dd11, dd12]);

    let res = deduplicate_fastq(fq_in, fq_out).unwrap();

    println!("{:?}", res);




}