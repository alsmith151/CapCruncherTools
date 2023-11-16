use anyhow::Ok;
use bio::io;
use bio::pattern_matching::bom::BOM;
use crossbeam::channel;
use indicatif::{ProgressBar, ProgressIterator};
use log::{debug, error, info, warn};
use noodles::bam::record::cigar::Op;
use polars::prelude::*;
use rand::prelude::*;
use rayon::prelude::*;
use serde::ser::{SerializeStruct, Serializer};
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::ops::Add;
use std::{hash::Hash, thread};
use strum::{Display, EnumString};

use crate::digest::DigestibleRead;
use crate::utils::{get_fastq_writer_file_handles, get_file_handles, ReadNumber, ReadType};

#[derive(Serialize, Deserialize, Debug, Clone)]
struct SliceNumberStats {
    unfiltered: u64,
    filtered: u64,
}

#[derive(Deserialize, Debug, Clone)]
struct ReadPairStat<V> {
    read1: V,
    read2: Option<V>,
}

impl Serialize for ReadPairStat<u64> {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        let read2 = match &self.read2 {
            Some(v) => v,
            None => &0,
        };

        let mut s = serializer.serialize_struct("ReadPairStat", 2)?;
        s.serialize_field("read1", &self.read1)?;
        s.serialize_field("read2", &read2)?;
        s.end()
    }
}

impl Serialize for ReadPairStat<Histogram> {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        let read2 = match &self.read2 {
            Some(v) => v.to_owned(),
            None => Histogram::new("slice_number".to_string()),
        };

        let mut s = serializer.serialize_struct("ReadPairStat", 2)?;
        s.serialize_field("read1", &self.read1)?;
        s.serialize_field("read2", &read2)?;
        s.end()
    }
}

#[derive(Serialize, Deserialize, Debug, Clone)]
struct DigestionReadPairStats {
    unfiltered: ReadPairStat<u64>,
    filtered: ReadPairStat<u64>,
}

#[derive(Serialize, Deserialize, Debug, Clone)]
struct Histogram {
    name: String,
    hist: HashMap<u64, u64>,
}

impl Histogram {
    pub fn new(name: String) -> Self {
        Self {
            name,
            hist: HashMap::new(),
        }
    }

    pub fn add_entry(&mut self, entry: u64) {
        let count = self.hist.entry(entry).or_insert(0);
        *count += 1;
    }
}

#[derive(Serialize, Deserialize, Debug, Clone)]
struct DigestionHistograms {
    unfiltered: ReadPairStat<Histogram>,
    filtered: ReadPairStat<Histogram>,
    lengths: ReadPairStat<Histogram>,
}

#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct DigestionStats {
    sample: String,
    read_type: ReadType,
    read_stats: DigestionReadPairStats,
    slice_stats: SliceNumberStats,
    histograms: DigestionHistograms,
}

impl DigestionStats {
    pub fn new(sample: String, read_type: ReadType) -> Self {
        let is_paired = match read_type {
            ReadType::Flashed => false,
            ReadType::Pe => true,
        };

        Self {
            sample,
            read_type,
            read_stats: DigestionReadPairStats {
                unfiltered: ReadPairStat {
                    read1: 0,
                    read2: is_paired.then(|| 0),
                },
                filtered: ReadPairStat {
                    read1: 0,
                    read2: is_paired.then(|| 0),
                },
            },
            slice_stats: SliceNumberStats {
                unfiltered: 0,
                filtered: 0,
            },
            histograms: DigestionHistograms {
                unfiltered: ReadPairStat {
                    read1: Histogram::new("slice_number".to_string()),
                    read2: is_paired.then(|| Histogram::new("slice_number".to_string())),
                },
                filtered: ReadPairStat {
                    read1: Histogram::new("slice_number".to_string()),
                    read2: is_paired.then(|| Histogram::new("slice_number".to_string())),
                },
                lengths: ReadPairStat {
                    read1: Histogram::new("slice_length".to_string()),
                    read2: is_paired.then(|| Histogram::new("slice_length".to_string())),
                },
            },
        }
    }
}

#[derive(Debug, Clone, EnumString, Display, PartialEq, Serialize, Deserialize)]
pub enum ValidDigestionParams {
    Valid,
    Invalid,
    InvalidNoEnzyme,
}

impl ValidDigestionParams {
    pub fn validate(n_files: usize, read_type: ReadType) -> Self {
        match (n_files, read_type) {
            (1, ReadType::Flashed) => ValidDigestionParams::Valid,
            (2, ReadType::Pe) => ValidDigestionParams::Valid,
            _ => {
                error!("Invalid combination of input files and read type");
                ValidDigestionParams::Invalid
            }
        }
    }
}

pub fn digest_fastq(
    fastqs: Vec<String>,
    output: String,
    restriction_site: String,
    read_type: ReadType,
    min_slice_length: Option<usize>,
    sample: Option<String>,
) -> anyhow::Result<DigestionStats> {
    let (slice_tx, slice_rx) = channel::unbounded();
    let (stats_tx, stats_rx) = channel::unbounded();

    let mut handles_writer = get_fastq_writer_file_handles(
        vec![output.clone()],
        niffler::Format::Gzip,
        Some(niffler::Level::One),
    )?;
    let mut fastq_writer = bio::io::fastq::Writer::new(handles_writer[0].as_mut());

    let restriction_site = restriction_site.as_bytes().to_vec();
    let allow_undigested = match read_type {
        ReadType::Flashed => false,
        ReadType::Pe => true,
    };

    let sample = sample.unwrap_or("digested.fastq.gz".to_string());
    let mut digestion_stats = DigestionStats::new(sample, read_type.clone());

    let digestion_thread = std::thread::spawn(move || {
        let mut handles_reader =
            get_file_handles(fastqs.clone()).expect("Error getting file handles");
        match fastqs.len() {
            1 => {
                let reader = io::fastq::Reader::new(handles_reader[0].as_mut());

                reader.records().for_each(|res| {
                    let r = res.expect("Error reading record");
                    let mut digestible_read = DigestibleRead::new(
                        &r,
                        restriction_site.clone(),
                        allow_undigested,
                        &read_type,
                        min_slice_length,
                        None,
                    );

                    // Update number of reads
                    digestion_stats.read_stats.unfiltered.read1 += 1;

                    // Digest the read
                    let slices = digestible_read.digest();

                    // Update stats
                    digestion_stats.read_stats.filtered.read1 += match slices.len() > 0 {
                        true => 1,
                        false => 0,
                    };

                    // Update slice stats
                    digestion_stats.slice_stats.unfiltered +=
                        digestible_read.n_slices_unfiltered as u64;

                    digestion_stats.slice_stats.filtered +=
                        digestible_read.n_slices_filtered as u64;

                    // Update histograms

                    digestion_stats
                        .histograms
                        .unfiltered
                        .read1
                        .add_entry(digestible_read.n_slices_unfiltered as u64);
                    digestion_stats
                        .histograms
                        .filtered
                        .read1
                        .add_entry(digestible_read.n_slices_filtered as u64);

                    for slice in slices {
                        digestion_stats
                            .histograms
                            .lengths
                            .read1
                            .add_entry(slice.seq().len() as u64);
                        slice_tx.send(slice).unwrap();
                    }
                })
            }
            2 => {
                let reader_1 = io::fastq::Reader::new(handles_reader.pop().unwrap());
                let reader_2 = io::fastq::Reader::new(handles_reader.pop().unwrap());

                reader_1
                    .records()
                    .zip(reader_2.records())
                    .for_each(|(res1, res2)| {
                        let r1 = res1.expect("Error getting read 1");
                        let r2 = res2.expect("Error getting read 2");

                        // Update number of reads
                        digestion_stats.read_stats.unfiltered.read1 += 1;
                        digestion_stats.read_stats.unfiltered.read2 = Some(
                            digestion_stats
                                .read_stats
                                .unfiltered
                                .read2
                                .expect("Read 2 not present")
                                + 1,
                        );

                        let mut digestible_read_1 = DigestibleRead::new(
                            &r1,
                            restriction_site.clone(),
                            allow_undigested,
                            &read_type,
                            min_slice_length,
                            None,
                        );

                        let slices_1 = digestible_read_1.digest();
                        let mut digestible_read_2 = DigestibleRead::new(
                            &r2,
                            restriction_site.clone(),
                            allow_undigested,
                            &read_type,
                            min_slice_length,
                            Some(digestible_read_1.n_slices_filtered),
                        );

                        let slices_2 = digestible_read_2.digest();

                        digestion_stats.read_stats.filtered.read1 +=
                            match (slices_1.len(), slices_2.len()) {
                                (0, 0) => 0,
                                _ => 1,
                            };

                        digestion_stats.read_stats.filtered.read2 =
                            Some(match (slices_1.len(), slices_2.len()) {
                                (0, 0) => digestion_stats.read_stats.filtered.read2.unwrap(),
                                _ => digestion_stats.read_stats.filtered.read2.unwrap() + 1,
                            });

                        // Update slice stats
                        digestion_stats.slice_stats.unfiltered +=
                            digestible_read_1.n_slices_unfiltered as u64
                                + digestible_read_2.n_slices_unfiltered as u64;

                        digestion_stats.slice_stats.filtered += digestible_read_1.n_slices_filtered
                            as u64
                            + digestible_read_2.n_slices_filtered as u64;

                        // Update histograms

                        // Unfiltered
                        digestion_stats
                            .histograms
                            .unfiltered
                            .read1
                            .add_entry(digestible_read_1.n_slices_unfiltered as u64);

                        digestion_stats
                            .histograms
                            .unfiltered
                            .read2
                            .as_mut()
                            .unwrap()
                            .add_entry(digestible_read_2.n_slices_unfiltered as u64);

                        // Filtered
                        digestion_stats
                            .histograms
                            .filtered
                            .read1
                            .add_entry(digestible_read_1.n_slices_filtered as u64);

                        digestion_stats
                            .histograms
                            .filtered
                            .read2
                            .as_mut()
                            .unwrap()
                            .add_entry(digestible_read_2.n_slices_filtered as u64);

                        for slice in slices_1.into_iter() {
                            digestion_stats
                                .histograms
                                .lengths
                                .read1
                                .add_entry(slice.seq().len() as u64);
                            slice_tx.send(slice).unwrap();
                        }

                        for slice in slices_2.into_iter() {
                            digestion_stats
                                .histograms
                                .lengths
                                .read2
                                .as_mut()
                                .unwrap()
                                .add_entry(slice.seq().len() as u64);
                            slice_tx.send(slice).unwrap();
                        }
                    })
            }
            _ => {}
        }
        stats_tx
            .send((digestion_stats,))
            .expect("Error sending stats");
        drop(stats_tx);
        drop(slice_tx);
    });

    for slice in slice_rx {
        fastq_writer.write_record(&slice).unwrap();
    }

    digestion_thread.join().expect("Error finalizing thread");

    let stats = stats_rx.recv().expect("Error receiving stats").0;

    Ok(stats)
}

// Tests
#[cfg(test)]
mod tests {
    use super::*;
    use anyhow::Result;
    use bio::io::fastq::Record;

    #[test]
    fn test_digestible_read() -> Result<()> {
        let read = Record::with_attrs(
            "@test_read",
            None,
            b"AAAAAAAAAAAAAAAAAAAAAAAAGATCAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA",
            b"############################################################",
        );
        let mut digestible_read = DigestibleRead::new(
            &read,
            b"GATC".to_ascii_lowercase().to_vec(),
            false,
            &ReadType::Flashed,
            None,
            None,
        );
        let slices = digestible_read.digest();

        println!("{:?}", slices);

        assert_eq!(slices.len(), 2);

        Ok(())
    }

    #[test]
    fn test_digest_fastq_pe() -> Result<()> {
        let fq1 = "tests/fastq_digest/digest_1.fastq.gz";
        let fq2 = "tests/fastq_digest/digest_2.fastq.gz";
        let output = "tests/fastq_digest/digest_output.fastq.gz";

        let (stats) = digest_fastq(
            vec![fq1.to_string(), fq2.to_string()],
            output.to_string(),
            "GATC".to_lowercase().to_string(),
            ReadType::Pe,
            None,
            Some("test_sample".to_string()),
        )?;

        println!("{:?}", stats);

        assert!(std::path::Path::new(output).exists());
        // assert_eq!(
        //     read_stats.number_of_read_pairs_filtered,
        //     read_stats.number_of_read_pairs_unfiltered
        // );

        Ok(())
    }

    #[test]
    fn test_digest_fastq_flashed() -> Result<()> {
        let fq1 = "tests/fastq_digest/digest_1.fastq.gz";
        let output = "tests/fastq_digest/digest_output.fastq.gz";

        let (stats) = digest_fastq(
            vec![fq1.to_string()],
            output.to_string(),
            "GATC".to_lowercase().to_string(),
            ReadType::Flashed,
            None,
            Some("test_sample".to_string()),
        )?;

        assert!(std::path::Path::new(output).exists());
        // assert_eq!(read_stats.number_of_read_pairs_filtered, 876);

        Ok(())
    }
}
