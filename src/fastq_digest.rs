use anyhow::Ok;
use bio::io;
use bio::pattern_matching::bom::BOM;
use crossbeam::channel;
use indicatif::{ProgressBar, ProgressIterator};
use log::{debug, info, warn};
use polars::prelude::*;
use rand::prelude::*;
use rayon::prelude::*;
use std::collections::HashMap;
use std::ops::Add;
use std::{hash::Hash, thread};
use strum::{Display, EnumString};

use crate::utils::{get_fastq_writer_file_handles, get_file_handles};

#[derive(Debug, Clone, EnumString, Display, PartialEq)]
pub enum ReadType {
    Flashed,
    Pe,
}

#[derive(Debug, Clone, EnumString, Display)]
pub enum ReadNumber {
    One,
    Two,
}

#[derive(Debug, Clone)]
pub struct DigestionHistogram {
    sample: String,
    read_type: ReadType,
    read_number: ReadNumber,
    histogram: HashMap<u64, u64>,
}

impl DigestionHistogram {
    pub fn new(sample: String, read_type: ReadType, read_number: ReadNumber) -> Self {
        Self {
            sample,
            read_type,
            read_number,
            histogram: HashMap::new(),
        }
    }

    pub fn add_entry(&mut self, slice_number: u64) {
        *self.histogram.entry(slice_number).or_insert(0) += 1;
    }

    pub fn to_dataframe(&self, description: Option<&str>) -> DataFrame {
        let sample = vec![self.sample.clone(); self.histogram.len()];
        let read_type = vec![self.read_type.to_string(); self.histogram.len()];
        let read_number = vec![self.read_number.to_string(); self.histogram.len()];

        let df = DataFrame::new(vec![
            Series::new("sample", &sample),
            Series::new("read_type", &read_type),
            Series::new("read_number", &read_number),
            Series::new(
                &description.unwrap_or("slice_number"),
                &self.histogram.keys().map(|x| *x as i64).collect::<Vec<_>>(),
            ),
            Series::new(
                "count",
                &self
                    .histogram
                    .values()
                    .map(|x| *x as i64)
                    .collect::<Vec<_>>(),
            ),
        ])
        .expect("Error creating dataframe");

        df
    }
}

impl Add for DigestionHistogram {
    type Output = Self;

    fn add(self, other: Self) -> Self {
        assert_eq!(self.read_type, other.read_type);

        let mut histogram = self.histogram.clone();
        for (k, v) in other.histogram {
            *histogram.entry(k).or_insert(0) += v;
        }
        Self {
            sample: self.sample,
            read_type: self.read_type,
            read_number: self.read_number,
            histogram,
        }
    }
}

struct DigestibleRead<'a, R> {
    read: &'a R,
    restriction_site: Vec<u8>,
    min_slice_length: usize,
    slice_number_start: usize,
    allow_undigested: bool,
    read_type: &'a ReadType,
    n_slices_unfiltered: usize,
    n_slices_filtered: usize,
}

impl<'a, R> DigestibleRead<'a, R> {
    pub fn new(
        read: &'a R,
        restriction_site: Vec<u8>,
        allow_undigested: bool,
        read_type: &'a ReadType,
        min_slice_length: Option<usize>,
        slice_number_start: Option<usize>,
    ) -> Self {
        let min_slice_length = min_slice_length.unwrap_or(18);
        let slice_number_start = slice_number_start.unwrap_or(0);

        Self {
            read,
            restriction_site,
            min_slice_length,
            slice_number_start,
            allow_undigested,
            read_type,
            n_slices_unfiltered: 0,
            n_slices_filtered: 0,
        }
    }

    fn validate_slice(&self, slice: &bio::io::fastq::Record) -> bool {
        slice.seq().len() >= self.min_slice_length
    }
}

impl DigestibleRead<'_, bio::io::fastq::Record> {
    fn digest_indicies(&self) -> Vec<usize> {
        let sequence = self.read.seq().to_ascii_lowercase().to_vec();
        let matcher = BOM::new(&self.restriction_site);
        let mut slice_indexes: Vec<_> = matcher.find_all(&sequence).collect();

        // Add the start and end to slices
        slice_indexes.insert(0, 0);
        slice_indexes.push(sequence.len());
        slice_indexes
    }

    pub fn digest(&mut self) -> Vec<bio::io::fastq::Record> {
        let mut slices = Vec::new();
        let mut slice_no = self.slice_number_start;

        let slice_indexes = self.digest_indicies();
        if slice_indexes.len() > 2 || self.allow_undigested {
            for (i, (slice_start, slice_end)) in slice_indexes
                .iter()
                .zip(slice_indexes.iter().skip(1))
                .enumerate()
            {
                let slice_start = match *slice_start {
                    0 => *slice_start,
                    _ => slice_start + self.restriction_site.len(),
                };

                let slice = self.prepare_slice(slice_start, *slice_end, slice_no);

                // Check if the slice passes the filters and add it to the vector if it does
                // If it doesn't pass the filter, increment the unfiltered counter
                match self.validate_slice(&slice) {
                    true => {
                        self.n_slices_filtered += 1;
                        self.n_slices_unfiltered += 1;
                        slices.push(slice);
                    }
                    false => self.n_slices_unfiltered += 1,
                };

                slice_no += 1;
            }
        }

        slices
    }

    fn prepare_slice(&self, start: usize, end: usize, slice_no: usize) -> bio::io::fastq::Record {
        let read_name = format!(
            "{}|{}|{}|{}",
            self.read.id(),
            self.read_type,
            slice_no,
            thread_rng().gen_range(0..100)
        );

        let sequence = self.read.seq()[start..end].to_ascii_uppercase();
        let quality = self.read.qual()[start..end].to_owned();

        bio::io::fastq::Record::with_attrs(&read_name, None, &sequence, &quality)
    }
}

#[derive(Debug, Clone)]
pub struct DigestionStats {
    sample: String,
    read_type: ReadType,
    number_of_reads: u64,
    number_of_read_pairs_unfiltered: u64,
    number_of_read_pairs_filtered: u64,
    number_of_slices_unfiltered: u64,
    number_of_slices_filtered: u64,
}

impl DigestionStats {
    pub fn new(sample: String, read_type: ReadType) -> Self {
        Self {
            sample: sample,
            read_type: read_type,
            number_of_reads: 0,
            number_of_read_pairs_unfiltered: 0,
            number_of_read_pairs_filtered: 0,
            number_of_slices_unfiltered: 0,
            number_of_slices_filtered: 0,
        }
    }
    pub fn to_dataframe(&self) -> DataFrame {
        let df = df!(
            "sample" => &[self.sample.clone()],
            "read_type" => &[self.read_type.to_string()],
            "number_of_reads" => &[self.number_of_reads],
            "number_of_read_pairs_unfiltered" => &[self.number_of_read_pairs_unfiltered],
            "number_of_read_pairs_filtered" => &[self.number_of_read_pairs_filtered],
            "number_of_slices_unfiltered" => &[self.number_of_slices_unfiltered],
            "number_of_slices_filtered" => &[self.number_of_slices_filtered]
        )
        .expect("Error creating dataframe");
        df
    }
}

pub fn digest_fastq(
    fastqs: Vec<String>,
    output: String,
    restriction_site: String,
    read_type: ReadType,
    min_slice_length: Option<usize>,
    sample: Option<String>,
) -> anyhow::Result<(
    DigestionStats,
    Vec<DigestionHistogram>,
    Vec<DigestionHistogram>,
    Vec<DigestionHistogram>,
)> {
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
    let mut stats_read_level = DigestionStats::new(sample.clone(), read_type.clone());

    let mut stats_slice_number_unfiltered = match &fastqs.len() {
        1 => vec![DigestionHistogram::new(
            sample.clone(),
            read_type.clone(),
            ReadNumber::One,
        )],
        2 => vec![
            DigestionHistogram::new(sample.clone(), read_type.clone(), ReadNumber::One),
            DigestionHistogram::new(sample.clone(), read_type.clone(), ReadNumber::Two),
        ],
        _ => panic!("Only 1 or 2 fastq files are supported"),
    };

    let mut stats_slice_number_filtered = stats_slice_number_unfiltered.clone();
    let mut stats_slice_length = stats_slice_number_unfiltered.clone();

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

                    stats_read_level.number_of_reads += 1;
                    let slices = digestible_read.digest();

                    stats_read_level.number_of_read_pairs_unfiltered += 1;
                    stats_read_level.number_of_read_pairs_filtered += match slices.len() > 0 {
                        true => 1,
                        false => 0,
                    };

                    stats_read_level.number_of_slices_unfiltered +=
                        digestible_read.n_slices_unfiltered as u64;

                    stats_read_level.number_of_slices_filtered +=
                        digestible_read.n_slices_filtered as u64;

                    stats_slice_number_unfiltered[0]
                        .add_entry(digestible_read.n_slices_unfiltered as u64);
                    stats_slice_number_filtered[0]
                        .add_entry(digestible_read.n_slices_filtered as u64);

                    for slice in slices {
                        stats_slice_length[0].add_entry(slice.seq().len() as u64);
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

                        stats_read_level.number_of_reads += 2;
                        stats_read_level.number_of_read_pairs_unfiltered += 1;

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

                        stats_read_level.number_of_read_pairs_filtered +=
                            match (slices_1.len(), slices_2.len()) {
                                (0, 0) => 0,
                                _ => 1,
                            };

                        stats_read_level.number_of_slices_unfiltered +=
                            digestible_read_1.n_slices_unfiltered as u64
                                + digestible_read_2.n_slices_unfiltered as u64;

                        stats_read_level.number_of_slices_filtered +=
                            digestible_read_1.n_slices_filtered as u64
                                + digestible_read_2.n_slices_filtered as u64;

                        stats_slice_number_unfiltered[0]
                            .add_entry(digestible_read_1.n_slices_unfiltered as u64);
                        stats_slice_number_filtered[0]
                            .add_entry(digestible_read_1.n_slices_filtered as u64);

                        stats_slice_number_unfiltered[1]
                            .add_entry(digestible_read_2.n_slices_unfiltered as u64);
                        stats_slice_number_filtered[1]
                            .add_entry(digestible_read_2.n_slices_filtered as u64);

                        for slice in slices_1.into_iter() {
                            stats_slice_length[0].add_entry(slice.seq().len() as u64);
                            slice_tx.send(slice).unwrap();
                        }

                        for slice in slices_2.into_iter() {
                            stats_slice_length[1].add_entry(slice.seq().len() as u64);
                            slice_tx.send(slice).unwrap();
                        }
                    })
            }
            _ => {}
        }
        stats_tx
            .send((
                stats_read_level,
                stats_slice_number_unfiltered,
                stats_slice_number_filtered,
                stats_slice_length,
            ))
            .expect("Error sending stats");
        drop(stats_tx);
        drop(slice_tx);
    });

    for slice in slice_rx {
        fastq_writer.write_record(&slice).unwrap();
    }

    digestion_thread.join().unwrap();

    let stats = stats_rx.recv().unwrap();

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

        let (read_stats, hist_unfilt, hist_filt, hist_len) = digest_fastq(
            vec![fq1.to_string(), fq2.to_string()],
            output.to_string(),
            "GATC".to_lowercase().to_string(),
            ReadType::Pe,
            None,
            Some("test_sample".to_string()),
        )?;

        assert!(std::path::Path::new(output).exists());
        assert_eq!(
            read_stats.number_of_read_pairs_filtered,
            read_stats.number_of_read_pairs_unfiltered
        );

        Ok(())
    }

    #[test]
    fn test_digest_fastq_flashed() -> Result<()> {
        let fq1 = "tests/fastq_digest/digest_1.fastq.gz";
        let output = "tests/fastq_digest/digest_output.fastq.gz";

        let (read_stats, hist_unfilt, hist_filt, hist_len) = digest_fastq(
            vec![fq1.to_string()],
            output.to_string(),
            "GATC".to_lowercase().to_string(),
            ReadType::Flashed,
            None,
            Some("test_sample".to_string()),
        )?;

        assert!(std::path::Path::new(output).exists());
        assert_eq!(read_stats.number_of_read_pairs_filtered, 876);

        Ok(())
    }
}
