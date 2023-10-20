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
use serde::{Serialize, Deserialize};


use crate::utils::{ReadType};

pub struct DigestibleRead<'a, R> {
    read: &'a R,
    restriction_site: Vec<u8>,
    min_slice_length: usize,
    slice_number_start: usize,
    allow_undigested: bool,
    read_type: &'a ReadType,
    pub n_slices_unfiltered: usize,
    pub n_slices_filtered: usize,
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
