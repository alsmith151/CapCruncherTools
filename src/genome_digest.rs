use std::fs::File;
use std::sync::Arc;

use bio::io::bed;
use bio::io::fasta;
use bio::pattern_matching::bom::BOM;
use log::{debug, error, info, warn};
use memchr::memmem;
use rayon::prelude::*;
use std::prelude::rust_2021::*;
use std::rc::Rc;

struct DigestedFastaEntry<'a> {
    name: String,
    sequence: &'a [u8],
    restriction_site: &'a [u8],
    remove_recognition_site: bool,
    min_slice_size: Option<usize>,
    slices: Vec<usize>,
}

impl<'a> DigestedFastaEntry<'a> {
    fn new(
        name: String,
        sequence: &'a &[u8],
        restriction_site: &'a [u8],
        remove_recognition_site: bool,
        min_slice_size: Option<usize>,
    ) -> Self {
        Self {
            name,
            sequence,
            restriction_site,
            remove_recognition_site,
            min_slice_size,
            slices: Vec::new(),
        }
    }

    fn digest(&mut self) {
        let matcher = BOM::new(self.restriction_site);
        let mut slice_indexes: Vec<_> = matcher.find_all(self.sequence).collect();

        // Add the start and end to slices
        slice_indexes.insert(0, 0);
        slice_indexes.push(self.sequence.len());

        // Remove all slices that are smaller than the minimum slice size
        if let Some(min_slice_size) = self.min_slice_size {
            slice_indexes = slice_indexes
                .into_iter()
                .filter(|&x| x >= min_slice_size)
                .collect();
        }

        self.slices = slice_indexes;
    }

    fn to_bed_records(&mut self) -> Vec<bed::Record> {
        self.digest();
        let mut bed_records = Vec::with_capacity(self.slices.len());

        for (start, end) in self.slices.iter().zip(self.slices.iter().skip(1)) {
            let mut bed_rec = bed::Record::new();
            bed_rec.set_chrom(&self.name);

            // Remove the recognition site from the slice if specified
            if self.remove_recognition_site {
                bed_rec.set_start((*start + self.restriction_site.len()) as u64);
            } else {
                bed_rec.set_start(*start as u64);
            }

            bed_rec.set_end(*end as u64);
            bed_records.push(bed_rec);
        }
        bed_records
    }
}

pub fn digest_fasta(
    fasta: String,
    restriction_site: String,
    output: String,
    remove_recognition_site: bool,
    min_slice_length: Option<usize>,
) -> Result<(), std::io::Error> {
    let mut reader = fasta::Reader::from_file(fasta).expect("Error opening FASTA file");
    let mut writer = bed::Writer::to_file(output).expect("Error opening BED output file");

    let min_slice_length = min_slice_length.unwrap_or(0);
    let restriction_site = restriction_site.as_bytes();

    // Iterate over the records in parallel with rayon and digest each one
    reader
        .records()
        .par_bridge()
        .map(|record| {
            let record = record.unwrap();
            let seq = record.seq();

            let mut digested_entry = DigestedFastaEntry::new(
                record.id().to_string(),
                &seq,
                restriction_site,
                remove_recognition_site,
                Some(min_slice_length),
            );

            digested_entry.to_bed_records()
        })
        .collect::<Vec<Vec<bed::Record>>>()
        .iter()
        .for_each(|bed_records| {
            for bed_rec in bed_records {
                writer.write(bed_rec).unwrap();
            }
        });

    Ok(())
}

// Test fasa digestion
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_digest_fasta() {
        let fasta = "tests/chr14.fa".to_string();
        let restriction_site = "GATC".to_string();
        let output = "test.bed".to_string();
        let remove_recognition_site = true;

        digest_fasta(
            fasta,
            restriction_site,
            output,
            remove_recognition_site,
            None,
        )
        .unwrap();
    }
}

// pub fn digest_fasta_file(fasta_file: String, recognition_site: &[u8]) -> Result<(), std::io::Error>{

//     let mut reader = fasta::Reader::from_file(fasta_file).expect("Error opening file");
//     let mut writer = fasta::Writer::to_file("test.bed")?;

//     let mut bed_records = Vec::new();

//     // Iterate over the records in parallel with rayon and digest each one
//     reader.records().par_bridge().for_each(|record| {
//         let record = record.unwrap();
//         let sequence = record.seq();
//         let matcher = BOM::new(recognition_site);
//         let mut slice_indexes: Vec<_> = matcher.find_all(sequence).collect();

//         // Add the start and end to slices
//         slice_indexes.insert(0, 0);
//         slice_indexes.push(sequence.len());

//         // Create a bed record for each slice
//         for (i, slice) in slice_indexes.iter().enumerate(){
//             let bed_rec = bed::Record::new();
//             bed_records.push(bed_rec);
//         }
//     });

//     Ok(())
// }

// fn digest_entry(entry: fasta::Record, recognition_site: &[u8]) -> Result<Vec<usize>, std::io::Error>{

//     let sequence = entry.seq();
//     let matcher = BOM::new(recognition_site);
//     let mut slice_indexes: Vec<_> = matcher.find_all(sequence).collect();

//     // Add the start and end to slices
//     slice_indexes.insert(0, 0);
//     slice_indexes.push(sequence.len());

//     Ok(slice_indexes)

// }

// fn convert_fasta_to_digested_bed(fasta_file: String, recognition_site: &[u8]) -> Result<(), std::io::Error>{

//     let mut reader = fasta::Reader::from_file(fasta_file).expect("Error opening file");
//     let mut writer = fasta::Writer::to_file("test.bed")?;

//     let bed_records = Vec::new();

//     // Iterate over the records in parallel with rayon and digest each one
//     reader.records().par_bridge().for_each(|record| {
//         let record = record.unwrap();
//         let slice_indexes = digest_entry(record, recognition_site).unwrap();
//         let bed_rec = bed::Record::new();
//     });

//     Ok(())
// }

// #[test]
// fn test_digest_entry_with_rsite(){

//     let sequence_parts = ["N".repeat(10), "GATC".to_string(),"N".repeat(10)];
//     let sequence_full = sequence_parts.join("");
//     let record = fasta::Record::with_attrs("record_1", None, sequence_full.as_bytes());
//     let res = digest_entry(record, b"GATC");
//     println!("{:?}", res)

// }

// #[test]
// fn test_digest_fasta_file(){

//         let fasta_file = "tests/chr14.fa";
//         let recognition_site = b"GATC";
//         convert_fasta_to_digested_bed(fasta_file.to_string(), recognition_site).unwrap();

// }
