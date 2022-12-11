// use std::prelude::rust_2021::*;
// use log::{debug, info, warn};
// use bio::io::fasta;
// use bio::io::bed;
// use bio::pattern_matching::bom::BOM;
// use rayon::prelude::*;


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