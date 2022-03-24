use std::prelude::rust_2021::*;
use log::{debug, info, warn};
use bio::io::fasta;
use bio::pattern_matching::bom::BOM;


fn digest_entry(entry: fasta::Record, recognition_site: &[u8]) -> Result<Vec<usize>, std::io::Error>{

    let sequence = entry.seq();
    let matcher = BOM::new(recognition_site);
    let mut slice_indexes: Vec<_> = matcher.find_all(sequence).collect();

    // Add the start and end to slices
    slice_indexes.insert(0, 0);
    slice_indexes.push(sequence.len());

    Ok(slice_indexes)

}





#[test]
fn test_digest_entry(){

    let sequence_parts = ["N".repeat(10), "GATC".to_string(),"N".repeat(10)];
    let sequence_full = sequence_parts.join("");
    let record = fasta::Record::with_attrs("record_1", None, sequence_full.as_bytes());
    let res = digest_entry(record, b"GATC");
    println!("{:?}", res)

}