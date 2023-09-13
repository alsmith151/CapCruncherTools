use fastq::{each_zipped, Parser, Record};
use log::{debug, info, warn};
use std::iter::Iterator;
use std::path::Path;
use std::prelude::rust_2021::*;



pub fn get_fastq_reader_file_handles<P>(
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

pub fn get_file_handles<P>(
    paths: Vec<P>,
) -> Result<Vec<Box<dyn std::io::Read>>, std::io::Error>
where
    P: AsRef<Path>,
{
    let file_handles: Vec<_> = paths
        .iter()
        .map(|p| niffler::from_path(p.as_ref()).expect("Compression not recognised"))
        .map(|(fh, _fmt)| fh)
        .collect();

    Ok(file_handles)
}


pub fn get_fastq_writer_file_handles<P>(
    paths: Vec<P>,
    compression_format: niffler::compression::Format,
    compression_level: Option<niffler::Level>,
) -> Result<Vec<Box<dyn std::io::Write>>, std::io::Error>
where
    P: AsRef<Path>,
{

    let compression_level = match compression_level {
        Some(l) => l,
        None => niffler::Level::Five,
    };

    let file_handles: Vec<_> = paths
        .iter()
        .map(|p| {
            niffler::to_path(p.as_ref(), compression_format, compression_level)
                .expect("Failed to create output file")
        })
        .collect();

    Ok(file_handles)
}

pub fn write_records<R, W>(records: Vec<R>, handles: &mut Vec<W>) -> (bool, bool)
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
