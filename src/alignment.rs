// use noodles::bam;
// use log::{debug, info, warn};
// use std::iter::Iterator;
// use std::path::Path;
// use std::collections::HashMap;
// use std::prelude::rust_2021::*;




// // CCAlignment = namedtuple(
// //     "CCAlignment",
// //     field_names=[
// //         "slice_id",
// //         "slice_name",
// //         "parent_id",
// //         "parent_read",
// //         "pe",
// //         "slice",
// //         "uid",
// //         "mapped",
// //         "multimapped",
// //         "chrom",
// //         "start",
// //         "end",
// //         "coordinates",
// //     ],
// // )


// pub struct CCAlignment {
//     slice_id: u64,
//     slice_name: String,
//     parent_id: u64,
//     parent_read: String,
//     pe: String,
//     slice: u64,
//     uid: u64,
//     mapped: bool,
//     multimapped: bool,
//     chrom: String,
//     start: u64,
//     end: u64,
//     coordinates: String,
// }

// impl CCAlignment {
//     pub fn new(
//         slice_id: u64,
//         slice_name: String,
//         parent_id: u64,
//         parent_read: String,
//         pe: String,
//         slice: u64,
//         uid: u64,
//         mapped: bool,
//         multimapped: bool,
//         chrom: String,
//         start: u64,
//         end: u64,
//         coordinates: String,
//     ) -> Self {
//         Self {
//             slice_id,
//             slice_name,
//             parent_id,
//             parent_read,
//             pe,
//             slice,
//             uid,
//             mapped,
//             multimapped,
//             chrom,
//             start,
//             end,
//             coordinates,
//         }
//     }

//     pub fn values(&self) -> Vec<String> {
//         vec![
//             self.slice_id.to_string(),
//             self.slice_name.clone(),
//             self.parent_id.to_string(),
//             self.parent_read.clone(),
//             self.pe.clone(),
//             self.slice.to_string(),
//             self.uid.to_string(),
//             self.mapped.to_string(),
//             self.multimapped.to_string(),
//             self.chrom.clone(),
//             self.start.to_string(),
//             self.end.to_string(),
//             self.coordinates.clone(),
//         ]
//     }

//     pub fn keys(&self) -> Vec<String> {
//         vec![
//             "slice_id".to_string(),
//             "slice_name".to_string(),
//             "parent_id".to_string(),
//             "parent_read".to_string(),
//             "pe".to_string(),
//             "slice".to_string(),
//             "uid".to_string(),
//             "mapped".to_string(),
//             "multimapped".to_string(),
//             "chrom".to_string(),
//             "start".to_string(),
//             "end".to_string(),
//             "coordinates".to_string(),
//         ]
//     }

// }

// impl From <bam::Record> for CCAlignment {
//     fn from(record: bam::Record) -> Self {
//         let slice_name = String::from_utf8_lossy(record.read_name()).into_owned();

//         let n = &slice_name.split("|").collect::<Vec<&str>>();
//         let parent_read = n[0].to_string();
//         let pe = n[1].to_string();
//         let slice_number = n[2].into();

//         let slice_id = xxhash::xxh3::xxh3_64(&slice_name.as_bytes());
//         let flags = record.flags();

//         let aln = match flags.is_unmapped() {
//             true => {
//                 let mapped = false;
//                 let multimapped = false;
//                 let chrom = "".to_string();
//                 let start = 0;
//                 let end = 0;
//                 let coordinates = "".to_string();

//                 Self::new(
//                     slice_id,
//                     slice_name,
//                     0,
//                     parent_read,
//                     pe,
//                     slice_number,
//                     0,
//                     mapped,
//                     multimapped,
//                     chrom,
//                     start,
//                     end,
//                     coordinates,
//                 )
//             },
//             false => {
//                 let mapped = true;
//                 let multimapped = flags.is_secondary();
//                 let chrom = record.reference_sequence_id().unwrap();
//                 let start = record.pos() as u64;
//                 let end = record.cigar().end_pos() as u64;
//                 let coordinates = format!("{}:{}-{}", chrom, start, end);

//                 Self::new(
//                     slice_id,
//                     slice_name,
//                     0,
//                     parent_read,
//                     pe,
//                     slice_number,
//                     0,
//                     mapped,
//                     multimapped,
//                     chrom,
//                     start,
//                     end,
//                     coordinates,
//                 )
//             }
//         };
        
        




//         let mapped = record.is_unmapped();
//         let multimapped = record.is_secondary();
//         let chrom = record.tid().to_string();
//         let start = record.pos() as u64;
//         let end = record.cigar().end_pos() as u64;
//         let coordinates = format!("{}:{}-{}", chrom, start, end);

//         Self::new(
//             slice_id,
//             slice_name,
//             parent_read,
//             pe,
//             slice_number,
//             uid,
//             mapped,
//             multimapped,
//             chrom,
//             start,
//             end,
//             coordinates,
//         )
//     }
// }


// // def parse_alignment(aln: pysam.AlignmentFile) -> CCAlignment:
// //     """Parses reads from a bam file into a list.

// //     Extracts:
// //      -read name
// //      -parent reads
// //      -flashed status
// //      -slice number
// //      -mapped status
// //      -multimapping status
// //      -chromosome number (e.g. chr10)
// //      -start (e.g. 1000)
// //      -end (e.g. 2000)
// //      -coords e.g. (chr10:1000-2000)


// //     Args:
// //      aln: pysam.AlignmentFile.
// //     Returns:
// //      list: Containing the attributes extracted.

// //     """

// //     slice_name = aln.query_name
// //     parent_read, pe, slice_number, uid = slice_name.split("|")
// //     parent_id = xxhash.xxh3_64_intdigest(parent_read, seed=42)
// //     slice_id = xxhash.xxh3_64_intdigest(slice_name, seed=42)
// //     ref_name = aln.reference_name
// //     ref_start = aln.reference_start
// //     ref_end = aln.reference_end
// //     # Check if read mapped
// //     if aln.is_unmapped:
// //         mapped = 0
// //         multimapped = 0
// //         ref_name = ""
// //         ref_start = 0
// //         ref_end = 0
// //         coords = ""
// //     else:
// //         mapped = 1
// //         coords = f"{ref_name}:{ref_start}-{ref_end}"
// //         # Check if multimapped
// //         if aln.is_secondary:
// //             multimapped = 1
// //         else:
// //             multimapped = 0

// //     return CCAlignment(
// //         slice_id=slice_id,
// //         slice_name=slice_name,
// //         parent_id=parent_id,
// //         parent_read=parent_read,
// //         pe=pe,
// //         slice=int(slice_number),
// //         uid=int(uid),
// //         mapped=mapped,
// //         multimapped=multimapped,
// //         chrom=ref_name,
// //         start=int(ref_start),
// //         end=int(ref_end),
// //         coordinates=coords,
// //     )


// // @get_timing(task_name="processing BAM file")
// // def parse_bam(bam: Union[str, pathlib.Path]) -> pd.DataFrame:
// //     """Uses parse_alignment function convert bam file to a dataframe.

// //     Extracts:
// //      -'slice_name'
// //      -'parent_read'
// //      -'pe'
// //      -'slice'
// //      -'mapped'
// //      -'multimapped'
// //      -'chrom'
// //      -'start'
// //      -'end'
// //      -'coordinates'

// //     Args:
// //         bam: Path to bam file.

// //     Returns:
// //      pd.Dataframe: DataFrame with the columns listed above.

// //     """

// //     # Load reads into dataframe
// //     logger.info("Parsing BAM file")
// //     df_bam = pd.DataFrame(
// //         [
// //             parse_alignment(aln)
// //             for aln in pysam.AlignmentFile(bam, "rb").fetch(until_eof=True)
// //         ],
// //     )

// //     # Perform dtype conversions
// //     logger.info("Converting dtypes")
// //     df_bam["chrom"] = df_bam["chrom"].astype("category")
// //     pe_category = pd.CategoricalDtype(["flashed", "pe"])
// //     df_bam["pe"] = df_bam["pe"].astype(
// //         pe_category
// //     )  # Only the one type present so need to include both

// //     df_bam.set_index(["slice_name", "chrom", "start"], inplace=True)

// //     logger.info("Finished parsing BAM file")
// //     return df_bam