use bio::io;
use bio::pattern_matching::bom::BOM;
use fastq::{each_zipped, Parser, Record};
use hashbrown::{HashMap, HashSet};
use indicatif::ParallelProgressIterator;
use log::{debug, info, warn};
use rand::prelude::*;
use rayon::prelude::*;

use crate::utils::get_fastq_reader_file_handles;

// class ReadDigestionProcess(multiprocessing.Process):
//     """
//     Process subclass for multiprocessing fastq digestion.

//     """

//     def __init__(
//         self,
//         inq: multiprocessing.Queue,
//         outq: multiprocessing.Queue,
//         stats_pipe: multiprocessing.Pipe = None,
//         **digestion_kwargs,
//     ) -> None:

//         """
//         Args:
//          inq (multiprocessing.SimpleQueue): Queue to hold list of reads to digest.
//          outq (multiprocessing.SimpleQueue): Queue to hold list of digested reads.
//          stats_pipe (multiprocessing.Pipe, optional): Pipe to send statistics to. Defaults to None.
//          **digestion_kwargs Dict[Any, Any]: Kwargs passed to DigestedRead.

//         Raises:
//             KeyError: digestion_kwargs must contain: cutsite
//         """

//         super(ReadDigestionProcess, self).__init__()

//         self.inq = inq
//         self.outq = outq
//         self.stats_pipe = stats_pipe
//         self.digestion_kwargs = digestion_kwargs
//         self.read_type = digestion_kwargs.get("read_type", "flashed")
//         self._stat_container = DigestionStats

//         if "cutsite" not in digestion_kwargs:
//             raise KeyError("Cutsite is required to be present in digestion arguments")

//     def _digest_reads(self, reads, **digestion_kwargs):
//         digested = []
//         for i, read in enumerate(reads):
//             if i == 0:
//                 digested.append(DigestedRead(read, **digestion_kwargs))
//             else:
//                 digestion_kwargs["slice_number_start"] = digested[i - 1].slices_filtered
//                 digested.append(DigestedRead(read, **digestion_kwargs))

//         return digested

//     def _digestion_statistics_to_dataframe(self, stats: List[DigestionStats]):

//         if stats:

//             df = pd.DataFrame(stats)
//             return (
//                 df.groupby(["read_type", "read_number", "unfiltered", "filtered"])
//                 .size()
//                 .to_frame("count")
//                 .reset_index()
//             )

//     def run(self):
//         """
//         Performs read digestion.

//         Reads to digest are pulled from inq, digested with the DigestedRead class
//         and the results placed on outq for writing.

//         If a statq is provided, read digestion stats are placed into this queue for
//         aggregation.

//         """

//         buffer_stats = list()
//         buffer_reads = list()
//         dframes_stats = list()

//         while True:
//             try:
//                 reads = self.inq.get(block=True, timeout=0.01)

//                 # Make sure that we don't need to terminate
//                 if reads:

//                     # Accounts for PE as well as flashed
//                     for read in reads:

//                         # Digest the read
//                         digested = self._digest_reads(read, **self.digestion_kwargs)

//                         # Parse the digestion results
//                         for read_number, digested_read in enumerate(digested):

//                             # Only write if there are valid slices present
//                             if digested_read.has_valid_slices:
//                                 buffer_reads.append(str(digested_read))

//                             # Will record all reads even if these do not digest
//                             digested_read_stats = self._stat_container(
//                                 read_type=self.read_type,
//                                 read_number=(
//                                     read_number + 1
//                                     if not self.read_type == "flashed"
//                                     else read_number
//                                 ),
//                                 unfiltered=digested_read.slices_unfiltered,
//                                 filtered=digested_read.slices_filtered,
//                             )

//                             # Append stats to the stats buffer
//                             buffer_stats.append(digested_read_stats)

//                     # Aggregate individual read stats into a dataframe
//                     df_stat_batch = self._digestion_statistics_to_dataframe(
//                         buffer_stats
//                     )

//                     # Add this summary to the overall stats
//                     dframes_stats.append(df_stat_batch)

//                     # Add the digested reads to the output queue
//                     self.outq.put("".join(buffer_reads.copy()))
//                     buffer_reads.clear()
//                     buffer_stats.clear()

//                 else:
//                     break

//             except queue.Empty:
//                 continue

//         if self.stats_pipe:
//             # Merge all dataframes together

//             try:
//                 df_stats = pd.concat(dframes_stats)
//                 df_stats_aggregated = (
//                     df_stats.groupby(
//                         ["read_type", "read_number", "unfiltered", "filtered"]
//                     )
//                     .sum()
//                     .reset_index()
//                 )

//                 # Send the statistics to the main process
//                 self.stats_pipe.send(df_stats_aggregated)
//             except ValueError:
//                 # Might not actually have got any reads, just return none
//                 self.stats_pipe.send(None)

pub enum ReadType {
    Flashed,
    Pe,
}

impl std::fmt::Display for ReadType {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            ReadType::Flashed => write!(f, "flashed"),
            ReadType::Pe => write!(f, "pe"),
        }
    }
}

struct DigestableRead<R> {
    read: R,
    restriction_site: String,
    min_slice_length: usize,
    slice_number_start: usize,
    allow_undigested: bool,
    read_type: ReadType,
    n_slices_unfiltered: usize,
    n_slices_filtered: usize,
}

impl<R> DigestableRead<R> {
    pub fn new(
        read: R,
        restriction_site: String,
        allow_undigested: bool,
        read_type: ReadType,
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
}
impl<R: Record> DigestableRead<R> {
    fn digest_indicies(&mut self) -> Vec<usize> {
        let sequence = self.read.seq().to_ascii_lowercase().to_vec();
        let matcher = BOM::new(self.restriction_site.to_ascii_lowercase().as_bytes());
        let mut slice_indexes: Vec<_> = matcher.find_all(&sequence).collect();

        // Add the start and end to slices
        slice_indexes.insert(0, 0);
        slice_indexes.push(sequence.len());

        // Remove all slices that are smaller than the minimum slice size

        slice_indexes = slice_indexes
            .into_iter()
            .filter(|&x| x >= self.min_slice_length)
            .collect();

        slice_indexes
    }

    pub fn digest(&mut self) -> Vec<bio::io::fastq::Record> {
        let mut slices = Vec::new();
        let mut slice_no = self.slice_number_start;

        let slice_indexes = self.digest_indicies();

        if slice_indexes.len() > 1 || self.allow_undigested {
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

    fn validate_slice(&self, slice: &bio::io::fastq::Record) -> bool {
        slice.seq().len() >= self.min_slice_length
    }

    fn prepare_slice(&self, start: usize, end: usize, slice_no: usize) -> bio::io::fastq::Record {
        let read_name = format!(
            "@{}|{}|{}|{}",
            std::str::from_utf8(self.read.head()).expect("Invalid UTF-8 in read name"),
            self.read_type,
            slice_no,
            thread_rng().gen_range(0..100)
        );

        let sequence = self.read.seq()[start..end].to_ascii_uppercase();
        let quality = self.read.qual()[start..end].to_owned();

        bio::io::fastq::Record::with_attrs(&read_name, None, &sequence, &quality)
    }
}

struct DigestableFastq {
    fastq: String,
    restriction_site: String,
    allow_undigested: bool,
    read_type: ReadType,
    min_slice_length: Option<usize>,
}

fn digest_fastq(
    fastqs: Vec<String>,
    restriction_site: String,
    read_type: ReadType,
    min_slice_length: Option<usize>,
) -> Result<(), std::io::Error> {



    Ok(())





}
