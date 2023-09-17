use anyhow::Ok;
use polars::prelude::*;
use pyo3::prelude::*;
use pyo3_polars::PyDataFrame;
use pythonize::pythonize;
use std::str::FromStr;

mod alignment;
mod fastq_deduplication;
mod fastq_digest;
mod genome_digest;
mod interactions_count;
mod utils;

// Rust based. Deduplicate FASTQ files based on exact sequence matches. Returns a dictionary with statistics."
#[pyfunction]
#[pyo3(
    name = "fastq_deduplicate",
    text_signature = "(fq_in, fq_out, shuffle=False)",
    signature = (fq_in, fq_out, shuffle=false)
)]
fn deduplicate_fastq_py(
    fq_in: Vec<(String, String)>,
    fq_out: Option<Vec<(String, String)>>,
    shuffle: bool,
) -> Py<PyAny> {
    // Set up ctrl-c handler
    ctrlc::set_handler(|| std::process::exit(2)).unwrap_or_default();

    // // Get the Python GIL
    // let gil = Python::acquire_gil();
    // let py = gil.python();

    let mut deduplicator = fastq_deduplication::FastqDeduplicator::new(fq_in, fq_out, shuffle);

    // Run the deduplication
    let deduplication_results = deduplicator
        .write_unique_reads()
        .expect("Error during deduplication");

    // Convert statistics to Python
    let py_deduplication_results =
        Python::with_gil(|py| pythonize(py, &deduplication_results).unwrap());
    py_deduplication_results
}

// Rust based. Digest a FASTA file with a restriction enzyme. Returns a BED file with the digested fragments.
#[pyfunction]
#[pyo3(
    name = "digest_fasta",
    text_signature = "(fasta, restriction_site, output, remove_recognition_site, min_slice_length)"
)]
fn digest_fasta_py(
    fasta: String,
    restriction_site: String,
    output: String,
    remove_recognition_site: bool,
    min_slice_length: Option<usize>,
    n_threads: Option<usize>,
) -> PyResult<()> {
    // Set up ctrl-c handler
    ctrlc::set_handler(|| std::process::exit(2)).unwrap_or_default();

    // Run the digest
    genome_digest::digest_fasta(
        fasta,
        restriction_site,
        output,
        remove_recognition_site,
        min_slice_length,
        n_threads,
    )?;

    Result::Ok(())
}

#[pyfunction]
#[pyo3(
    name = "digest_fastq",
    text_signature = "(fastq, restriction_site, output, read_type, sample, min_slice_length)"
)]
fn digest_fastq_py(
    fastqs: Vec<String>,
    output: String,
    restriction_site: String,
    read_type: String,
    sample: String,
    min_slice_length: Option<usize>,
) -> PyResult<(PyDataFrame, PyDataFrame, PyDataFrame, PyDataFrame)> {
    // Set up ctrl-c handler
    ctrlc::set_handler(|| std::process::exit(2)).unwrap_or_default();

    // Run the digest
    let res = fastq_digest::digest_fastq(
        fastqs,
        output,
        restriction_site.to_lowercase(),
        fastq_digest::ReadType::from_str(&read_type).expect("Invalid read type"),
        min_slice_length,
        Some(sample),
    );

    match res {
        Result::Ok((read_stats, hist_unfilt, hist_filt, hist_len)) => {
            // Create a DataFrame with the read statistics
            let df_read = read_stats.to_dataframe();
            let df_read_py = PyDataFrame(df_read);

            // Create a DataFrame with the unfiltered histogram
            let df_hist_unfilt = hist_unfilt
                .iter()
                .map(|h| h.to_dataframe(Some("slice_number")))
                .reduce(|a, b| a.vstack(&b).unwrap())
                .expect("Error during histogram creation");
            let df_hist_unfilt_py = PyDataFrame(df_hist_unfilt);

            // Create a DataFrame with the filtered histogram
            let df_hist_filt = hist_filt
                .iter()
                .map(|h| h.to_dataframe(Some("slice_number")))
                .reduce(|a, b| a.vstack(&b).unwrap())
                .expect("Error during histogram creation");
            let df_hist_filt_py = PyDataFrame(df_hist_filt);

            // Create a DataFrame with the length histogram
            let df_hist_len = hist_len
                .iter()
                .map(|h| h.to_dataframe(Some("slice_length")))
                .reduce(|a, b| a.vstack(&b).unwrap())
                .expect("Error during histogram creation");
            let df_hist_len_py = PyDataFrame(df_hist_len);

            // Return a tuple of the DataFrames
            Result::Ok((
                df_read_py,
                df_hist_unfilt_py,
                df_hist_filt_py,
                df_hist_len_py,
            ))
        }
        Err(e) => {
            println!("Error: {}", e);
            std::process::exit(1);
        }
    }
}

#[pyfunction]
#[pyo3(name = "count_interactions", text_signature = "(df: DataFrame)")]
fn count_interactions(df: PyDataFrame) -> PyDataFrame {
    ctrlc::set_handler(|| std::process::exit(2)).unwrap_or_default();
    let df = interactions_count::count(df.into());
    df
}

#[pymodule]
#[pyo3(name = "capcruncher_tools")]
fn capcruncher_tools(_py: Python, m: &PyModule) -> PyResult<()> {
    // Initialize the logger
    pyo3_log::init();

    // Create a submodule
    let deduplicate = PyModule::new(_py, "deduplicate")?;
    deduplicate.add_function(wrap_pyfunction!(deduplicate_fastq_py, m)?)?;
    m.add_submodule(deduplicate)?;

    // Create a submodule
    let digest = PyModule::new(_py, "digest")?;
    digest.add_function(wrap_pyfunction!(digest_fasta_py, m)?)?;
    digest.add_function(wrap_pyfunction!(digest_fastq_py, m)?)?;
    m.add_submodule(digest)?;

    // Create a submodule
    let interactions = PyModule::new(_py, "interactions")?;
    interactions.add_function(wrap_pyfunction!(count_interactions, m)?)?;
    m.add_submodule(interactions)?;

    Result::Ok(())
}
