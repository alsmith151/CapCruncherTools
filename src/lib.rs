use pyo3::prelude::*;
use pythonize::pythonize;

mod genome_digest;
mod utils;
mod fastq_deduplication;

// Rust based. Deduplicate FASTQ files based on exact sequence matches. Returns a dictionary with statistics."
#[pyfunction]
#[pyo3(
    name = "fastq_deduplicate",
    text_signature = "(fq_in, fq_out, shuffle)",
)]
fn deduplicate_fastq_py(
    fq_in: Vec<(String, String)>,
    fq_out: Option<Vec<(String, String)>>,
    shuffle: bool,
) -> Py<PyAny> {
    // Set up ctrl-c handler
    ctrlc::set_handler(|| std::process::exit(2)).unwrap_or_default();

    // Get the Python GIL
    let gil = Python::acquire_gil();
    let py = gil.python();

    let mut deduplicator = fastq_deduplication::FastqDeduplicator::new(
        fq_in,
        fq_out,
        shuffle,
    );

    // Run the deduplication
    let deduplication_results = deduplicator.write_unique_reads().expect("Error during deduplication");

    // Convert statistics to Python
    let py_deduplication_results = pythonize(py, &deduplication_results).unwrap();
    py_deduplication_results
}

// Rust based. Digest a FASTA file with a restriction enzyme. Returns a BED file with the digested fragments.
#[pyfunction]
#[pyo3(
    name = "digest_fasta",
    text_signature = "(fasta, restriction_site, output, remove_recognition_site, min_slice_length)",
)]
fn digest_fasta_py(
    fasta: String,
    restriction_site: String,
    output: String,
    remove_recognition_site: bool,
    min_slice_length: Option<usize>,
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
    )?;

    Ok(())
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
    m.add_submodule(digest)?;
    

    Ok(())
}
