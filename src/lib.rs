use pyo3::prelude::*;
use pythonize::pythonize;

mod fastq_deduplication;
mod genome_digest;
mod utils;

// Deduplicate FASTQ files
#[pyfunction]
#[pyo3(name = "fastq_deduplicate", text_signature = "(fq_in, fq_out, shuffle, error_rate/)")]
fn deduplicate_fastq_py(
    fq_in: Vec<(String, String)>,
    fq_out: Option<Vec<(String, String)>>,
    shuffle: bool,
    error_rate: Option<f32>,
) -> Py<PyAny> {

    // Set up ctrl-c handler
    ctrlc::set_handler(|| std::process::exit(2)).unwrap_or_default();

    // Get the Python GIL
    let gil = Python::acquire_gil();
    let py = gil.python();

    let mut deduplicator = fastq_deduplication::FastqDeduplicator::new(
        fq_in,
        fq_out,
        error_rate,
        shuffle,
    );

    // Identify duplicates
    deduplicator.identify_duplicates().expect("Error identifying duplicates");

    // Remove duplicates
    let deduplication_results = deduplicator.remove_duplicate_sequences().expect("Error removing duplicate sequences");

    // Convert statistics to Python
    let py_deduplication_results = pythonize(py, &deduplication_results).unwrap();
    py_deduplication_results


}

#[pymodule]
#[pyo3(name = "capcruncher_tools")]
fn capcruncher_tools(_py: Python, m: &PyModule) -> PyResult<()> {
    pyo3_log::init();
    m.add_function(wrap_pyfunction!(deduplicate_fastq_py, m)?)?;
    Ok(())
}
