use hashbrown::HashMap;
use pyo3::prelude::*;
use pythonize::pythonize;

mod fastq_deduplication;

#[pyfunction]
#[pyo3(name = "deduplicate_fastq")]
fn deduplicate_fastq_py(
    fq_in: HashMap<String, Vec<String>>,
    fq_out: HashMap<String, Vec<String>>,
    shuffle: bool,
) -> Py<PyAny> {

    ctrlc::set_handler(|| std::process::exit(2)).unwrap_or_default();
    let gil = Python::acquire_gil();
    let py = gil.python();

    let deduplication_results = fastq_deduplication::deduplicate_fastq(fq_in, fq_out, shuffle).unwrap();
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
