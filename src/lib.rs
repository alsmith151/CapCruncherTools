use hashbrown::HashMap;
use pyo3::prelude::*;
use pythonize::pythonize;

mod fastq_deduplication;

#[pyfunction]
#[pyo3(name = "deduplicate_fastq")]
fn deduplicate_fastq_py(fq_in: HashMap<String, Vec<String>>, fq_out: HashMap<String, Vec<String>>) -> Py<PyAny>{
    
    let gil = Python::acquire_gil();
    let py = gil.python();

    let deduplication_results = fastq_deduplication::deduplicate_fastq(fq_in, fq_out).unwrap();
    let py_deduplication_results = pythonize(py, &deduplication_results).unwrap();
    py_deduplication_results
}

#[pymodule]
#[pyo3(name = "capcruncher_tools")]
fn rscapcruncher(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(deduplicate_fastq_py, m)?)?;
    Ok(())
}
