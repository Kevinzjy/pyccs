use pyo3::prelude::*;
use pyo3::wrap_pyfunction;

use crate::scan::find_concensus;

#[pyfunction]
fn ccs(
    seq: &[u8],
) -> Result<(String, String), PyErr> {
    match find_concensus(seq) {
        Err(_) => return Ok((String::new(), String::new())),//
        Ok(x) => {
            let ccs_seq = x.1;
            let mut segments: Vec<String> = Vec::default();
            for (s, e) in x.0 {
                segments.push(String::from(vec![s.to_string(), e.to_string()].join("-")));
            }
            let segment_str = segments.join(";");
            return Ok((segment_str, ccs_seq))
        }
    };
}

#[pymodule]
fn pyccs(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(ccs, m)?)?;
    Ok(())
}
