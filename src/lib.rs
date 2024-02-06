mod strkit;

pub use crate::strkit::cigar;
pub use crate::strkit::consensus;
pub use crate::strkit::snvs;

use pyo3::prelude::*;

#[pymodule]
fn strkit_rust_ext(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(snvs::shannon_entropy, m)?)?;
    m.add_function(wrap_pyfunction!(snvs::get_snvs_dbsnp, m)?)?;
    m.add_function(wrap_pyfunction!(snvs::get_snvs_meticulous, m)?)?;
    m.add_function(wrap_pyfunction!(snvs::get_snvs_simple, m)?)?;
    m.add_function(wrap_pyfunction!(snvs::get_read_snvs, m)?)?;
    m.add_function(wrap_pyfunction!(cigar::get_aligned_pair_matches, m)?)?;
    m.add_function(wrap_pyfunction!(consensus::consensus_seq, m)?)?;
    Ok(())
}
