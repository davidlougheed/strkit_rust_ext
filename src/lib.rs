mod strkit;

pub use crate::strkit::cigar;
pub use crate::strkit::consensus;
pub use crate::strkit::locus;
pub use crate::strkit::snvs;

use pyo3::prelude::*;

#[pymodule]
fn strkit_rust_ext(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(cigar::get_aligned_pair_matches, m)?)?;
    m.add_function(wrap_pyfunction!(consensus::consensus_seq, m)?)?;
    m.add_function(wrap_pyfunction!(locus::get_pairs_and_tr_read_coords, m)?)?;
    m.add_function(wrap_pyfunction!(locus::process_read_snvs_for_locus_and_calculate_useful_snvs, m)?)?;
    m.add_function(wrap_pyfunction!(snvs::shannon_entropy, m)?)?;
    m.add_function(wrap_pyfunction!(snvs::get_read_snvs, m)?)?;
    Ok(())
}
