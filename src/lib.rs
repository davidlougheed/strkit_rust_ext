mod strkit;

pub use crate::strkit::cigar;
pub use crate::strkit::consensus;
pub use crate::strkit::locus;
pub use crate::strkit::reads;
pub use crate::strkit::repeats;
pub use crate::strkit::snvs;

use pyo3::prelude::*;

#[pymodule]
#[pyo3(name = "strkit_rust_ext")]
fn strkit_rust_ext(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(cigar::get_aligned_pair_matches, m)?)?;

    m.add_function(wrap_pyfunction!(consensus::consensus_seq, m)?)?;

    m.add_function(wrap_pyfunction!(locus::get_read_coords_from_matched_pairs, m)?)?;
    m.add_function(wrap_pyfunction!(locus::get_pairs_and_tr_read_coords, m)?)?;
    m.add_function(wrap_pyfunction!(locus::process_read_snvs_for_locus_and_calculate_useful_snvs, m)?)?;

    m.add_class::<reads::STRkitBAMReader>()?;
    m.add_class::<reads::STRkitAlignedSegment>()?;

    m.add_function(wrap_pyfunction!(repeats::get_repeat_count, m)?)?;

    m.add_class::<snvs::CandidateSNVs>()?;
    m.add_class::<snvs::STRkitVCFReader>()?;
    m.add_function(wrap_pyfunction!(snvs::shannon_entropy, m)?)?;
    m.add_function(wrap_pyfunction!(snvs::get_read_snvs, m)?)?;

    Ok(())
}
