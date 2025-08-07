mod strkit;

pub use crate::strkit::aligned_coords;
pub use crate::strkit::cigar;
pub use crate::strkit::consensus;
pub use crate::strkit::locus;
pub use crate::strkit::reads;
pub use crate::strkit::repeats;
pub use crate::strkit::snvs;
pub use crate::strkit::utils;

use pyo3::prelude::*;

#[pymodule]
#[pyo3(name = "strkit_rust_ext")]
fn strkit_rust_ext(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<aligned_coords::STRkitAlignedCoords>()?;

    m.add_function(wrap_pyfunction!(cigar::get_aligned_pair_matches, m)?)?;

    m.add_function(wrap_pyfunction!(consensus::consensus_seq, m)?)?;

    m.add_class::<locus::STRkitLocus>()?;
    m.add_class::<locus::STRkitLocusWithRefData>()?;
    m.add_function(wrap_pyfunction!(locus::get_read_coords_from_matched_pairs, m)?)?;
    m.add_function(wrap_pyfunction!(locus::get_pairs_and_tr_read_coords, m)?)?;
    m.add_function(wrap_pyfunction!(locus::process_read_snvs_for_locus_and_calculate_useful_snvs, m)?)?;

    m.add_class::<reads::STRkitBAMReader>()?;
    m.add_class::<reads::STRkitAlignedSegment>()?;
    m.add_class::<reads::STRkitLocusBlockSegments>()?;

    m.add_function(wrap_pyfunction!(repeats::get_repeat_count, m)?)?;

    m.add_class::<snvs::CandidateSNVs>()?;
    m.add_class::<snvs::STRkitVCFReader>()?;
    m.add_function(wrap_pyfunction!(snvs::shannon_entropy, m)?)?;

    m.add_function(wrap_pyfunction!(utils::find_coord_idx_by_ref_pos_py, m)?)?;
    m.add_function(wrap_pyfunction!(utils::calculate_seq_with_wildcards_py, m)?)?;

    Ok(())
}
