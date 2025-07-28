use numpy::{PyArray2, PyArrayMethods};
use pyo3::prelude::*;

use crate::aligned_coords::STRkitAlignedCoords;

pub fn get_aligned_pair_matches_rs(
    cigar: &Bound<'_, PyArray2<u32>>,
    query_start: u64,
    ref_start: u64,
) -> STRkitAlignedCoords {
    let mut qi = query_start;
    let mut di = ref_start;

    // Give these vectors an initial capacity that roughly matches a bit below normal HiFi read length.
    // I'm sure there's a smarter way to do this, but this is easy and fast...
    let mut qi_vec: Vec<u64> = Vec::with_capacity(11000);
    let mut di_vec: Vec<u64> = Vec::with_capacity(11000);

    let cigar_ro = cigar.readonly();
    let cigar_arr = cigar_ro.as_array();

    for cigar_op_idx in 0..cigar_arr.shape()[0] {
        let dco0 = cigar_arr[[cigar_op_idx, 0]];

        match dco0 {
            0 | 7 | 8 => {  // MATCH | SEQ_MATCH | SEQ_MISMATCH
                let dco1 = cigar_arr[[cigar_op_idx, 1]] as u64;

                qi_vec.extend(qi..qi+dco1);
                di_vec.extend(di..di+dco1);

                qi += dco1;
                di += dco1;
            },
            1 | 4 => {  // INSERTION | SOFT_CLIPPED
                let dco1 = cigar_arr[[cigar_op_idx, 1]] as u64;
                qi += dco1;
            },
            2 | 3 => {  // DELETION | SKIPPED
                let dco1 = cigar_arr[[cigar_op_idx, 1]] as u64;
                di += dco1;
            },
            _ => {  // HARD_CLIPPED | PADDING | (unknown cigar op)
                // Do nothing
            }
        }
    }

    qi_vec.shrink_to_fit();
    di_vec.shrink_to_fit();

    STRkitAlignedCoords { query_coords: qi_vec, ref_coords: di_vec }
}


#[pyfunction]
pub fn get_aligned_pair_matches<'py>(
    py: Python<'py>,
    cigar: &Bound<'py, PyArray2<u32>>,
    query_start: u64,
    ref_start: u64,
) -> PyResult<Py<STRkitAlignedCoords>> {
    let aligned_coords = get_aligned_pair_matches_rs(cigar, query_start, ref_start);
    Py::new(py, aligned_coords)
}
