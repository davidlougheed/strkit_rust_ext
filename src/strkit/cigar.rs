use numpy::{PyArray1, PyArrayMethods};
use numpy::ndarray::ArrayView1;
use pyo3::prelude::*;

use crate::aligned_coords::STRkitAlignedCoords;

pub fn decode_cigar_item(cigar_item: u32) -> (u32, u64) {
    (cigar_item & 0b1111, (cigar_item >> 4).into())
}

pub fn get_aligned_pair_matches_rs(
    cigar: ArrayView1<u32>,
    query_start: u64,
    ref_start: u64,
) -> STRkitAlignedCoords {
    let mut qi = query_start;
    let mut di = ref_start;

    // Give these vectors an initial capacity that roughly matches a bit below normal HiFi read length.
    // I'm sure there's a smarter way to do this, but this is easy and fast...
    let mut qi_vec: Vec<u64> = Vec::with_capacity(11000);
    let mut di_vec: Vec<u64> = Vec::with_capacity(11000);

    for &cigar_item in cigar.iter() {
        let (cigar_op, len) = decode_cigar_item(cigar_item);

        match cigar_op {
            0 | 7 | 8 => {  // MATCH | SEQ_MATCH | SEQ_MISMATCH
                qi_vec.extend(qi..qi+len);
                di_vec.extend(di..di+len);

                qi += len;
                di += len;
            },
            1 | 4 => {  // INSERTION | SOFT_CLIPPED
                qi += len;
            },
            2 | 3 => {  // DELETION | SKIPPED
                di += len;
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
    cigar: &Bound<'py, PyArray1<u32>>,
    query_start: u64,
    ref_start: u64,
) -> PyResult<Py<STRkitAlignedCoords>> {
    let cigar_ro = cigar.readonly();
    let cigar_view = cigar_ro.as_array();
    let aligned_coords = get_aligned_pair_matches_rs(cigar_view, query_start, ref_start);
    Py::new(py, aligned_coords)
}
