use numpy::{ndarray::Array1, PyArray1, PyArrayMethods};
use pyo3::{prelude::*, pybacked::PyBackedStr};
use std::iter::zip;
use crate::aligned_coords::STRkitAlignedCoords;

pub fn find_coord_idx_by_ref_pos(
    aligned_coords: &STRkitAlignedCoords,
    target: usize,
    start_left: usize,
) -> (usize, bool) {
    let t = target as u64;
    let idx = start_left + aligned_coords.ref_coords[start_left..].partition_point(|&x| x < t);
    let found = idx < aligned_coords.ref_coords.len() && aligned_coords.ref_coords[idx] == t;
    (idx, found)
}

#[pyfunction]
#[pyo3(name = "find_coord_idx_by_ref_pos")]
pub fn find_coord_idx_by_ref_pos_py(
    aligned_coords: &Bound<'_, STRkitAlignedCoords>,
    target: usize,
    start_left: usize,
) -> (usize, bool) {
    find_coord_idx_by_ref_pos(&aligned_coords.borrow(), target, start_left)
}

// TODO: use this
pub fn calculate_seq_with_wildcards(qs: &str, quals: Option<Array1<u8>>, base_wildcard_threshold: u8) -> String {
    if let Some(qls) = quals {
        zip(qs.chars(), qls).map(|(c, q)| {
            if q > base_wildcard_threshold { c } else { 'X' }
        }).collect::<String>()
    } else {
        qs.to_owned()
    }
}

#[pyfunction]
#[pyo3(name = "calculate_seq_with_wildcards")]
pub fn calculate_seq_with_wildcards_py(
    qs: PyBackedStr,
    quals: Option<&Bound<'_, PyArray1<u8>>>,
    base_wildcard_threshold: u8,
) -> String {
    if let Some(qls) = quals {
        let qls_ro = qls.readonly();
        zip(qs.chars(), qls_ro.as_array()).map(|(c, &q)| {
            if q > base_wildcard_threshold { c } else { 'X' }
        }).collect::<String>()
    } else {
        qs.to_string()
    }
}
