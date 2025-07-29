use pyo3::prelude::*;
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
