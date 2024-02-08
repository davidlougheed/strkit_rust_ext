pub fn find_coord_idx_by_ref_pos(r_coords: &[u64], target: usize, start_left: usize) -> (usize, bool) {
    let t = target as u64;
    let idx = start_left + r_coords[start_left..].partition_point(|&x| x < t);
    let found = idx < r_coords.len() && r_coords[idx] == t; 
    (idx, found)
}
