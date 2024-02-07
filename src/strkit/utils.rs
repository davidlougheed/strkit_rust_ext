pub fn find_coord_idx_by_ref_pos(r_coords: &Vec<usize>, target: usize, start_left: usize) -> (usize, bool) {
    let idx = start_left + r_coords[start_left..].partition_point(|&x| x < target);
    let found = idx < r_coords.len() && r_coords[idx] == target; 
    (idx, found)
}
