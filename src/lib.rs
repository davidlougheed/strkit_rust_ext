use pyo3::prelude::*;
use std::collections::HashMap;

#[pyfunction]
fn find_pair_by_ref_pos(pairs: Vec<(usize, usize)>, target: usize) -> (usize, bool) {
    let res = &pairs[..].binary_search_by_key(&target, |&(_qc, rc)| rc);
    return match *res {
        Ok(i)  => (i, true),
        Err(i) => (i, false),
    };
}

#[pyfunction]
fn get_snvs_simple(
    query_sequence: &str,
    pairs: Vec<(usize, usize)>,
    ref_seq: &str,
    ref_coord_start: usize,
    tr_start_pos: usize,
    tr_end_pos: usize,
) -> HashMap<usize, char> {
    let qry_seq_bytes = query_sequence.as_bytes();
    let ref_seq_bytes = ref_seq.as_bytes();

    pairs
        .iter()
        .filter_map(|&(read_pos, ref_pos)| {
            (!(tr_start_pos <= ref_pos && ref_pos < tr_end_pos)
                && (qry_seq_bytes[read_pos] != ref_seq_bytes[ref_pos - ref_coord_start]))
                .then(|| (ref_pos, qry_seq_bytes[read_pos] as char))
        })
        .collect()
}

#[pymodule]
fn strkit_rust_ext(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(find_pair_by_ref_pos, m)?)?;
    m.add_function(wrap_pyfunction!(get_snvs_simple, m)?)?;
    Ok(())
}
