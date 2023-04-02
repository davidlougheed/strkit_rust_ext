use pyo3::prelude::*;
use pyo3::types::PyDict;
use std::collections::HashMap;

#[pyfunction]
fn mk_snvs_dict<'a>(
    py: Python<'a>,
    query_sequence: &str,
    pairs: Vec<(usize, usize)>,
    ref_seq: &str,
    ref_coord_start: usize,
    tr_start_pos: usize,
    tr_end_pos: usize,
) -> PyResult<&'a PyDict> {
    let qry_seq_bytes = query_sequence.as_bytes();
    let ref_seq_bytes = ref_seq.as_bytes();

    let out = PyDict::new(py);

    for (read_pos, ref_pos) in pairs {
        if !(tr_start_pos <= ref_pos && ref_pos < tr_end_pos)
            && (qry_seq_bytes[read_pos] != ref_seq_bytes[ref_pos - ref_coord_start])
        {
            out.set_item(ref_pos, qry_seq_bytes[read_pos] as char)?;
        }
    }

    Ok(out)
}

#[pyfunction]
fn mk_snvs_hash(
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
    m.add_function(wrap_pyfunction!(mk_snvs_dict, m)?)?;
    m.add_function(wrap_pyfunction!(mk_snvs_hash, m)?)?;
    Ok(())
}
