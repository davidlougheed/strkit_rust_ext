use pyo3::prelude::*;
use std::cmp;
use std::collections::HashMap;

fn shannon_entropy(seq: &[u8]) -> f32 {
    let base_counts: HashMap<u8, i32> = seq
        .iter()
        .fold(HashMap::new(), |mut map, &b| {
            *map.entry(b).or_insert(0) += 1;
            map
        });

    -1.0 * base_counts.values().map(|&c| {
        let p = c as f32 / seq.len() as f32;
        p * p.log2()
    }).sum::<f32>()
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
    m.add_function(wrap_pyfunction!(get_snvs_simple, m)?)?;
    Ok(())
}
