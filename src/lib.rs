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
fn get_snvs_meticulous(
    query_sequence: &str,
    pairs: Vec<(usize, usize)>,
    ref_seq: &str,
    ref_coord_start: usize,
    tr_start_pos: usize,
    tr_end_pos: usize,
    contiguous_threshold: usize,
    max_snv_group_size: usize,
    entropy_flank_size: usize,
    entropy_threshold: f32,
) -> HashMap<usize, char> {
    let qry_seq_len = query_sequence.len();

    let qry_seq_bytes = query_sequence.as_bytes();
    let ref_seq_bytes = ref_seq.as_bytes();

    let mut lhs_contiguous: usize = 0;
    let mut rhs_contiguous: usize = 0;
    let mut last_rp: Option<usize> = None;

    let mut snv_group = Vec::<(usize, char)>::new();
    let mut snvs = HashMap::<usize, char>::new();

    for (read_pos, ref_pos) in pairs {
        if tr_start_pos <= ref_pos && ref_pos < tr_end_pos {
            continue;
        }

        let read_base = qry_seq_bytes[read_pos];
        let ref_base = ref_seq_bytes[ref_pos - ref_coord_start];

        let contiguous_at_base = match last_rp {
            Some(lr) => contiguous_threshold == 0 || ref_pos - lr == 1,
            None => true,
        };

        if read_base == ref_base && contiguous_at_base {
            let snv_group_len = snv_group.len();

            if snv_group_len > 0 {
                rhs_contiguous += 1;
            } else {
                lhs_contiguous += 1;
            }

            if lhs_contiguous >= contiguous_threshold && rhs_contiguous >= contiguous_threshold {
                if snv_group_len <= max_snv_group_size {
                    snvs.extend(snv_group.iter().cloned());
                }

                // Otherwise, it might be a little mismapped area or a longer deletion vs reference, so ignore it.
                lhs_contiguous = 0;
                rhs_contiguous = 0;
                snv_group.clear();
            }

            last_rp = Some(ref_pos);
            continue;
        }

        if !contiguous_at_base {  // Non-contiguous jump; insertion in query
            lhs_contiguous = 0;
            last_rp = Some(ref_pos);
            continue;
        }

        // Otherwise, contiguous

        if read_base != ref_base {
            // If our entropy is ok, add this to the SNV group
            let seq = &qry_seq_bytes[read_pos - cmp::min(entropy_flank_size, read_pos)..cmp::min(read_pos + entropy_flank_size, qry_seq_len)];
            if shannon_entropy(seq) >= entropy_threshold {
                snv_group.push((ref_pos, read_base as char));
            }

            // Don't reset either contiguous variable; instead, take this as part of a SNP group
            last_rp = Some(ref_pos);
        }
    }

    // Special case: if we have stuff in the SNV group with no contiguous requirements,
    // add it to the SNV HashMap.
    let sgl = snv_group.len();
    if contiguous_threshold == 0 && sgl > 0 && sgl <= max_snv_group_size {
        snvs.extend(snv_group.iter().cloned());
    }

    snvs
}

#[pyfunction]
fn get_snvs_simple(
    query_sequence: &str,
    pairs: Vec<(usize, usize)>,
    ref_seq: &str,
    ref_coord_start: usize,
    tr_start_pos: usize,
    tr_end_pos: usize,
    entropy_flank_size: usize,
    entropy_threshold: f32,
) -> HashMap<usize, char> {
    let qry_seq_len = query_sequence.len();
    let qry_seq_bytes = query_sequence.as_bytes();
    let ref_seq_bytes = ref_seq.as_bytes();

    pairs
        .iter()
        .filter_map(|&(read_pos, ref_pos)| {
            let seq = &qry_seq_bytes[read_pos - cmp::min(entropy_flank_size, read_pos)..cmp::min(read_pos + entropy_flank_size, qry_seq_len)];
            (
                !(tr_start_pos <= ref_pos && ref_pos < tr_end_pos) && 
                (qry_seq_bytes[read_pos] != ref_seq_bytes[ref_pos - ref_coord_start]) && 
                (shannon_entropy(seq) >= entropy_threshold)
            ).then(|| (ref_pos, qry_seq_bytes[read_pos] as char))
        })
        .collect()
}

#[pymodule]
fn strkit_rust_ext(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(get_snvs_meticulous, m)?)?;
    m.add_function(wrap_pyfunction!(get_snvs_simple, m)?)?;
    Ok(())
}
