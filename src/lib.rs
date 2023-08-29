use pyo3::prelude::*;
use std::cmp;
use std::collections::HashMap;

const SNV_GAP_CHAR: char = '_' as char;

#[pyfunction]
fn shannon_entropy(seq: &[u8]) -> f32 {
    let seq_len = seq.len() as f32;

    let base_counts: HashMap<u8, i32> = seq
        .iter()
        .fold(HashMap::new(), |mut map, &b| {
            *map.entry(b).or_insert(0) += 1;
            map
        });

    -1.0 * base_counts.values().map(|&c| {
        let p = c as f32 / seq_len;
        p * p.log2()
    }).sum::<f32>()
}

// The below function is a rewritten version of code written on time paid for by
// McGill University, and is thus (c) McGill University 2023. 
#[pyfunction]
fn get_snvs_dbsnp(
    candidate_snv_dict_items_flat: Vec<(usize, &str, &str, Vec<&str>)>,
    query_sequence: &str,
    pairs: Vec<(usize, usize)>,
    tr_start_pos: usize,
    tr_end_pos: usize,
) -> HashMap<usize, char> {
    let query_by_ref: HashMap<usize, usize> = pairs.iter().cloned().map(|(a, b)| (b, a)).collect();
    let query_seq_bytes = query_sequence.as_bytes();

    let mapped_l = pairs.first().unwrap().1;
    let mapped_r = pairs.last().unwrap().1;

    let mut snvs = HashMap::<usize, char>::new();

    for (pos, _id, snv_ref, snv_alts) in candidate_snv_dict_items_flat {
        if pos < mapped_l {
            continue;
        } else if pos > mapped_r {
            break;
        } else if tr_start_pos <= pos && pos <= tr_end_pos {
            continue;
        }

        let read_base = match query_by_ref.get(&pos) {
            Some(&q_pos) => query_seq_bytes[q_pos] as char,
            None => SNV_GAP_CHAR,
        };

        if read_base == snv_ref.chars().next().unwrap() || snv_alts.iter().any(|&snv_alt| {
            read_base == snv_alt.chars().next().unwrap()
        }) {
            snvs.insert(pos, read_base);
        }
    }

    snvs
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

fn _get_snvs_simple(
    query_sequence: &str,
    pairs: &Vec<(usize, usize)>,
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
    // Wrapper function for _get_snvs_simple for both Python binding and to 
    // borrow the Vec of pairs (from ourselves) for the inner function
    _get_snvs_simple(
        query_sequence, 
        &pairs, 
        ref_seq, 
        ref_coord_start,
        tr_start_pos, 
        tr_end_pos, 
        entropy_flank_size,
        entropy_threshold,
    )
}

#[pyfunction]
fn get_read_snvs(
    query_sequence: &str,
    pairs: Vec<(usize, usize)>,
    ref_seq: &str,
    ref_coord_start: usize,
    tr_start_pos: usize,
    tr_end_pos: usize,
    contiguous_threshold: usize,
    max_snv_group_size: usize,
    too_many_snvs_threshold: usize,
    entropy_flank_size: usize,
    entropy_threshold: f32,
) -> HashMap<usize, char> {
    // Given a list of tuples of aligned (read pos, ref pos) pairs, this function finds non-reference SNVs which are
    // surrounded by a stretch of aligned bases of a specified size on either side.
    // Returns a hash map of <position, base>

    let snvs = _get_snvs_simple(
        query_sequence, 
        &pairs, 
        ref_seq, 
        ref_coord_start, 
        tr_start_pos, 
        tr_end_pos, 
        entropy_flank_size, 
        entropy_threshold,
    );

    if snvs.keys().len() >= too_many_snvs_threshold {  // TOO MANY, some kind of mismapping going on?
        get_snvs_meticulous(
            query_sequence, 
            pairs, 
            ref_seq, 
            ref_coord_start, 
            tr_start_pos, 
            tr_end_pos, 
            contiguous_threshold, 
            max_snv_group_size, 
            entropy_flank_size, 
            entropy_threshold,
        )
    } else {
        snvs
    }
}

#[pymodule]
fn strkit_rust_ext(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(shannon_entropy, m)?)?;
    m.add_function(wrap_pyfunction!(get_snvs_dbsnp, m)?)?;
    m.add_function(wrap_pyfunction!(get_snvs_meticulous, m)?)?;
    m.add_function(wrap_pyfunction!(get_snvs_simple, m)?)?;
    m.add_function(wrap_pyfunction!(get_read_snvs, m)?)?;
    Ok(())
}
