use pyo3::prelude::*;
use pyo3::types::{IntoPyDict, PyBytes, PyDict, PyList, PyString};
use std::cmp;
use std::collections::HashMap;
use entropy::shannon_entropy as _shannon_entropy;

const SNV_GAP_CHAR: char = '_' as char;

#[pyfunction]
pub fn shannon_entropy(data: &PyBytes) -> f32 {
    _shannon_entropy(data.as_bytes())
}

// The below function is a rewritten version of code written on time paid for by
// McGill University, and is thus (c) McGill University 2023. 
#[pyfunction]
pub fn get_snvs_dbsnp(
    candidate_snv_dict_items_flat: Vec<(usize, &str, &str, Vec<&str>)>,
    query_sequence: &PyString,
    pairs: Vec<(usize, usize)>,
    tr_start_pos: usize,
    tr_end_pos: usize,
) -> HashMap<usize, char> {
    let query_by_ref: HashMap<usize, usize> = pairs.iter().cloned().map(|(a, b)| (b, a)).collect();
    let query_seq_bytes = query_sequence.to_str().unwrap().as_bytes();

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
pub fn get_snvs_meticulous<'py>(
    py: Python<'py>,
    query_sequence: &PyString,
    ref_seq: &PyString,
    query_coords: &PyList,
    ref_coords: &PyList,
    ref_coord_start: usize,
    tr_start_pos: usize,
    tr_end_pos: usize,
    contiguous_threshold: usize,
    max_snv_group_size: usize,
    entropy_flank_size: usize,
    entropy_threshold: f32,
) -> &'py PyDict {
    let qry_seq_len = query_sequence.len().unwrap();

    let qry_seq_bytes = query_sequence.to_str().unwrap().as_bytes();
    let ref_seq_bytes = ref_seq.to_str().unwrap().as_bytes();

    let mut lhs_contiguous: usize = 0;
    let mut rhs_contiguous: usize = 0;
    let mut last_rp: Option<usize> = None;

    let mut snv_group = Vec::<(usize, char)>::new();
    let snvs = PyDict::new(py);

    for i in 0..(ref_coords.len()) {
        let ref_pos = ref_coords.get_item(i).unwrap().extract::<usize>().unwrap();

        if tr_start_pos <= ref_pos && ref_pos < tr_end_pos {
            continue;
        }

        let read_pos = query_coords.get_item(i).unwrap().extract::<usize>().unwrap();

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
                    for &(snv_pos, snv_a) in snv_group.iter() {
                        snvs.set_item(snv_pos, snv_a).unwrap();
                    }
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
            if _shannon_entropy(seq) >= entropy_threshold {
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
        for &(snv_pos, snv_a) in snv_group.iter() {
            snvs.set_item(snv_pos, snv_a).unwrap();
        }
    }

    snvs
}

#[pyfunction]
pub fn get_snvs_simple (
    query_sequence: &PyString,
    ref_seq: &PyString,
    query_coords: &PyList,
    ref_coords: &PyList,
    ref_coord_start: usize,
    tr_start_pos: usize,
    tr_end_pos: usize,
    too_many_snvs_threshold: usize,
    entropy_flank_size: usize,
    entropy_threshold: f32,
) -> HashMap<usize, char> {
    let qry_seq_bytes = query_sequence.to_str().unwrap().as_bytes();
    let qry_seq_len = qry_seq_bytes.len();
    let ref_seq_bytes = ref_seq.to_str().unwrap().as_bytes();

    let mut n_snvs = 0;
    let mut res = HashMap::new();

    for i in 0..query_coords.len() {
        let ref_pos = ref_coords.get_item(i).unwrap().extract::<usize>().unwrap();

        if tr_start_pos <= ref_pos && ref_pos < tr_end_pos {
            continue;
        }

        let read_pos = query_coords.get_item(i).unwrap().extract::<usize>().unwrap();

        if qry_seq_bytes[read_pos] == ref_seq_bytes[ref_pos - ref_coord_start] {
            continue;
        }

        let seq = &qry_seq_bytes[read_pos - cmp::min(entropy_flank_size, read_pos)..cmp::min(read_pos + entropy_flank_size, qry_seq_len)];
        if _shannon_entropy(seq) >= entropy_threshold {
            n_snvs += 1;

            // Below it's '>=' but we can skip one set_item by using '>'
            if too_many_snvs_threshold > 0 && n_snvs > too_many_snvs_threshold {
                return res;
            }

            res.insert(ref_pos, qry_seq_bytes[read_pos] as char);
        }
    }

    res
}

#[pyfunction]
pub fn get_read_snvs<'py>(
    py: Python<'py>,
    query_sequence: &PyString,
    ref_seq: &PyString,
    query_coords: &PyList,
    ref_coords: &PyList,
    ref_coord_start: usize,
    tr_start_pos: usize,
    tr_end_pos: usize,
    contiguous_threshold: usize,
    max_snv_group_size: usize,
    too_many_snvs_threshold: usize,
    entropy_flank_size: usize,
    entropy_threshold: f32,
) -> &'py PyDict {
    // Given a list of tuples of aligned (read pos, ref pos) pairs, this function finds non-reference SNVs which are
    // surrounded by a stretch of aligned bases of a specified size on either side.
    // Returns a hash map of <position, base>

    let snvs = get_snvs_simple(
        query_sequence, 
        ref_seq, 
        query_coords,
        ref_coords,
        ref_coord_start, 
        tr_start_pos, 
        tr_end_pos, 
        too_many_snvs_threshold,
        entropy_flank_size, 
        entropy_threshold,
    );

    if snvs.len() >= too_many_snvs_threshold {  // TOO MANY, some kind of mismapping going on?
        get_snvs_meticulous(
            py,
            query_sequence, 
            ref_seq, 
            query_coords,
            ref_coords,
            ref_coord_start, 
            tr_start_pos, 
            tr_end_pos, 
            contiguous_threshold, 
            max_snv_group_size, 
            entropy_flank_size, 
            entropy_threshold,
        )
    } else {
        snvs.into_py_dict(py)
    }
}
