use pyo3::prelude::*;
use pyo3::types::{IntoPyDict, PyBytes, PyDict, PyList, PyString, PyTuple};
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
fn get_snvs_dbsnp(
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
fn get_snvs_meticulous<'py>(
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
fn get_snvs_simple<'py> (
    py: Python<'py>,
    query_sequence: &PyString,
    ref_seq: &PyString,
    query_coords: &PyList,
    ref_coords: &PyList,
    ref_coord_start: usize,
    tr_start_pos: usize,
    tr_end_pos: usize,
    entropy_flank_size: usize,
    entropy_threshold: f32,
) -> &'py PyDict {
    let qry_seq_len = query_sequence.len().unwrap();
    let qry_seq_bytes = query_sequence.to_str().unwrap().as_bytes();
    let ref_seq_bytes = ref_seq.to_str().unwrap().as_bytes();

    (0..query_coords.len()).filter_map(|i| {
        let ref_pos = ref_coords.get_item(i).unwrap().extract::<usize>().unwrap();

        if tr_start_pos <= ref_pos && ref_pos < tr_end_pos {
            return None;
        }

        let read_pos = query_coords.get_item(i).unwrap().extract::<usize>().unwrap();

        if qry_seq_bytes[read_pos] == ref_seq_bytes[ref_pos - ref_coord_start] {
            return None;
        }

        let seq = &qry_seq_bytes[read_pos - cmp::min(entropy_flank_size, read_pos)..cmp::min(read_pos + entropy_flank_size, qry_seq_len)];
        (_shannon_entropy(seq) >= entropy_threshold).then(|| (ref_pos, qry_seq_bytes[read_pos] as char))
    }).into_py_dict(py)
}

#[pyfunction]
fn get_read_snvs<'py>(
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
        py,
        query_sequence, 
        ref_seq, 
        query_coords,
        ref_coords,
        ref_coord_start, 
        tr_start_pos, 
        tr_end_pos, 
        entropy_flank_size, 
        entropy_threshold,
    );

    if snvs.keys().len() >= too_many_snvs_threshold {  // TOO MANY, some kind of mismapping going on?
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
        snvs
    }
}

#[pyfunction]
fn get_aligned_pair_matches<'py>(
    py: Python<'py>, 
    cigar: &PyList, 
    query_start: usize, 
    ref_start: usize,
) -> &'py PyTuple {
    let mut qi = query_start;
    let mut di = ref_start;

    // let mut res_vec: Vec<&PyTuple> = vec![];
    let mut qi_vec: Vec<usize> = Vec::new();
    let mut di_vec: Vec<usize> = Vec::new();

    for cigar_op in cigar.iter() {
        let dco0 = cigar_op.get_item(0).unwrap().extract::<usize>().unwrap();

        match dco0 {
            0 | 7 | 8 => {  // MATCH | SEQ_MATCH | SEQ_MISMATCH
                let dco1 = cigar_op.get_item(1).unwrap().extract::<usize>().unwrap();
                qi_vec.extend((0..dco1).map(|i| qi + i));
                di_vec.extend((0..dco1).map(|i| di + i));
                // res_vec.extend((0..(dco1)).map(|i: usize| {
                //     PyTuple::new(py, [qi + i, di + i])
                // }));
                qi += dco1;
                di += dco1;
            },
            1 | 4 => {  // INSERTION | SOFT_CLIPPED
                let dco1 = cigar_op.get_item(1).unwrap().extract::<usize>().unwrap();
                qi += dco1;
            },
            2 | 3 => {  // DELETION | SKIPPED
                let dco1 = cigar_op.get_item(1).unwrap().extract::<usize>().unwrap();
                di += dco1;
            },
            5 | 6 => {  // HARD_CLIPPED | PADDING
                // Do nothing
            }
            _ => panic!("Invalid CIGAR operation")
        }
    }

    PyTuple::new(py, [PyList::new(py, qi_vec), PyList::new(py, di_vec)])
}

#[pymodule]
fn strkit_rust_ext(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(shannon_entropy, m)?)?;
    m.add_function(wrap_pyfunction!(get_snvs_dbsnp, m)?)?;
    m.add_function(wrap_pyfunction!(get_snvs_meticulous, m)?)?;
    m.add_function(wrap_pyfunction!(get_snvs_simple, m)?)?;
    m.add_function(wrap_pyfunction!(get_read_snvs, m)?)?;
    m.add_function(wrap_pyfunction!(get_aligned_pair_matches, m)?)?;
    Ok(())
}
