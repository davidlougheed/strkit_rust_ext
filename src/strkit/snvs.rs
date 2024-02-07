use pyo3::prelude::*;
use pyo3::types::{IntoPyDict, PyBytes, PyDict, PyList, PyString};
use std::cmp;
use std::collections::{HashMap, HashSet};
use entropy::shannon_entropy as _shannon_entropy;

use crate::strkit::utils::find_coord_idx_by_ref_pos;

static SNV_OUT_OF_RANGE_CHAR: char = '-';
static SNV_GAP_CHAR: char = '_';


// We check entropy against a threshold in order to make sure the SNVs we find are somewhat
// useful, i.e., surrounded by a nice mixture of different bases rather than some possibly
// mis-mappable base inside a homopolymer (I used to get results like AAAAAAAAAACAAAAAAAA
// where one of those As could be a gap instead... really useless stuff that made it
// through the filters.)

#[pyfunction]
pub fn shannon_entropy(data: &PyBytes) -> f32 {
    _shannon_entropy(data.as_bytes())
}

pub fn get_snvs_meticulous(
    query_sequence: &str,
    ref_seq: &str,
    query_coords: &Vec<usize>,
    ref_coords: &Vec<usize>,
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

    let mut snvs = HashMap::new();

    for i in 0..(ref_coords.len()) {
        let ref_pos = ref_coords[i];

        if tr_start_pos <= ref_pos && ref_pos < tr_end_pos {
            continue;
        }

        let read_pos = query_coords[i];

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
                        snvs.insert(snv_pos, snv_a);
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
            snvs.insert(snv_pos, snv_a);
        }
    }

    snvs
}

pub fn get_snvs_simple (
    query_sequence: &str,
    ref_seq: &str,
    query_coords: &Vec<usize>,
    ref_coords: &Vec<usize>,
    ref_coord_start: usize,
    tr_start_pos: usize,
    tr_end_pos: usize,
    too_many_snvs_threshold: usize,
    entropy_flank_size: usize,
    entropy_threshold: f32,
) -> HashMap<usize, char> {
    let qry_seq_bytes = query_sequence.as_bytes();
    let qry_seq_len = qry_seq_bytes.len();
    let ref_seq_bytes = ref_seq.as_bytes();

    let mut n_snvs = 0;
    let mut res = HashMap::new();

    for i in 0..query_coords.len() {
        let &ref_pos = ref_coords.get(i).unwrap();

        if tr_start_pos <= ref_pos && ref_pos < tr_end_pos {
            continue;
        }

        let &read_pos = query_coords.get(i).unwrap();

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

pub fn get_read_snvs_rs(
    query_sequence: &str,
    ref_seq: &str,
    query_coords: Vec<usize>,
    ref_coords: Vec<usize>,
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

    let snvs = get_snvs_simple(
        query_sequence, 
        ref_seq, 
        &query_coords,
        &ref_coords,
        ref_coord_start, 
        tr_start_pos, 
        tr_end_pos, 
        too_many_snvs_threshold,
        entropy_flank_size, 
        entropy_threshold,
    );

    if snvs.len() >= too_many_snvs_threshold {  // TOO MANY, some kind of mismapping going on?
        get_snvs_meticulous(
            query_sequence, 
            ref_seq, 
            &query_coords,
            &ref_coords,
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

    let qc = query_coords.extract::<Vec<usize>>().unwrap();
    let rc = ref_coords.extract::<Vec<usize>>().unwrap();

    let snvs = get_snvs_simple(
        query_sequence.to_str().unwrap(), 
        ref_seq.to_str().unwrap(), 
        &qc,
        &rc,
        ref_coord_start, 
        tr_start_pos, 
        tr_end_pos, 
        too_many_snvs_threshold,
        entropy_flank_size, 
        entropy_threshold,
    );

    if snvs.len() >= too_many_snvs_threshold {  // TOO MANY, some kind of mismapping going on?
        get_snvs_meticulous(
            query_sequence.to_str().unwrap(), 
            ref_seq.to_str().unwrap(), 
            &qc,
            &rc,
            ref_coord_start, 
            tr_start_pos, 
            tr_end_pos, 
            contiguous_threshold, 
            max_snv_group_size, 
            entropy_flank_size, 
            entropy_threshold,
        ).into_py_dict(py)
    } else {
        snvs.into_py_dict(py)
    }
}

fn find_base_at_pos(
    query_sequence: &str, 
    q_coords: &Vec<usize>, 
    r_coords: &Vec<usize>, 
    t: usize,
    start_left: usize,
) -> (char, usize) {
    let (idx, found) = find_coord_idx_by_ref_pos(r_coords, t, start_left);

    if found {
        // Even if not in SNV set, it is not guaranteed to be a reference base, since
        // it's possible it was surrounded by too much other variation during the original
        // SNV getter algorithm.
        let qc = query_sequence.chars().nth(q_coords[idx]).unwrap();
        (qc, idx)
    } else {
        // Nothing found, so must have been a gap
        (SNV_GAP_CHAR, idx)
    }
}

pub fn calculate_useful_snvs(
    read_dict_extra: HashMap<&str, &PyDict>,
    read_q_coords: HashMap<&str, Vec<usize>>,
    read_r_coords: HashMap<&str, Vec<usize>>,
    read_snvs: HashMap<&str, HashMap<usize, char>>,
    locus_snvs: HashSet<usize>,
    min_allele_reads: usize,
) -> Vec<(usize, usize)> {
    // Mutates read_dict_extra - adds snv_bases keys to read entries

    let n_reads = read_dict_extra.len();

    let mut sorted_snvs: Vec<usize> = locus_snvs.into_iter().collect();
    sorted_snvs.sort();

    let mut snv_counters: HashMap<usize, HashMap<char, usize>> = 
        sorted_snvs
            .iter()
            .map(|&s| (s, HashMap::new()))
            .collect();

    for rn in read_dict_extra.keys() {
        let &read_dict_extra_for_read = read_dict_extra.get(rn).unwrap();
        let snvs = read_snvs.get(rn).unwrap();

        // Know this to not be None since we were passed only segments with non-None strings earlier
        let qs = read_dict_extra_for_read.get_item("_qs").unwrap().unwrap().extract::<&str>().unwrap();

        let q_coords = read_q_coords.get(rn).unwrap();
        let r_coords = read_r_coords.get(rn).unwrap();

        let segment_start = read_dict_extra_for_read.get_item("_ref_start")
            .unwrap().unwrap().extract::<usize>().unwrap();
        let segment_end = read_dict_extra_for_read.get_item("_ref_end")
            .unwrap().unwrap().extract::<usize>().unwrap();

        let mut snv_list: Vec<char> = Vec::new();
        let mut last_pair_idx: usize = 0;

        for &snv_pos in sorted_snvs.iter() {
            let mut base: char = SNV_OUT_OF_RANGE_CHAR;
            if segment_start <= snv_pos && snv_pos <= segment_end {
                if let Some(&bb) = snvs.get(&snv_pos) {
                    base = bb;
                } else {
                    // Binary search for base from correct pair
                    //  - We go in order, so we don't have to search left of the last pair index we tried.
                    let (sb, idx) = find_base_at_pos(qs, q_coords, r_coords, snv_pos, last_pair_idx);
                    // eprintln!("snv_pos {}, sb {} idx {}", snv_pos, sb, idx);
                    base = sb;
                    last_pair_idx = idx;
                }
            }

            // Otherwise, leave as out-of-range

            snv_list.push(base);

            if base != SNV_OUT_OF_RANGE_CHAR && base != SNV_GAP_CHAR {
                // Only count SNV bases where we've actually found the base for the read.
                let counter = snv_counters.get_mut(&snv_pos).unwrap();
                let count = counter.entry(base).or_insert_with(|| 0);
                *count += 1;
            }
        }

        // TODO: set snv_bases as tuple
        read_dict_extra_for_read.set_item("snv_bases", snv_list).unwrap();
    }

    // Enough reads to try for SNV based separation

    // require 2 alleles for the SNV, both with at least 1/5 of the reads, in order to differentiate alleles.
    // require over ~55% of the reads to have the SNV; otherwise it becomes too easy to occasionally get cases with
    // disjoint sets of SNVs.

    // TODO: parametrize proportion:
    let allele_read_threshold = cmp::max((n_reads as f32 / 5.0).round() as usize, min_allele_reads);
    let total_read_threshold = cmp::max((n_reads as f32 * 0.55).round() as usize, 5);  // TODO: parametrize

    // snv_counters is guaranteed by the previous inner loop to not have SNV_OUT_OF_RANGE_CHAR or SNV_GAP_CHAR

    let res: Vec<(usize, usize)> = sorted_snvs.iter().enumerate().filter_map(|(s_idx, s_pos)| {
        let snv_counter = snv_counters.get(s_pos).unwrap();

        let n_alleles_meeting_threshold: usize = snv_counter.values().map(|&v| (v >= allele_read_threshold) as usize).sum();
        let n_reads_with_this_snv_called: usize = snv_counter.values().sum();

        (n_alleles_meeting_threshold >= 2 && n_reads_with_this_snv_called >= total_read_threshold).then(|| (s_idx, sorted_snvs[s_idx]))
    }).collect();

    res
}
