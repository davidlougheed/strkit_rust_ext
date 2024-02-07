use pyo3::prelude::*;
use pyo3::types::{IntoPyDict, PyAny, PyDict, PyList, PySet, PyString};
use std::collections::{HashMap, HashSet};
use std::cmp;

use crate::strkit::snvs::get_read_snvs_rs;

static SNV_OUT_OF_RANGE_CHAR: char = '-';
static SNV_GAP_CHAR: char = '_';

fn find_base_at_pos(
    query_sequence: &str, 
    q_coords: Vec<usize>, 
    r_coords: Vec<usize>, 
    t: usize,
    start_left: usize,
) -> (char, usize) {
    let idx = r_coords.partition_point(|&x| x <= t);
    let found = idx < r_coords.len() && r_coords[idx] == t;

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

fn calculate_useful_snvs(
    n_reads: usize,
    read_dict_items: HashMap<&str, &PyDict>,
    read_dict_extra: HashMap<&str, &PyDict>,
    read_match_pairs: HashMap<&str, (Vec<usize>, Vec<usize>)>,
    read_snvs: HashMap<&str, HashMap<usize, char>>,
    locus_snvs: HashSet<usize>,
    min_allele_reads: usize,
) -> Vec<(usize, usize)> {
    // Mutates read_dict_extra - adds snv_bases keys to read entries

    let mut sorted_snvs: Vec<usize> = locus_snvs.into_iter().collect();
    sorted_snvs.sort();

    let mut snv_counters: HashMap<usize, HashMap<char, usize>> = sorted_snvs.iter().map(|&s| (s, HashMap::new())).collect();

    for (rn, read) in read_dict_items {
        let &read_dict_extra_for_read = read_dict_extra.get(rn).unwrap();
        let snvs = read_snvs.get(rn).unwrap();

        // Know this to not be None since we were passed only segments with non-None strings earlier
        let qs = read_dict_extra_for_read.get_item("_qs").unwrap().unwrap().extract::<&str>().unwrap();

        let &(q_coords, r_coords) = read_match_pairs.get(rn).unwrap();

        let segment_start = read_dict_extra_for_read.get_item("_ref_start")
            .unwrap().unwrap().extract::<usize>().unwrap();
        let segment_end = read_dict_extra_for_read.get_item("_ref_end")
            .unwrap().unwrap().extract::<usize>().unwrap();

        let mut snv_list: Vec<char> = Vec::new();
        let mut last_pair_idx: usize = 0;

        for snv_pos in sorted_snvs {
            let mut base: char = SNV_OUT_OF_RANGE_CHAR;
            if segment_start <= snv_pos && snv_pos <= segment_end {
                if let Some(&bb) = snvs.get(&snv_pos) {
                    base = bb;
                } else {
                    // Binary search for base from correct pair
                    //  - We go in order, so we don't have to search left of the last pair index we tried.
                    let (sb, idx) = find_base_at_pos(qs, q_coords, r_coords, snv_pos, last_pair_idx);
                    last_pair_idx = idx;
                }
            }

            // Otherwise, leave as out-of-range

            snv_list.push(base);

            if base != SNV_OUT_OF_RANGE_CHAR && base != SNV_GAP_CHAR {
                // Only count SNV bases where we've actually found the base for the read.
                *snv_counters.get(&snv_pos).unwrap().entry(base).or_insert_with(|| 0) += 1;
            }
        }

        // TODO: set snv_bases as tuple
        read_dict_extra_for_read.set_item("snv_bases", snv_list);
    }

    // Enough reads to try for SNV based separation

    // require 2 alleles for the SNV, both with at least 1/5 of the reads, in order to differentiate alleles.
    // require over ~55% of the reads to have the SNV; otherwise it becomes too easy to occasionally get cases with
    // disjoint sets of SNVs.

    let mut useful_snvs: Vec<usize> = Vec::new();

    // TODO: parametrize proportion:
    let allele_read_threshold = cmp::max((n_reads as f32 / 5.0).round() as usize, min_allele_reads);
    let total_read_threshold = cmp::max((n_reads as f32 * 0.55).round() as usize, 5);  // TODO: parametrize

    // snv_counters is guaranteed by the previous inner loop to not have SNV_OUT_OF_RANGE_CHAR or SNV_GAP_CHAR

    // TODO

    useful_snvs.into_iter().map(|s_idx| (s_idx, sorted_snvs[s_idx])).collect()
}

#[pyfunction]
pub fn process_read_snvs_for_locus_and_calculate_useful_snvs(
    contig: &str,
    left_coord_adj: usize,
    right_coord_adj: usize,
    left_most_coord: usize,
    right_most_coord: usize,
    ref_fasta: &PyAny,
    read_dict_items: HashMap<&str, &PyDict>,
    read_dict_extra: HashMap<&str, &PyDict>,
    read_match_pairs: HashMap<&str, (Vec<usize>, Vec<usize>)>,
    candidate_snvs_dict: HashMap<usize, &PyDict>,
    significant_clip_snv_take_in: usize,
    only_known_snvs: bool,
    logger_: &PyAny,
    locus_log_str: &str,
) -> Vec<(usize, usize)> {
    // Loop through a second time if we are using SNVs. We do a second loop rather than just using the first loop
    // in order to have collected the edges of the reference sequence we can cache for faster SNV calculation.

    // Mutates: read_dict_extra

    let mut locus_snvs: HashSet<usize> = HashSet::new();
    let mut read_snvs: HashMap<&str, HashMap<usize, char>> = HashMap::new();

    let ref_cache = ref_fasta
        .call_method("fetch", (contig, left_most_coord, right_most_coord+1), None)
        .unwrap()
        .extract::<&str>()
        .unwrap()
        .to_ascii_uppercase()
        .as_str();

    for (rn, read) in read_dict_items {
        let &read_dict_extra_for_read = read_dict_extra.get(rn).unwrap();

        let scl = read_dict_extra_for_read.get_item("sig_clip_left").unwrap().unwrap().extract::<usize>().unwrap();
        let scr = read_dict_extra_for_read.get_item("sig_clip_right").unwrap().unwrap().extract::<usize>().unwrap();

        if scl == 0 || scr == 0 {
            logger_.call_method(
                "debug", 
                (
                    format!(
                        "{} - {} has significant clipping; trimming pairs by {} bp per side for SNV-finding",
                        locus_log_str,
                        rn,
                        significant_clip_snv_take_in
                    ),
                ),
                None
            ).unwrap();
        }

        let query_coords = read_match_pairs.get(rn).unwrap().0;

        let twox_takein = significant_clip_snv_take_in * 2;
        if query_coords.len() < twox_takein {
            logger_.call_method(
                "warning", 
                (
                    format!(
                        "{} - skipping SNV calculation for '{}' (<{} pairs)",
                        locus_log_str,
                        rn,
                        twox_takein
                    ),
                ),
                None
            ).unwrap();
        }

        let mut q_coords = query_coords.clone();
        let mut r_coords = read_match_pairs.get(rn).unwrap().0;

        if scl != 0 {
            q_coords.drain(0..scl);
        }
        if scr != 0 {
            r_coords.drain(0..scr);
        }

        let query_sequence = read_dict_extra_for_read.get_item("_qs").unwrap().unwrap().extract::<&str>().unwrap();

        let snvs = get_read_snvs_rs(
            query_sequence,
            ref_cache,
            q_coords,
            r_coords,
            left_most_coord,
            left_coord_adj,
            right_coord_adj,
            // below: magic values for skipping false positives / weird 'SNVs' that aren't helpful
            5,  // contiguous_threshold
            5,  // max_snv_group_size
            20,  // too_many_snvs_threshold
            10,  // entropy_flank_size
            1.8  // entropy_threshold
        );

        for &p in snvs.keys() {
            if !only_known_snvs || candidate_snvs_dict.contains_key(&p) {
                locus_snvs.insert(p);
            }
        }

        read_snvs.insert(rn, snvs);
    }

    // --------------------------------------------------------------------------------------------

    locus_snvs
}
