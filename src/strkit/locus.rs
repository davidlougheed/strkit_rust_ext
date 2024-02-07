use pyo3::prelude::*;
use pyo3::types::{PyAny, PyDict, PyList};
use std::collections::{HashMap, HashSet};

use crate::cigar::get_aligned_pair_matches;
use crate::strkit::snvs::{get_read_snvs_rs, calculate_useful_snvs};
use crate::strkit::utils::find_coord_idx_by_ref_pos;

fn get_read_coords_from_matched_pairs(
    left_flank_coord: i32,
    left_coord: i32,
    right_coord: i32,
    right_flank_coord: i32,
    motif: &str,
    motif_size: i32,
    query_seq: &str,
    q_coords: &Vec<usize>,
    r_coords: &Vec<usize>,
) -> (i32, i32, i32, i32) {
    let mut left_flank_end: i32 = -1;
    let mut right_flank_start: i32 = -1;
    let mut right_flank_end: i32 = -1;

    let mut last_idx: i32 = -1;

    // Skip gaps on either side to find mapped flank indices

    // Binary search for left flank start ------------------------------------------------------------------------------
    
    let (mut lhs, found) = find_coord_idx_by_ref_pos(r_coords, left_flank_coord as usize, 0);

    // lhs now contains the index for the closest starting coordinate to left_flank_coord

    if !found && (lhs == 0 || lhs == r_coords.len()) {
        // Completely out of bounds; either right at the start or inserting after the end
        return (-1, -1, -1, -1);
    }

    if !found {
        // Choose pair to the left of where we'd insert the pair to maintain sorted order, since we want the closest
        // starting coordinate to left_flank_coord which gives us enough flanking material.
        lhs -=1;
    }

    let left_flank_start: i32 = q_coords[lhs] as i32;

    // -----------------------------------------------------------------------------------------------------------------

    for i in lhs+1..q_coords.len() {
        let query_coord = q_coords[i] as i32;
        let ref_coord = r_coords[i] as i32;

        // Skip gaps on either side to find mapped flank indices

        if ref_coord < left_coord {
            // Coordinate here is exclusive - we don't want to include a gap between the flanking region and
            // the STR; if we include the left-most base of the STR, we will have a giant flanking region which
            // will include part of the tandem repeat itself.
            left_flank_end = query_coord + 1;  // Add 1 to make it exclusive
        } else if ref_coord >= right_coord && (
            // Reached end of TR region and haven't set end of TR region yet, or there was an indel with the motif
            // in it right after we finished due to a subtle mis-alignment - this can be seen in the HTT alignments
            // in bc1018
            // TODO: do the same thing for the left side
            right_flank_start == -1 ||
            (
                query_coord - last_idx >= motif_size && 
                (ref_coord - right_coord <= motif_size * 2) &&
                (query_seq[(last_idx as usize)..(query_coord as usize)].matches(motif).count() as f64 / ((query_coord - last_idx) / motif_size) as f64) >= 0.5
            )
        ) {
            right_flank_start = query_coord;
        } else if ref_coord >= right_flank_coord {
            right_flank_end = query_coord;
            break;
        }

        last_idx = query_coord;
    }

    (left_flank_start, left_flank_end, right_flank_start, right_flank_end)
}

#[pyfunction]
pub fn get_pairs_and_tr_read_coords(
    cigar: &PyList,
    segment_start: usize,
    left_flank_coord: i32,
    left_coord: i32,
    right_coord: i32,
    right_flank_coord: i32,
    motif: &str,
    motif_size: i32,
    query_seq: &str,
) -> (Option<(Vec<usize>, Vec<usize>)>, i32, i32, i32, i32) {
    let (q_coords, r_coords) = get_aligned_pair_matches(cigar, 0, segment_start);
    let (left_flank_start, left_flank_end, right_flank_start, right_flank_end) = get_read_coords_from_matched_pairs(
        left_flank_coord, 
        left_coord, 
        right_coord, 
        right_flank_coord, 
        motif,
        motif_size, 
        query_seq, 
        &q_coords, 
        &r_coords
    );

    if left_flank_start == -1 || left_flank_end == -1 || right_flank_start == -1 || right_flank_end == -1 {
        // Avoid passing large vectors over Python-Rust boundary, return a None instead
        (None, left_flank_start, left_flank_end, right_flank_start, right_flank_end)
    } else {
        (Some((q_coords, r_coords)), left_flank_start, left_flank_end, right_flank_start, right_flank_end)
    }
}

#[pyfunction]
pub fn process_read_snvs_for_locus_and_calculate_useful_snvs(
    left_coord_adj: usize,
    right_coord_adj: usize,
    left_most_coord: usize,
    ref_cache: &str,
    read_dict_extra: HashMap<&str, &PyDict>,
    read_q_coords: HashMap<&str, Vec<usize>>,
    read_r_coords: HashMap<&str, Vec<usize>>,
    candidate_snvs_dict: &PyDict,
    min_allele_reads: usize,
    significant_clip_snv_take_in: usize,
    only_known_snvs: bool,
    logger: &PyAny,
    locus_log_str: &str,
) -> Vec<(usize, usize)> {
    // Loop through a second time if we are using SNVs. We do a second loop rather than just using the first loop
    // in order to have collected the edges of the reference sequence we can cache for faster SNV calculation.

    // Mutates: read_dict_extra

    let mut locus_snvs: HashSet<usize> = HashSet::new();
    let mut read_snvs: HashMap<&str, HashMap<usize, char>> = HashMap::new();

    for rn in read_dict_extra.keys() {
        let &read_dict_extra_for_read = read_dict_extra.get(rn).unwrap();

        let scl = read_dict_extra_for_read.get_item("sig_clip_left").unwrap().unwrap().extract::<usize>().unwrap();
        let scr = read_dict_extra_for_read.get_item("sig_clip_right").unwrap().unwrap().extract::<usize>().unwrap();

        if scl > 0 || scr > 0 {
            logger.call_method(
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

        let query_coords = read_q_coords.get(rn).unwrap();

        let twox_takein = significant_clip_snv_take_in * 2;
        if query_coords.len() < twox_takein {
            logger.call_method(
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
            continue;
        }

        let mut q_coords = query_coords.clone();
        let mut r_coords = read_r_coords.get(rn).unwrap().clone();

        if scl > 0 {
            q_coords.drain(0..scl);
            r_coords.drain(0..scl);
        }
        if scr > 0 {
            q_coords.drain((q_coords.len() - scr)..q_coords.len());
            r_coords.drain((r_coords.len() - scr)..r_coords.len());
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
            if !only_known_snvs || candidate_snvs_dict.contains(&p).unwrap() {
                locus_snvs.insert(p);
            }
        }

        read_snvs.insert(rn, snvs);
    }

    // --------------------------------------------------------------------------------------------

    calculate_useful_snvs(read_dict_extra, read_q_coords, read_r_coords, read_snvs, locus_snvs, min_allele_reads)
}