use pyo3::prelude::*;
use pyo3::types::{PyAny, PyDict};
use std::collections::{HashMap, HashSet};

use crate::strkit::snvs::{get_read_snvs_rs, calculate_useful_snvs};

#[pyfunction]
pub fn process_read_snvs_for_locus_and_calculate_useful_snvs(
    left_coord_adj: usize,
    right_coord_adj: usize,
    left_most_coord: usize,
    ref_cache: &str,
    read_dict_extra: HashMap<&str, &PyDict>,
    read_q_coords: HashMap<&str, Vec<usize>>,
    read_r_coords: HashMap<&str, Vec<usize>>,
    candidate_snvs_dict: HashMap<usize, &PyDict>,
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
            if !only_known_snvs || candidate_snvs_dict.contains_key(&p) {
                locus_snvs.insert(p);
            }
        }

        read_snvs.insert(rn, snvs);
    }

    // --------------------------------------------------------------------------------------------

    calculate_useful_snvs(read_dict_extra, read_q_coords, read_r_coords, read_snvs, locus_snvs, min_allele_reads)
}
