use bincode;
use numpy::{PyArray1, PyArray2, PyArrayMethods};
use pyo3::intern;
use pyo3::prelude::*;
use pyo3::types::{PyAny, PyBytes, PyDict, PyString};
use serde::{Deserialize, Serialize};
use std::collections::{HashMap, HashSet};

use crate::aligned_coords::STRkitAlignedCoords;
use crate::strkit::cigar::get_aligned_pair_matches_rs;
use crate::strkit::snvs::{CandidateSNVs, calculate_useful_snvs, get_read_snvs};
use crate::strkit::utils::find_coord_idx_by_ref_pos;

use super::snvs::UsefulSNVsParams;


#[derive(Clone, Deserialize, Serialize)]
#[pyclass(module = "strkit_rust_ext")]
pub struct STRkitLocus {
    #[pyo3(get)]
    pub t_idx: usize,
    #[pyo3(get)]
    pub locus_id: String,

    #[pyo3(get)]
    pub contig: String,

    #[pyo3(get)]
    pub left_coord: i32,
    #[pyo3(get)]
    pub left_flank_coord: i32,
    #[pyo3(get)]
    pub right_coord: i32,
    #[pyo3(get)]
    pub right_flank_coord: i32,
    #[pyo3(get)]
    pub ref_size: i32, // reference size, in terms of coordinates

    #[pyo3(get)]
    pub motif: String,
    #[pyo3(get)]
    pub motif_size: usize,

    #[pyo3(get)]
    pub n_alleles: usize,

    _flank_size: i32,
    _log_str: String,
}

#[pymethods]
impl STRkitLocus {
    #[new]
    fn py_new(
        t_idx: usize,
        locus_id: &str,
        contig: &str,
        left_coord: i32,
        right_coord: i32,
        motif: &str,
        n_alleles: usize,
        flank_size: i32,
    ) -> PyResult<Self> {
        let log_str = format!(
            "locus {} (id={}): {}:{}-{} [{}]", t_idx, locus_id, contig, left_coord, right_coord, motif
        );

        Ok(
            STRkitLocus {
                t_idx,
                locus_id: locus_id.to_string(),

                contig: contig.to_string(),

                left_coord,
                left_flank_coord: left_coord - flank_size,
                right_coord,
                right_flank_coord: right_coord + flank_size,
                ref_size: right_coord - left_coord,

                motif: motif.to_string(),
                motif_size: motif.len(),

                n_alleles,

                _flank_size: flank_size,
                _log_str: log_str,
            }
        )
    }

    pub fn log_str(&self) -> &str { &self._log_str }

    fn to_dict<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        let res = PyDict::new(py);
        res.set_item("locus_index", self.t_idx).unwrap();
        res.set_item("locus_id", self.locus_id.clone()).unwrap();
        res.set_item("contig", self.contig.clone()).unwrap();
        res.set_item("start", self.left_coord).unwrap();
        res.set_item("end", self.right_coord).unwrap();
        res.set_item("motif", self.motif.clone()).unwrap();
        Ok(res)
    }

    fn with_ref_data(
        &self,
        left_coord_adj: i32,
        right_coord_adj: i32,
        ref_contig: String,
        ref_cn: i32,
        ref_seq: String,
        ref_left_flank_seq: String,
        ref_right_flank_seq: String,
        ref_total_seq: String,
        ref_time: f64,
    ) -> PyResult<STRkitLocusWithRefData> {
        // Hacky Python-boundary-crossing version of a builder pattern
        Ok(
            STRkitLocusWithRefData {
                locus_def: self.clone(), // If this becomes Rust-only, this can be a move instead.
                left_coord_adj,
                right_coord_adj,
                ref_contig,
                ref_cn,
                ref_seq,
                ref_left_flank_seq,
                ref_right_flank_seq,
                ref_total_seq,
                ref_time,
            }
        )
    }

    pub fn __setstate__(&mut self, state: Bound<'_, PyBytes>) -> PyResult<()> {
        (*self, _) = bincode::serde::decode_from_slice(state.as_bytes(), bincode::config::standard()).unwrap();
        Ok(())
    }

    pub fn __getstate__<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyBytes>> {
        Ok(PyBytes::new(py, &bincode::serde::encode_to_vec(self, bincode::config::standard()).unwrap()))
    }

    pub fn __getnewargs__(&self) -> PyResult<(usize, String, String, i32, i32, String, usize, i32)> {
        Ok((
            self.t_idx,
            self.locus_id.clone(),
            self.contig.clone(),
            self.left_coord,
            self.right_coord,
            self.motif.clone(),
            self.n_alleles,
            self._flank_size,
        ))
    }
}

#[derive(Clone)]
#[pyclass(frozen)]
pub struct STRkitLocusWithRefData {
    pub locus_def: STRkitLocus,

    #[pyo3(get)]
    pub left_coord_adj: i32,
    #[pyo3(get)]
    pub right_coord_adj: i32,

    #[pyo3(get)]
    pub ref_contig: String,
    #[pyo3(get)]
    pub ref_cn: i32,
    #[pyo3(get)]
    pub ref_seq: String,
    #[pyo3(get)]
    pub ref_left_flank_seq: String,
    #[pyo3(get)]
    pub ref_right_flank_seq: String,
    #[pyo3(get)]
    pub ref_total_seq: String,
    #[pyo3(get)]
    pub ref_time: f64,

    // This struct must be built with the STRkitLocus.with_ref_data builder function above.
}

fn _get_read_coords_from_matched_pairs(
    locus_with_ref_data: &STRkitLocusWithRefData,
    query_seq: &str,
    aligned_coords: &STRkitAlignedCoords,
) -> (i32, i32, i32, i32) {
    // Skip gaps on either side to find mapped flank indices

    // Binary search for left flank start ------------------------------------------------------------------------------

    let (mut lhs, found) = find_coord_idx_by_ref_pos(
        &aligned_coords, locus_with_ref_data.locus_def.left_flank_coord as usize, 0
    );

    // lhs now contains the index for the closest starting coordinate to left_flank_coord

    if !found && (lhs == 0 || lhs == aligned_coords.ref_coords.len()) {
        // Completely out of bounds; either right at the start or inserting after the end
        return (-1, -1, -1, -1);
    }

    if !found {
        // Choose pair to the left of where we'd insert the pair to maintain sorted order, since we want the closest
        // starting coordinate to left_flank_coord which gives us enough flanking material.
        lhs -= 1;
    }

    let left_flank_start: i32 = aligned_coords.query_coords[lhs] as i32;

    // -----------------------------------------------------------------------------------------------------------------

    // Binary search for left flank end (may "fail" if we don't have a pair for left_coord - 1 (i.e., there's a gap in
    // the reference to the direct left of left_coord), in which case we can do it the slow way in O(n) time).

    let mut left_flank_end: i32 = -1;

    let mut loop_start = lhs + 1;

    let (lhs_end, lhs_end_found) = find_coord_idx_by_ref_pos(
        &aligned_coords,
        (locus_with_ref_data.left_coord_adj - 1) as usize,
        loop_start,
    );
    if lhs_end_found {
        left_flank_end = aligned_coords.query_coords[lhs_end] as i32 + 1;
        loop_start = lhs_end + 1;
    } else {
        // eprintln!("lhs_end q_coord={}, found={}", q_coords[lhs_end] + 1, lhs_end_found);
    }

    // Binary search for right flank end -------------------------------------------------------------------------------

    let motif_size = locus_with_ref_data.locus_def.motif_size as i32;

    let mut right_flank_start: i32 = -1;
    let mut right_flank_end: i32 = -1;

    let mut last_idx: i32 = -1;

    for i in loop_start..aligned_coords.query_coords.len() {
        let query_coord = aligned_coords.query_coords[i] as i32;
        let ref_coord = aligned_coords.ref_coords[i] as i32;

        // Skip gaps on either side to find mapped flank indices

        if ref_coord < locus_with_ref_data.left_coord_adj {
            // Coordinate here is exclusive - we don't want to include a gap between the flanking region and
            // the STR; if we include the left-most base of the STR, we will have a giant flanking region which
            // will include part of the tandem repeat itself.
            left_flank_end = query_coord + 1; // Add 1 to make it exclusive
        } else if ref_coord >= locus_with_ref_data.right_coord_adj
            && (
                // Reached end of TR region and haven't set end of TR region yet, or there was an indel with the motif
                // in it right after we finished due to a subtle mis-alignment - this can be seen in the HTT alignments
                // in bc1018
                // TODO: do the same thing for the left side
                right_flank_start == -1
                    || (query_coord - last_idx >= motif_size
                        && (ref_coord - locus_with_ref_data.right_coord_adj <= motif_size * 2)
                        && (query_seq[(last_idx as usize)..(query_coord as usize)]
                            .matches(&locus_with_ref_data.locus_def.motif)
                            .count() as f64
                            / ((query_coord - last_idx) / motif_size) as f64)
                            >= 0.5)
            )
        {
            right_flank_start = query_coord;
        } else if ref_coord >= locus_with_ref_data.locus_def.right_flank_coord {
            right_flank_end = query_coord;
            break;
        }

        last_idx = query_coord;
    }

    (left_flank_start, left_flank_end, right_flank_start, right_flank_end)
}

#[pyfunction]
pub fn get_read_coords_from_matched_pairs(
    locus_with_ref_data: &STRkitLocusWithRefData,
    query_seq: &str,
    aligned_coords: &Bound<'_, STRkitAlignedCoords>,
) -> (i32, i32, i32, i32) {
    _get_read_coords_from_matched_pairs(locus_with_ref_data, query_seq, &aligned_coords.borrow())
}

#[pyfunction]
pub fn get_pairs_and_tr_read_coords<'py>(
    py: Python<'py>,
    cigar: &Bound<'py, PyArray2<u32>>,
    segment_start: u64,
    locus_with_ref_data: &Bound<'py, STRkitLocusWithRefData>,
    query_seq: &str,
) -> (Option<Py<STRkitAlignedCoords>>, i32, i32, i32, i32) {
    let aligned_coords = get_aligned_pair_matches_rs(cigar, 0, segment_start);
    let (left_flank_start, left_flank_end, right_flank_start, right_flank_end) =
        _get_read_coords_from_matched_pairs(&locus_with_ref_data.borrow(), query_seq, &aligned_coords);

    if left_flank_start == -1 || left_flank_end == -1 || right_flank_start == -1 || right_flank_end == -1 {
        // Avoid converting to Python objects / passing over Python-Rust boundary, return a None instead
        (None, left_flank_start, left_flank_end, right_flank_start, right_flank_end)
    } else {
        (
            Some(Py::new(py, aligned_coords).unwrap()),
            left_flank_start,
            left_flank_end,
            right_flank_start,
            right_flank_end,
        )
    }
}

#[pyfunction]
pub fn process_read_snvs_for_locus_and_calculate_useful_snvs(
    py: Python<'_>,
    left_coord_adj: usize,
    right_coord_adj: usize,
    left_most_coord: usize,
    ref_cache: &str,
    read_dict_extra: Bound<PyDict>,
    read_aligned_coords: &Bound<PyDict>,
    candidate_snvs: &Bound<'_, CandidateSNVs>,
    min_allele_reads: usize,
    significant_clip_snv_take_in: usize,
    only_known_snvs: bool,
    logger: Bound<PyAny>,
    locus_log_str: &str,
) -> Result<Vec<(usize, usize)>, PyErr> {
    // Loop through a second time if we are using SNVs. We do a second loop rather than just using the first loop
    // in order to have collected the edges of the reference sequence we can cache for faster SNV calculation.

    // Mutates: read_dict_extra

    let mut locus_snvs: HashSet<usize> = HashSet::new();
    let mut read_snvs: HashMap<String, HashMap<usize, (char, u8)>> = HashMap::new();

    // below: magic values for skipping false positives / weird 'SNVs' that aren't helpful
    let useful_snvs_params = UsefulSNVsParams {
        contiguous_threshold: 5,
        max_snv_group_size: 5,
        too_many_snvs_threshold: 20,
        entropy_flank_size: 10,
        entropy_threshold: 1.8,
    };

    let candidate_snvs_b = candidate_snvs.borrow();

    for rn in read_dict_extra.keys().into_iter().map(|x| x.downcast_into::<PyString>().unwrap()) {
        let read_dict_extra_for_read = read_dict_extra.get_item(&rn).unwrap().unwrap().downcast_into::<PyDict>().unwrap();

        let scl = read_dict_extra_for_read.get_item(intern!(py, "sig_clip_left")).unwrap().unwrap().extract::<usize>().unwrap();
        let scr = read_dict_extra_for_read.get_item(intern!(py, "sig_clip_right")).unwrap().unwrap().extract::<usize>().unwrap();

        if scl > 0 || scr > 0 {
            logger.call_method1(
                intern!(py, "debug"),
                (
                    intern!(py, "%s - %s has significant clipping; trimming pairs by %d bp per side for SNV-finding"),
                    locus_log_str,
                    &rn,
                    significant_clip_snv_take_in,
                ),
            )?;
        }

        let aligned_coords =
            read_aligned_coords
                .get_item(&rn)?
                .unwrap()
                .downcast_into::<STRkitAlignedCoords>()?
                .borrow();
        let coords_len = aligned_coords.query_coords.len();

        let twox_takein = significant_clip_snv_take_in * 2;
        if coords_len < twox_takein {
            logger.call_method1(
                intern!(py, "warning"),
                (
                    intern!(py, "%s - skipping SNV calculation for '%s' (<%d pairs)"),
                    locus_log_str,
                    &rn,
                    twox_takein,
                ),
            )?;
            continue;
        }

        let qsi = read_dict_extra_for_read.get_item(intern!(py, "_qs"))?.unwrap();
        let query_sequence = qsi.extract::<&str>()?;
        let fqqs_i1 = read_dict_extra_for_read.get_item(intern!(py, "_fqqs"))?.unwrap();
        let fqqs_i2 = fqqs_i1.downcast::<PyArray1<u8>>()?.readonly();
        let fqqs = fqqs_i2.as_slice()?;

        let snvs = get_read_snvs(
            query_sequence,
            fqqs,
            ref_cache,
            // aligned coords without clipping at ends
            &aligned_coords.query_coords[scl..(coords_len - scr)],
            &aligned_coords.ref_coords[scl..(coords_len - scr)],
            left_most_coord,
            left_coord_adj,
            right_coord_adj,
            // below: magic values for skipping false positives / weird 'SNVs' that aren't helpful
            &useful_snvs_params,
        );

        // Add all locus SNVs that are known candidates (if required)
        locus_snvs.extend(snvs.keys().filter(|p| !only_known_snvs || candidate_snvs_b.snvs.contains_key(p)));
        // Add read SNVs
        read_snvs.insert(rn.to_string(), snvs);
    }

    // --------------------------------------------------------------------------------------------

    calculate_useful_snvs(py, read_dict_extra, read_aligned_coords, read_snvs, locus_snvs, min_allele_reads)
}
