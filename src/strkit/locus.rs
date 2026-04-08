use bincode;
use pyo3::intern;
use pyo3::prelude::*;
use pyo3::types::{PyAny, PyBytes, PyDict, PyString};
use serde::{Deserialize, Serialize};
use std::collections::{HashMap, HashSet};
use std::hash::{DefaultHasher, Hash, Hasher};

use crate::aligned_coords::STRkitAlignedCoords;
use crate::reads::STRkitAlignedSegment;
use crate::reads::STRkitLocusBlockSegments;
use crate::reads::STRkitSegmentAlignmentDataForLocus;
use crate::snvs::{CandidateSNVs, GetReadSNVs, UsefulSNVsParams, calculate_useful_snvs};
use crate::coords::{QueryCoord, RefCoord};


#[derive(Clone, Deserialize, Hash, Serialize)]
#[pyclass(from_py_object, module = "strkit_rust_ext")]
pub struct STRkitLocus {
    #[pyo3(get)]
    pub t_idx: usize,
    #[pyo3(get)]
    pub locus_id: String,

    #[pyo3(get)]
    pub contig: String,

    #[pyo3(get)]
    pub left_coord: RefCoord,
    #[pyo3(get)]
    pub left_flank_coord: RefCoord,
    #[pyo3(get)]
    pub right_coord: RefCoord,
    #[pyo3(get)]
    pub right_flank_coord: RefCoord,
    #[pyo3(get)]
    pub ref_size: u64, // reference size, in terms of coordinates

    #[pyo3(get)]
    pub motif: String,
    #[pyo3(get)]
    pub motif_size: u16,

    #[pyo3(get)]
    pub n_alleles: usize,

    #[pyo3(get)]
    pub flank_size: u64,

    #[pyo3(get)]
    pub annotations: Vec<String>,

    _log_str: String,
}

#[pymethods]
impl STRkitLocus {
    #[new]
    fn py_new(
        t_idx: usize,
        locus_id: &str,
        contig: &str,
        left_coord: RefCoord,
        right_coord: RefCoord,
        motif: &str,
        n_alleles: usize,
        flank_size: u64,
        annotations: Vec<String>,
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
                motif_size: motif.len() as u16,

                n_alleles,

                flank_size,

                annotations,

                _log_str: log_str,
            }
        )
    }

    pub fn log_str(&self) -> &str { &self._log_str }

    fn to_dict<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        let res = PyDict::new(py);
        res.set_item("locus_index", self.t_idx)?;
        res.set_item("locus_id", self.locus_id.clone())?;
        res.set_item("contig", self.contig.clone())?;
        res.set_item("start", self.left_coord)?;
        res.set_item("end", self.right_coord)?;
        res.set_item("motif", self.motif.clone())?;
        res.set_item("annotations", self.annotations.clone())?;
        Ok(res)
    }

    fn __hash__(&self) -> u64 {
        // Implement hash by hand, as we cannot have a frozen PyO3 class here due to the pickling stuff.
        let mut hasher = DefaultHasher::new();
        self.hash(&mut hasher);
        hasher.finish()
    }

    fn __repr__(&self) -> String {
        let repr = format!(
            "<STRkitLocus t_idx={} locus_id={} contig={} left_coord={} left_flank_coord={} right_coord={} \
            right_flank_coord={} ref_size={} motif={} motif_size={} n_alleles={} flank_size={} _log_str={}>",
            &self.t_idx,
            &self.locus_id,
            &self.contig,
            &self.left_coord,
            &self.left_flank_coord,
            &self.right_coord,
            &self.right_flank_coord,
            &self.ref_size,
            &self.motif,
            &self.motif_size,
            &self.n_alleles,
            &self.flank_size,
            &self._log_str,
        );
        repr
    }

    fn with_ref_data(
        &self,
        left_coord_adj: u64,
        right_coord_adj: u64,
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

    // --- below are functions which make this class pickle-able ---

    pub fn __setstate__(&mut self, state: Bound<'_, PyBytes>) -> PyResult<()> {
        // TODO: replace unwrap with actual PyResult error
        (*self, _) = bincode::serde::decode_from_slice(state.as_bytes(), bincode::config::standard()).unwrap();
        Ok(())
    }

    pub fn __getstate__<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyBytes>> {
        // TODO: replace unwrap with actual PyResult error
        Ok(PyBytes::new(py, &bincode::serde::encode_to_vec(self, bincode::config::standard()).unwrap()))
    }

    pub fn __getnewargs__(&self) -> PyResult<(usize, String, String, u64, u64, String, usize, u64, Vec<String>)> {
        Ok((
            self.t_idx,
            self.locus_id.clone(),
            self.contig.clone(),
            self.left_coord,
            self.right_coord,
            self.motif.clone(),
            self.n_alleles,
            self.flank_size,
            self.annotations.clone(),
        ))
    }
}

#[derive(Clone)]
#[pyclass(skip_from_py_object, frozen)]
pub struct STRkitLocusWithRefData {
    #[pyo3(get)]
    pub locus_def: STRkitLocus,

    #[pyo3(get)]
    pub left_coord_adj: RefCoord,
    #[pyo3(get)]
    pub right_coord_adj: RefCoord,

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

#[pymethods]
impl STRkitLocusWithRefData {
    #[getter]
    fn left_flank_coord(&self) -> RefCoord {
        self.locus_def.left_flank_coord
    }

    #[getter]
    fn right_flank_coord(&self) -> RefCoord {
        self.locus_def.right_flank_coord
    }

    fn log_str(&self) -> &str {
        self.locus_def.log_str()
    }
}


#[pyclass]
pub struct STRkitLocusBlockIter {
    inner: std::vec::IntoIter<STRkitLocus>,
}

#[pymethods]
impl STRkitLocusBlockIter {
    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    fn __next__(mut slf: PyRefMut<'_, Self>) -> Option<STRkitLocus> {
        slf.inner.next()
    }
}

#[derive(Clone, Serialize, Deserialize)]
#[pyclass(from_py_object, module = "strkit_rust_ext")]
pub struct STRkitLocusBlock {
    pub loci: Vec<STRkitLocus>,
    pub left: RefCoord, // Left-most coordinate of all loci
    pub right: RefCoord, // Right-most coordinate of all loci
    pub log_str: String,
}

#[pymethods]
impl STRkitLocusBlock {
    #[new]
    fn py_new(loci: Vec<STRkitLocus>, left: u64, right: u64) -> PyResult<Self> {
        let log_str = format!("[block {}:{}-{}]", loci[0].contig, left, right);
        Ok(STRkitLocusBlock { loci, left, right, log_str })
    }

    #[getter]
    fn contig(&self) -> String {
        self.loci[0].contig.clone()
    }

    pub fn __len__(&self) -> usize {
        self.loci.len()
    }

    fn __iter__(slf: PyRef<'_, Self>) -> PyResult<Py<STRkitLocusBlockIter>> {
        let iter = STRkitLocusBlockIter { inner: slf.loci.clone().into_iter() };
        Py::new(slf.py(), iter)
    }

    // --- below are functions which make this class pickle-able ---

    pub fn __setstate__(&mut self, state: Bound<'_, PyBytes>) -> PyResult<()> {
        // TODO: replace unwrap with actual PyResult error
        (*self, _) = bincode::serde::decode_from_slice(state.as_bytes(), bincode::config::standard()).unwrap();
        Ok(())
    }

    pub fn __getstate__<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyBytes>> {
        // TODO: replace unwrap with actual PyResult error
        Ok(PyBytes::new(py, &bincode::serde::encode_to_vec(self, bincode::config::standard()).unwrap()))
    }

    pub fn __getnewargs__(&self) -> PyResult<(Vec<STRkitLocus>, u64, u64)> {
        Ok((
            self.loci.clone(), // TODO: perf hack with None / empty vec?
            self.left,
            self.right,
        ))
    }
}


#[pyclass]
pub struct LocusReadCoords {
    #[pyo3(get)]
    pub left_flank_start: Option<QueryCoord>,
    #[pyo3(get)]
    pub left_flank_end: Option<QueryCoord>,
    #[pyo3(get)]
    pub right_flank_start: Option<QueryCoord>,
    #[pyo3(get)]
    pub right_flank_end: Option<QueryCoord>,
    #[pyo3(get)]
    pub full_left_flank: bool,
    #[pyo3(get)]
    pub full_right_flank: bool,
}

#[pymethods]
impl LocusReadCoords {
    #[staticmethod]
    pub fn new_all_incomplete() -> Self {
        return LocusReadCoords {
            left_flank_start: None,
            left_flank_end: None,
            right_flank_start: None,
            right_flank_end: None,
            full_left_flank: false,
            full_right_flank: false,
        };
    }

    pub fn is_incomplete(&self) -> bool {
        return self.left_flank_start.is_none()
            || self.left_flank_end.is_none()
            || self.right_flank_start.is_none()
            || self.right_flank_end.is_none();
    }

    pub fn __repr__(&self) -> String {
        format!(
            "<LocusReadCoords {:?} {:?} {:?} {:?} {} {}>",
            self.left_flank_start,
            self.left_flank_end,
            self.right_flank_start,
            self.right_flank_end,
            self.full_left_flank,
            self.full_right_flank,
        )
    }
}


/// Given some locus definition + reference data, a read query sequence, and a set of aligned coordinates between the
/// query and the reference, we want to find the query coordinates that encompass the locus + the flanking region.
/// To do this, we need to examine the aligned coordinates and use the locus definition to find flank/TR boundaries.
///
/// The aligned coordinates will be either from the original algined read (query) directly or a realignment of the read.
///
/// If allow_only_one_full_flank is set to True, we only require one full flanking sequence rather than both. This can
/// help recover reads, which is especially valuable for expansion discovery and/or at low coverage levels.
fn _get_read_coords_from_matched_pairs(
    locus_with_ref_data: &STRkitLocusWithRefData,
    query_seq: &str,
    aligned_coords: PyRef<'_, STRkitAlignedCoords>,
    vcf_anchor_size: u64,
    allow_only_one_full_flank: bool,
) -> LocusReadCoords {
    // Skip gaps on either side to find mapped flank indices

    // Binary search for left flank start ------------------------------------------------------------------------------

    let (mut lhs, mut full_left_flank) = aligned_coords.find_coord_idx_by_ref_pos(
        locus_with_ref_data.locus_def.left_flank_coord, 0
    );

    // lhs now contains the index for the closest starting coordinate to left_flank_coord (or is out-of-bounds)

    if !full_left_flank {
        if lhs == 0 || lhs == aligned_coords.ref_coords.len() {
            if !allow_only_one_full_flank {
                // Completely out of bounds (either right at the start or inserting after the end) and we require both
                // flanking sequences in their entirety.
                return LocusReadCoords::new_all_incomplete();
            }
        } else {
            // Choose pair to the left of where we'd insert the pair to maintain sorted order, since we want the closest
            // starting coordinate to left_flank_coord which gives us enough flanking material. This won't be the same
            // as the original target coordinate, but it'll encompass a full flank worth of sequence.
            lhs -= 1;
            full_left_flank = true;
        }
    }

    // If we're allowed to have a partial flank on one side, we can keep going even if we don't have the full left
    // flank. In that case, we will require a full right flank later on.
    let mut left_flank_start: QueryCoord = if full_left_flank { aligned_coords.query_coords[lhs] } else { 0 };

    // -----------------------------------------------------------------------------------------------------------------

    // Binary search for left flank end (may "fail" if we don't have a pair for left_coord - 1 (i.e., there's a gap in
    // the reference to the direct left of left_coord), in which case we can do it the slow way in O(n) time).

    let mut left_flank_end: Option<QueryCoord> = None;

    let mut loop_start = if full_left_flank { lhs + 1 } else { 0 };

    let (lhs_end, lhs_end_found) = aligned_coords.find_coord_idx_by_ref_pos(
        locus_with_ref_data.left_coord_adj - 1,
        loop_start,
    );
    if lhs_end_found {
        left_flank_end = Some(aligned_coords.query_coords[lhs_end] + 1);
        loop_start = lhs_end + 1;
    } else {
        // eprintln!("lhs_end q_coord={}, found={}", q_coords[lhs_end] + 1, lhs_end_found);
    }

    // Linear search for left flank end (if not found via binary search) and right flank start/end ---------------------

    let motif_size = locus_with_ref_data.locus_def.motif_size as u64;

    let mut right_flank_start: Option<QueryCoord> = None;
    let mut right_flank_end: Option<QueryCoord> = None;
    let mut full_right_flank = false;

    let mut last_idx: QueryCoord = 0;

    for i in loop_start..aligned_coords.query_coords.len() {
        let query_coord = aligned_coords.query_coords[i];
        let ref_coord = aligned_coords.ref_coords[i];

        // Skip gaps on either side to find mapped flank indices

        if ref_coord < locus_with_ref_data.left_coord_adj {
            // Coordinate here is exclusive - we don't want to include a gap between the flanking region and
            // the STR; if we include the left-most base of the STR, we will have a giant flanking region which
            // will include part of the tandem repeat itself.
            let lfe = query_coord + 1; // Add 1 to make it exclusive

            // Even if we're allowed to have a small flank (allow_only_one_full_flank), we still require a flanking
            // region large enough to anchor the sequence in the VCF output.
            if lfe - left_flank_start < vcf_anchor_size {
                return LocusReadCoords::new_all_incomplete();
            }

            left_flank_end = Some(lfe);
        } else if ref_coord >= locus_with_ref_data.right_coord_adj
            && (
                // Reached end of TR region and haven't set end of TR region yet, or there was an indel with the motif
                // in it right after we finished due to a subtle mis-alignment - this can be seen in the HTT alignments
                // in bc1018
                // TODO: do the same thing for the left side
                right_flank_start.is_none()
                    || (last_idx > 0
                        && query_coord - last_idx >= motif_size as u64
                        && (ref_coord - locus_with_ref_data.right_coord_adj <= motif_size * 2)
                        && (query_seq[(last_idx as usize)..(query_coord as usize)]
                            .matches(&locus_with_ref_data.locus_def.motif)
                            .count() as f64
                            / ((query_coord - last_idx) / motif_size) as f64)
                            >= 0.5)
            )
        {
            if left_flank_end.is_none() {
                // no left flank at all, early return error
                // TODO: in future: partial supporting read instead
                return LocusReadCoords::new_all_incomplete();
            }
            right_flank_start = Some(query_coord);
        } else if ref_coord >= locus_with_ref_data.locus_def.right_flank_coord || (
            allow_only_one_full_flank && full_left_flank && !right_flank_start.is_none()
        ) {
            right_flank_end = Some(query_coord);
            if ref_coord >= locus_with_ref_data.locus_def.right_flank_coord {
                full_right_flank = true;
                break;
            }
        }

        last_idx = query_coord;
    }

    // -----------------------------------------------------------------------------------------------------------------

    // Truncate to flank_size (plus some leeway for small indels in flanking region) to stop relatively distant
    // expansion sequences from accidentally being included in the flanking region; e.g. if the insert gets mapped
    // onto bases outside the definition coordinates.  TODO: better to do something else here?
    // The +10 here won't include any real TR region if the mapping is solid, since the flank coordinates will
    // contain a correctly-sized sequence.

    let flank_size = locus_with_ref_data.locus_def.flank_size;
    if let Some(lfe) = left_flank_end && lfe > flank_size + 10 && (left_flank_start < lfe - flank_size - 10) {
        left_flank_start = lfe - flank_size - 10;
    }
    if let Some(rfs) = right_flank_start && let Some(rfe) = right_flank_end && (rfe > rfs + flank_size + 10) {
        right_flank_end = Some(rfs + flank_size + 10);
    }

    // -----------------------------------------------------------------------------------------------------------------

    LocusReadCoords {
        left_flank_start: Some(left_flank_start),
        left_flank_end,
        right_flank_start,
        right_flank_end,
        full_left_flank,
        full_right_flank,
    }
}

#[pyfunction]
pub fn get_read_coords_from_matched_pairs(
    py: Python<'_>,
    locus_with_ref_data: &STRkitLocusWithRefData,
    segment: &mut STRkitAlignedSegment,
    aligned_coords: Option<Py<STRkitAlignedCoords>>,
    vcf_anchor_size: u64,
    allow_only_one_full_flank: bool,
) -> PyResult<Py<LocusReadCoords>> {
    if aligned_coords.is_none() {
        segment.cache_cigar_aligned_coords(py)?;
    }

    let ac = aligned_coords.as_ref()
        .or(segment.cigar_aligned_coords.as_ref())
        .expect("should have aligned coordinates from function or segment");

    let ac_ref = ac.borrow(py);

    Py::new(py, _get_read_coords_from_matched_pairs(
        locus_with_ref_data, &segment.query_sequence, ac_ref, vcf_anchor_size, allow_only_one_full_flank
    ))
}

#[pyfunction]
pub fn process_read_snvs_for_locus_and_calculate_useful_snvs(
    py: Python<'_>,
    block_segments: &STRkitLocusBlockSegments,
    locus_with_ref_data: &STRkitLocusWithRefData,
    left_most_coord: RefCoord,
    ref_cache: &str,
    read_dict_extra: Bound<PyDict>,
    read_locus_alignment_data: &Bound<PyDict>,
    candidate_snvs: &Bound<'_, CandidateSNVs>,
    min_allele_reads: usize,
    significant_clip_snv_take_in: usize,
    only_known_snvs: bool,
    logger: Bound<PyAny>,
    locus_log_str: &str,
) -> Result<Vec<(usize, RefCoord)>, PyErr> {
    // Loop through a second time if we are using SNVs. We do a second loop rather than just using the first loop
    // in order to have collected the edges of the reference sequence we can cache for faster SNV calculation.

    let left_coord_adj = locus_with_ref_data.left_coord_adj;
    let right_coord_adj = locus_with_ref_data.left_coord_adj;

    // Mutates: read_dict_extra (via call to calculate_useful_snvs)

    let mut locus_snvs: HashSet<RefCoord> = HashSet::new();
    let mut read_snvs: HashMap<String, HashMap<RefCoord, (char, u8)>> = HashMap::new();

    // below: magic values for skipping false positives / weird 'SNVs' that aren't helpful
    let useful_snvs_params = UsefulSNVsParams {
        contiguous_threshold: 5,
        max_snv_group_size: 5,
        too_many_snvs_threshold: 20,
        entropy_flank_size: 10,
        entropy_threshold: 1.8,
    };

    let twox_takein = significant_clip_snv_take_in * 2;

    let candidate_snvs_b = candidate_snvs.borrow();

    for rn in read_dict_extra.keys().into_iter().map(|x| x.cast_into::<PyString>().unwrap()) {
        let segment = block_segments.get_segment_by_name(rn.as_borrowed().to_str()?).unwrap().borrow(py);

        if segment.sig_clip_left || segment.sig_clip_right {
            logger.call_method1(
                intern!(py, "debug"),
                (
                    intern!(
                        py,
                        "%s - %s has significant clipping; trimming pairs by %d bp clipped per side for SNV-finding"
                    ),
                    locus_log_str,
                    &rn,
                    significant_clip_snv_take_in,
                ),
            )?;
        }

        let scl = if segment.sig_clip_left { significant_clip_snv_take_in } else { 0 };
        let scr = if segment.sig_clip_right { significant_clip_snv_take_in } else { 0 };

        let segment_alignment_data_for_locus =
            read_locus_alignment_data
                .get_item(&rn)?
                .unwrap()
                .cast_into::<STRkitSegmentAlignmentDataForLocus>()?
                .borrow();
        let aligned_coords = &segment_alignment_data_for_locus.aligned_coords.borrow(py);
        let coords_len = aligned_coords.query_coords.len();

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

        let snvs: HashMap<RefCoord, (char, u8)> = segment.get_read_snvs(
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

    calculate_useful_snvs(
        py, block_segments, read_dict_extra, read_locus_alignment_data, read_snvs, locus_snvs, min_allele_reads
    )
}
