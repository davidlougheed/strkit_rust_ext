use numpy::ndarray::{Array1, s};
use numpy::{PyArray, PyArray1, ToPyArray};
use pyo3::exceptions::PyValueError;
use pyo3::intern;
use pyo3::prelude::*;
use rust_htslib::bam::ext::BamRecordExtensions;
use rust_htslib::bam::record::{Aux, Record};
use rust_htslib::bam::{IndexedReader, Read};
use rust_htslib::errors::Error as RustHTSlibError;
use rust_lapper::{Interval, Lapper};
use std::cmp;
use std::collections::{HashMap, HashSet};
use std::sync::Mutex;

use crate::aligned_coords::{AlignedCoordsMethods, STRkitAlignedCoords};
use crate::cigar::{decode_cigar_item, get_aligned_pair_matches_rs};
use crate::exceptions::LowMeanBaseQual;
use crate::locus::{STRkitLocus, STRkitLocusBlock};
use crate::snvs::CandidateSNVs;
use crate::utils::{normalize_contig, starts_with_chr, calculate_seq_with_wildcards, calc_motif_size_kmers};

/// Locus-specific alignment data extracted from an aligned segment, possibly via realignment.
#[derive(Clone)]
#[pyclass(frozen)]
pub struct STRkitSegmentAlignmentDataForLocus {
    // Either cigar aligned coords from read, or realigned coords (which are locus-specific):
    pub aligned_coords: STRkitAlignedCoords,
    // ---
    pub left_flank_start: usize,
    #[pyo3(get)]
    pub left_flank_end: usize,
    pub right_flank_start: usize,
    pub right_flank_end: usize,
    // ---
    #[pyo3(get)]
    pub realigned: bool,
}

// TODO: remove this struct with future Rust porting in favour of just accessing .aligned_coords
impl AlignedCoordsMethods for STRkitSegmentAlignmentDataForLocus {
    fn query_coord_at_idx(&self, idx: usize) -> u64 {
        self.aligned_coords.query_coord_at_idx(idx)
    }
    fn find_coord_idx_by_ref_pos(&self, target: usize, start_left: usize) -> (usize, bool) {
        self.aligned_coords.find_coord_idx_by_ref_pos(target, start_left)
    }
}

#[pymethods]
impl STRkitSegmentAlignmentDataForLocus {
    #[new]
    fn py_new(
        aligned_coords: STRkitAlignedCoords,
        left_flank_start: usize,
        left_flank_end: usize,
        right_flank_start: usize,
        right_flank_end: usize,
        realigned: bool,
    ) -> PyResult<Self> {
        Ok(STRkitSegmentAlignmentDataForLocus {
            aligned_coords,
            left_flank_start,
            left_flank_end,
            right_flank_start,
            right_flank_end,
            realigned,
        })
    }

    pub fn query_coord_at_idx(&self, idx: usize) -> u64 {
        AlignedCoordsMethods::query_coord_at_idx(self, idx)
    }

    pub fn find_coord_idx_by_ref_pos(&self, target: usize, start_left: usize) -> (usize, bool) {
        AlignedCoordsMethods::find_coord_idx_by_ref_pos(self, target, start_left)
    }
}

#[pyclass(frozen)]
pub struct STRkitAlignedSegmentSequenceDataForLocus {
    // This struct must be created with STRkitAlignedSegment.get_sequence_data_for_locus.
    // ==================================================================================
    _motif_size: usize,
    _flank_size: usize,
    pub tr_qqs: Array1<u8>, // Tandem repeat region query sequence qualities
    pub mean_base_qual: Option<f32>, // Mean base quality of tandem repeat region query sequence
    // --
    #[pyo3(get)]
    pub flank_left_seq_wc: String,
    #[pyo3(get)]
    pub flank_right_seq_wc: String,
    #[pyo3(get)]
    pub tr_seq: String,
    #[pyo3(get)]
    pub tr_seq_wc: String,
    pub tr_len: usize,
    #[pyo3(get)]
    pub tr_len_with_flank: usize, // len(flank_left_seq) + len(tr_seq) + len(flank_right_seq)
}

#[pymethods]
impl STRkitAlignedSegmentSequenceDataForLocus {
    /// For updating read k-mers counter:
    ///  TODO: refact as iterator
    fn get_motif_size_kmers(&self) -> Vec<&str> {
        calc_motif_size_kmers(&self.tr_seq_wc, self.tr_len, self._motif_size)
    }

    /// Pretty basic - divides TR-aligned sequence length by the motif size to get the estimated copy number.
    /// Used as a starting point for optimizing.
    fn get_est_copy_num(&self) -> usize {
        (self.tr_len as f32 / self._motif_size as f32).round() as usize
    }

    /// Adjusts a pairwise alignment score (integer) from parasail into a floating point score (roughly) normalized to
    /// the length of the TR sequence + left/right flanking sequence.
    fn calc_adj_score(&self, read_cn_score: i32) -> Option<f32> {
        (self.tr_len > 0).then(|| read_cn_score as f32 / (self.tr_len + self._flank_size * 2) as f32)
    }
}

/// A wrapper struct for a read aligned to a reference genome, with some data pre-extracted from the aligned read file
/// and some STRkit-specific methods for the genotyping pipeline.
#[derive(Clone)]
#[pyclass]
pub struct STRkitAlignedSegment {
    #[pyo3(get)]
    name: String,
    #[pyo3(get)]
    length: usize,
    #[pyo3(get)]
    pub start: u64,
    #[pyo3(get)]
    pub end: u64,
    #[pyo3(get)]
    is_reverse: bool,
    #[pyo3(get)]
    is_supplementary: bool,
    #[pyo3(get)]
    is_secondary: bool,
    #[pyo3(get)]
    pub query_sequence: String,
    pub query_qualities_: Array1<u8>,
    pub raw_cigar: Array1<u32>,
    pub cigar_aligned_coords: Option<STRkitAlignedCoords>,
    cigar_first_op: u32, // first grabbed for significant clipping check, later used for soft clip overlap calc.
    cigar_first_len: u64,
    cigar_last_op: u32,  // idem
    cigar_last_len: u64,
    // ---
    // Observed significant increase in annoying, probably false SNVs near the edges of significantly clipped reads in
    // CCS data. Store if we have large clipping for later use here in the SNV finder.
    pub sig_clip_left: bool,
    pub sig_clip_right: bool,
    // ---
    #[pyo3(get)]
    hp: Option<i64>,
    #[pyo3(get)]
    ps: Option<i64>,
}

fn _extract_i64_tag_value(a: Result<Aux<'_>, RustHTSlibError>) -> Option<i64> {
    match a {
        Ok(value) => match value {
            Aux::U8(v) => Some(v as i64),
            Aux::U16(v) => Some(v as i64),
            Aux::U32(v) => Some(v as i64),
            Aux::I8(v) => Some(v as i64),
            Aux::I16(v) => Some(v as i64),
            Aux::I32(v) => Some(v as i64),
            _ => None,
        },
        Err(_) => None,
    }
}

impl STRkitAlignedSegment {
    pub fn query_qualities_slice(&self) -> &[u8] {
        self.query_qualities_.as_slice().unwrap()
    }

    pub fn cache_cigar_aligned_coords(&mut self) {
        match self.cigar_aligned_coords {
            Some(_) => (), // Already cached
            None => {
                let nac = get_aligned_pair_matches_rs(self.raw_cigar.view(), 0, self.start);
                self.cigar_aligned_coords = Some(nac);
            },
        }
    }
}

#[pymethods]
impl STRkitAlignedSegment {
    #[getter]
    fn query_sequence_bytes<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyArray1<u8>>> {
        Ok(PyArray1::from_array(py, &Array1::from_iter(self.query_sequence.clone().as_bytes().iter().copied())))
    }

    #[getter]
    fn query_qualities<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyArray1<u8>>> {
        Ok(PyArray1::from_array(py, &self.query_qualities_))
    }

    /// Whether we have soft-clipping (CIGAR op: 4) which overlaps our TR region + flank. Used when deciding whether to
    /// try to realign a read with parasail.
    fn soft_clip_overlaps_locus(&self, locus: &STRkitLocus) -> bool {
        let lfc = locus.left_flank_coord as u64;
        let rfc = locus.right_flank_coord as u64;

        (self.cigar_first_op == 4 && self.start > lfc && lfc >= self.start - self.cigar_first_len)
        || (self.cigar_last_op == 4 && self.end < rfc && rfc <= self.end + self.cigar_last_len)
    }

    fn get_sequence_data_for_locus<'py>(
        &self,
        py: Python<'py>,
        locus: &STRkitLocus,
        segment_alignment_data_for_locus: &STRkitSegmentAlignmentDataForLocus,
        min_avg_phred: f32,
        base_wildcard_threshold: u8,
    ) -> PyResult<STRkitAlignedSegmentSequenceDataForLocus> {
        let left_flank_start = segment_alignment_data_for_locus.left_flank_start;
        let left_flank_end = segment_alignment_data_for_locus.left_flank_end;
        let right_flank_start = segment_alignment_data_for_locus.right_flank_start;
        let right_flank_end = segment_alignment_data_for_locus.right_flank_end;

        // -------------------------------------------------------------------------------------------------------------

        let tr_qqs = self.query_qualities_.slice(s![left_flank_end..right_flank_start]).to_owned();
        // Built-in mean was panicking; DIY it
        let mean_base_qual = if right_flank_start > left_flank_end + 1 {
            Some(tr_qqs.iter().map(|&e| e as f32).sum::<f32>() / ((right_flank_start - left_flank_end) as f32))
        } else { None };

        if let Some(mbq) = mean_base_qual {
            if mbq < min_avg_phred {
                // Low average base quality, we will skip this read. Debug logging will be handled on the Python side.
                let err = LowMeanBaseQual { mean_base_qual: mbq };
                return Err(PyErr::from_value(err.into_pyobject(py)?));
            }
        }

        // -------------------------------------------------------------------------------------------------------------

        // Truncate to flank_size (plus some leeway for small indels in flanking region) to stop relatively distant
        // expansion sequences from accidentally being included in the flanking region; e.g. if the insert gets mapped
        // onto bases outside the definition coordinates.  TODO: better to do something else here?
        // The +10 here won't include any real TR region if the mapping is solid, since the flank coordinates will
        // contain a correctly-sized sequence.

        let flank_left_seq = &self.query_sequence[left_flank_start..left_flank_end];
        let flank_right_seq = &self.query_sequence[right_flank_start..right_flank_end];

        let tr_len = right_flank_start - left_flank_end;
        let tr_read_seq= String::from(&self.query_sequence[left_flank_end..right_flank_start]); // TODO: rust: zero-copy

        // Cache this now that we've done some adjustments (although tr_len_w_flank should not change)
        let tr_len_with_flank = tr_len + flank_left_seq.len() + flank_right_seq.len();

        // -------------------------------------------------------------------------------------------------------------

        // Extract qualities for our TR sequence and flanking regions, and calculate wildcard-subsituted sequences if
        // qualities are too low to be confident in base identity.

        let flank_left_qqs = &self.query_qualities_.slice(s![left_flank_start..left_flank_end]);
        let flank_right_qqs = &self.query_qualities_.slice(s![right_flank_start..right_flank_end]);

        let flank_left_seq_wc = calculate_seq_with_wildcards(
            flank_left_seq, Some(flank_left_qqs), base_wildcard_threshold
        );
        let flank_right_seq_wc = calculate_seq_with_wildcards(
            flank_right_seq, Some(flank_right_qqs), base_wildcard_threshold
        );
        let tr_read_seq_wc = calculate_seq_with_wildcards(
            &tr_read_seq, Some(&tr_qqs.view()), base_wildcard_threshold
        );

        // !!! right_flank_start and right_flank_end are now invalid !!!

        // -------------------------------------------------------------------------------------------------------------

        Ok(
            // TODO: rust: 0-copy + same lifetime
            STRkitAlignedSegmentSequenceDataForLocus {
                _motif_size: locus.motif_size,
                _flank_size: locus.flank_size as usize,
                tr_qqs: tr_qqs, // TODO: clone or view?
                mean_base_qual,
                // --
                flank_left_seq_wc,
                flank_right_seq_wc,
                tr_seq: tr_read_seq,
                tr_seq_wc: tr_read_seq_wc,
                tr_len,
                tr_len_with_flank,
            }
        )
    }
}


/// In STRkit, we extract read data for loci in blocks to minimize read file accesses (BAM/CRAM). This is a container
/// struct for these aligned reads/segments for a block of loci, with an interval tree for quickly extracting reads
/// overlapping a specific locus.
/// This struct must be initialized by STRkitBAMReader.get_overlapping_segments_and_related_data_for_block(...)
#[pyclass]
pub struct STRkitLocusBlockSegments {
    #[pyo3(get)]
    left_most_coord: u64,
    #[pyo3(get)]
    right_most_coord: u64,
    pub segments: Vec<STRkitAlignedSegment>,
    // hash map of read name to index for fast lookup of segment sequence/query qualities:
    pub name_index_lookup: HashMap<String, usize>,
    tree: Mutex<Lapper<usize, usize>>,
    logger: Option<Py<PyAny>>,
    // params:
    skip_supp: bool,
    skip_sec: bool,
    max_locus_reads: usize,
}

impl STRkitLocusBlockSegments {
    pub fn get_segment_by_name(&self, rn: &str) -> Option<&STRkitAlignedSegment> {
        self.name_index_lookup.get(rn).map(|&si| &self.segments[si])
    }
}

#[pymethods]
impl STRkitLocusBlockSegments {
    fn get_segments_for_locus<'py>(
        &mut self,
        py: Python<'py>,
        locus: Bound<'_, STRkitLocus>,
    ) -> PyResult<(
        Bound<'py, PyArray1<Py<PyAny>>>,
        usize,
        Bound<'py, PyArray1<usize>>,
        HashMap<String, u8>,
        u64,
        u64,
    )> {
        let tree = self.tree.lock().unwrap();
        let loc = locus.borrow();

        let mut left_most_coord = 999999999999u64;
        let mut right_most_coord = 0u64;

        let mut segments: Vec<Py<STRkitAlignedSegment>> = Vec::new();
        let mut read_lengths: Vec<usize> = Vec::new();
        let mut chimeric_read_status: HashMap<String, u8> = HashMap::new();
        let mut seen_reads: HashSet<String> = HashSet::new();

        // Fetch every segment (interval from the tree) which OVERLAPS [locus left flank coord, locus right flank coord]
        for i in tree.find(loc.left_flank_coord as usize, loc.right_flank_coord as usize) {
            let seg = &self.segments[i.val];

            // If we have two overlapping alignments for the same read, we have a chimeric read within the TR locus
            // (so probably a large expansion...)
            let crs: &mut u8 = chimeric_read_status.entry(seg.name.clone()).or_insert(0u8);
            *crs |= if seg.is_supplementary { 2u8 } else { 1u8 };

            if self.skip_supp && seg.is_supplementary {  // If configured, skip supplementary alignments
                if let Some(logger) = &self.logger {
                    // Keep debug log level check in Rust to avoid needless Python call
                    logger.call_method1(
                        py,
                        intern!(py, "debug"),
                        (
                            intern!(py, "%s - skipping entry for read %s (supplementary)"),
                            locus.borrow().log_str(),
                            &seg.name,
                        ),
                    )?;
                }
                continue;
            }

            if self.skip_sec && seg.is_secondary {  // If configured, skip secondary alignments
                if let Some(logger) = &self.logger {
                    // Keep debug log level check in Rust to avoid needless Python call
                    logger.call_method1(
                        py,
                        intern!(py, "debug"),
                        (
                            intern!(py, "%s - skipping entry for read %s (secondary)"),
                            locus.borrow().log_str(),
                            &seg.name,
                        ),
                    )?;
                }
                continue;
            }

            // Copy them into a new vector (--> array) with associated data, replacing the read-level fetches we used to do.
            segments.push(Py::new(py, seg.clone())?);
            read_lengths.push(seg.length);
            seen_reads.insert(seg.name.clone());

            left_most_coord = cmp::min(left_most_coord, seg.start);
            right_most_coord = cmp::max(right_most_coord, seg.end);

            if seen_reads.len() > self.max_locus_reads {
                // We specifically break when we're over the maximum so that Python can also see we're over
                // the threshold and log something.
                break;
            }
        }

        let n_segments = segments.len();

        Ok((
            PyArray::from_owned_object_array(py, Array1::from_vec(segments)),
            n_segments,
            read_lengths.to_pyarray(py),
            chimeric_read_status,
            left_most_coord,
            right_most_coord,
        ))
    }
}


#[pyclass]
pub struct STRkitBAMReader {
    reader: Mutex<IndexedReader>,
    max_locus_reads: usize,
    skip_supp: bool,
    skip_sec: bool,
    use_hp: bool,
    significant_clip_threshold: u64,
    contig_names: Vec<String>,
    any_contig_name_has_chr: bool,
    logger: Py<PyAny>,
    debug_logs: bool,
}

#[pymethods]
impl STRkitBAMReader {
    #[new]
    fn py_new(
        py: Python<'_>,
        path: &str,
        ref_path: &str,
        max_locus_reads: usize,
        skip_supp: bool,
        skip_sec: bool,
        use_hp: bool,
        significant_clip_threshold: u64,
        logger: Bound<PyAny>,
        debug_logs: bool,
    ) -> PyResult<Self> {
        let r = IndexedReader::from_path(path);

        if let Ok(mut rdr) = r {
            rdr.set_reference(ref_path).unwrap();

            let contig_names: Vec<String> = {
                let names = rdr.header().target_names();
                names.into_iter().map(|n| String::from_utf8_lossy(n).to_string()).collect()
            };

            let any_contig_name_has_chr = contig_names.iter().any(|c| starts_with_chr(c));

            let reader = Mutex::new(rdr);

            Ok(
                STRkitBAMReader {
                    reader,
                    max_locus_reads,
                    skip_supp,
                    skip_sec,
                    use_hp,
                    significant_clip_threshold,
                    contig_names,
                    any_contig_name_has_chr,
                    logger: logger.unbind().clone_ref(py),
                    debug_logs,
                }
            )
        } else {
            Err(PyErr::new::<PyValueError, _>(format!("Could not load BAM from path: {}", path)))
        }
    }

    fn get_overlapping_segments_and_related_data_for_block<'py>(
        &mut self,
        py: Python<'py>,
        locus_block: &STRkitLocusBlock,
    ) -> PyResult<Py<STRkitLocusBlockSegments>> {
        let mut reader = self.reader.lock().unwrap();

        let contig_norm = normalize_contig(&locus_block.loci[0].contig, self.any_contig_name_has_chr); // TODO: hash map
        reader.fetch((&contig_norm, locus_block.left, locus_block.right)).unwrap();

        let mut left_most_coord = 999999999999u64;
        let mut right_most_coord = 0u64;

        let mut segments: Vec<STRkitAlignedSegment> = Vec::new();
        let mut name_index_lookup: HashMap<String, usize> = HashMap::new(); // doubles as seen reads hashset
        let mut intervals: Vec<Interval<usize, usize>> = Vec::new();

        let mut record = Record::new();

        while let Some(r) = reader.read(&mut record) {
            match r {
                Ok(_) => {
                    let name = String::from_utf8_lossy(record.qname()).to_string();
                    let supp = record.is_supplementary();
                    let sec = record.is_secondary();

                    if name_index_lookup.contains_key(&name) {
                        if self.debug_logs {  // Keep debug log level check in Rust to avoid needless Python call
                            self.logger.call_method1(
                                py,
                                intern!(py, "debug"),
                                (
                                    intern!(py, "%s - skipping entry for read %s (already seen)"),
                                    &locus_block.log_str,
                                    name,
                                ),
                            )?;
                        }
                        continue;
                    }

                    let length = record.seq_len_from_cigar(false);

                    if length == 0 || record.seq_len() == 0 {
                        // No aligned segment, skip entry (used to pull query sequence, but that's extra work)
                        continue;
                    }

                    let start = record.reference_start() as u64;
                    let end = record.reference_end() as u64;
                    let is_reverse = record.is_reverse();
                    let query_sequence = String::from_utf8(record.seq().as_bytes())
                        .map_err(|e| PyErr::new::<PyValueError, _>(e.to_string()))?;

                    let raw_cigar_vec = record.raw_cigar().to_vec();

                    // -------------------------------------------------------------------------------------------------

                    // Observed significant increase in annoying, probably false SNVs near the edges of significantly
                    // clipped reads in CCS data. Figure out if we have large clipping for later use here in the SNV
                    // finder, by extracting the first/last CIGAR operations and checking for clipping.

                    let (cigar_first_op, cigar_first_len) = decode_cigar_item(raw_cigar_vec[0]);
                    let sig_clip_left =
                        (cigar_first_op == 4 || cigar_first_op == 5)
                        && cigar_first_len >= self.significant_clip_threshold;

                    let (cigar_last_op, cigar_last_len) = decode_cigar_item(raw_cigar_vec[raw_cigar_vec.len() - 1]);
                    let sig_clip_right =
                        (cigar_last_op == 4 || cigar_last_op == 5)
                        && cigar_last_len >= self.significant_clip_threshold;

                    // -------------------------------------------------------------------------------------------------

                    let aligned_segment = STRkitAlignedSegment {
                        name,
                        length,
                        start,
                        end,
                        is_reverse,
                        is_supplementary: supp,
                        is_secondary: sec,
                        query_sequence,
                        query_qualities_: Array1::from_vec(record.qual().to_vec()),
                        // --- begin cigar stuff
                        raw_cigar: Array1::from_vec(raw_cigar_vec),
                        cigar_aligned_coords: None,
                        cigar_first_op,
                        cigar_first_len,
                        cigar_last_op,
                        cigar_last_len,
                        sig_clip_left,
                        sig_clip_right,
                        // --- end cigar stuff
                        hp: if self.use_hp { _extract_i64_tag_value(record.aux(b"HP")) } else { None },
                        ps: if self.use_hp { _extract_i64_tag_value(record.aux(b"PS")) } else { None },
                    };

                    let nc = aligned_segment.name.clone();
                    name_index_lookup.insert(nc, segments.len());
                    segments.push(aligned_segment);
                    intervals.push(Interval { start: start as usize, stop: end as usize, val: segments.len() - 1 });

                    left_most_coord = cmp::min(left_most_coord, start);
                    right_most_coord = cmp::max(right_most_coord, end);
                },
                Err(_) => panic!("Error reading alignment record"),
            }
        }

        Py::new(py, STRkitLocusBlockSegments {
            left_most_coord,
            right_most_coord,
            segments,
            name_index_lookup,
            tree: Mutex::new(Lapper::new(intervals)),
            logger: self.debug_logs.then(|| self.logger.clone_ref(py)),
            // params:
            skip_supp: self.skip_supp,
            skip_sec: self.skip_sec,
            max_locus_reads: self.max_locus_reads,
        })
    }
}
