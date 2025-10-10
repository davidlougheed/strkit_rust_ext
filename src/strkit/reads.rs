use numpy::ndarray::Array1;
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

use crate::locus::STRkitLocus;

#[derive(Clone)]
#[pyclass]
pub struct STRkitAlignedSegment {
    #[pyo3(get)]
    name: String,
    #[pyo3(get)]
    length: usize,
    #[pyo3(get)]
    start: i64,
    #[pyo3(get)]
    end: i64,
    #[pyo3(get)]
    is_reverse: bool,
    #[pyo3(get)]
    is_supplementary: bool,
    #[pyo3(get)]
    is_secondary: bool,
    #[pyo3(get)]
    query_sequence: String,
    _query_qualities: Array1<u8>,
    _raw_cigar: Array1<u32>,
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

#[pymethods]
impl STRkitAlignedSegment {
    #[getter]
    fn query_sequence_bytes<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyArray1<u8>>> {
        Ok(PyArray1::from_array(py, &Array1::from_iter(self.query_sequence.clone().as_bytes().iter().copied())))
    }

    #[getter]
    fn query_qualities<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyArray1<u8>>> {
        Ok(PyArray1::from_array(py, &self._query_qualities))
    }

    #[getter]
    fn raw_cigar<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyArray1<u32>>> {
        Ok(PyArray1::from_array(py, &self._raw_cigar))
    }
}


#[pyclass]
pub struct STRkitLocusBlockSegments {
    #[pyo3(get)]
    left_most_coord: i64,
    #[pyo3(get)]
    right_most_coord: i64,
    segments: Vec<STRkitAlignedSegment>,
    tree: Mutex<Lapper<usize, usize>>,
    logger: Option<Py<PyAny>>,
    // params:
    skip_supp: bool,
    skip_sec: bool,
    max_locus_reads: usize,
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
        i64,
        i64,
    )> {
        let tree = self.tree.lock().unwrap();
        let loc = locus.borrow();

        let mut left_most_coord = 999999999999i64;
        let mut right_most_coord = 0i64;

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
        logger: Bound<PyAny>,
        debug_logs: bool,
    ) -> PyResult<Self> {
        let r = IndexedReader::from_path(path);

        if let Ok(mut reader) = r {
            reader.set_reference(ref_path).unwrap();
            Ok(
                STRkitBAMReader {
                    reader: Mutex::new(reader),
                    max_locus_reads,
                    skip_supp,
                    skip_sec,
                    use_hp,
                    logger: logger.unbind().clone_ref(py),
                    debug_logs,
                }
            )
        } else {
            Err(PyErr::new::<PyValueError, _>(format!("Could not load BAM from path: {}", path)))
        }
    }

    #[getter]
    fn references(&self) -> Vec<String> {
        let reader = self.reader.lock().unwrap();
        let names = reader.header().target_names();
        names.into_iter().map(|n| String::from_utf8_lossy(n).to_string()).collect()
    }

    fn get_overlapping_segments_and_related_data_for_block<'py>(
        &mut self,
        py: Python<'py>,
        contig: &str,
        left_coord: i64,
        right_coord: i64,
        log_str: &str,
    ) -> PyResult<Py<STRkitLocusBlockSegments>> {
        let mut reader = self.reader.lock().unwrap();

        reader.fetch((contig, left_coord, right_coord)).unwrap();

        let mut left_most_coord = 999999999999i64;
        let mut right_most_coord = 0i64;

        let mut segments: Vec<STRkitAlignedSegment> = Vec::new();
        let mut intervals: Vec<Interval<usize, usize>> = Vec::new();
        let mut seen_reads: HashSet<String> = HashSet::new();

        let mut record = Record::new();

        while let Some(r) = reader.read(&mut record) {
            match r {
                Ok(_) => {
                    let name = String::from_utf8_lossy(record.qname()).to_string();
                    let supp = record.is_supplementary();
                    let sec = record.is_secondary();

                    if seen_reads.contains(&name) {  // Skip already-seen reads
                        if self.debug_logs {  // Keep debug log level check in Rust to avoid needless Python call
                            self.logger.call_method1(
                                py,
                                intern!(py, "debug"),
                                (intern!(py, "%s - skipping entry for read %s (already seen)"), log_str, name),
                            )?;
                        }
                        continue;
                    }

                    let length = record.seq_len_from_cigar(false);

                    if length == 0 || record.seq_len() == 0 {
                        // No aligned segment, skip entry (used to pull query sequence, but that's extra work)
                        continue;
                    }

                    let start = record.reference_start();
                    let end = record.reference_end();
                    let is_reverse = record.is_reverse();
                    let query_sequence = String::from_utf8(record.seq().as_bytes())
                        .map_err(|e| PyErr::new::<PyValueError, _>(e.to_string()))?;

                    let aligned_segment = STRkitAlignedSegment {
                        name,
                        length,
                        start,
                        end,
                        is_reverse,
                        is_supplementary: supp,
                        is_secondary: sec,
                        query_sequence,
                        _query_qualities: Array1::from_vec(record.qual().to_vec()),
                        _raw_cigar: Array1::from_vec(record.raw_cigar().to_vec()),
                        hp: if self.use_hp { _extract_i64_tag_value(record.aux(b"HP")) } else { None },
                        ps: if self.use_hp { _extract_i64_tag_value(record.aux(b"PS")) } else { None },
                    };

                    let nc = aligned_segment.name.clone();
                    seen_reads.insert(nc);

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
            tree: Mutex::new(Lapper::new(intervals)),
            logger: self.debug_logs.then(|| self.logger.clone_ref(py)),
            // params:
            skip_supp: self.skip_supp,
            skip_sec: self.skip_sec,
            max_locus_reads: self.max_locus_reads,
        })
    }
}
