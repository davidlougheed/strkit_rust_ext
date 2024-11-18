use std::cmp;
use std::collections::{HashMap, HashSet};
use numpy::ndarray::Array1;
use numpy::{PyArray1, ToPyArray};
use pyo3::exceptions::PyValueError;
use pyo3::intern;
use pyo3::prelude::*;
use rust_htslib::bam::ext::BamRecordExtensions;
use rust_htslib::bam::{IndexedReader, Read};
use rust_htslib::bam::record::{Aux, Record};
use rust_htslib::errors::Error as RustHTSlibError;

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
        Ok(value) => {
            match value {
                Aux::U8(v) => Some(v as i64),
                Aux::U16(v) => Some(v as i64),
                Aux::U32(v) => Some(v as i64),
                Aux::I8(v) => Some(v as i64),
                Aux::I16(v) => Some(v as i64),
                Aux::I32(v) => Some(v as i64),
                _ => None
            }
        },
        Err(_) => None
    }
}

#[pymethods]
impl STRkitAlignedSegment {
    #[getter]
    fn query_qualities<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyArray1<u8>>> {
        Ok(PyArray1::from_array_bound(py, &self._query_qualities))
    }

    #[getter]
    fn raw_cigar<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyArray1<u32>>> {
        Ok(PyArray1::from_array_bound(py, &self._raw_cigar))
    }
}


#[pyclass]
pub struct STRkitBAMReader {
    reader: IndexedReader,
}

#[pymethods]
impl STRkitBAMReader {
    #[new]
    fn py_new(path: &str, ref_path: &str) -> PyResult<Self> {
        let r = IndexedReader::from_path(path);

        if let Ok(mut reader) = r {
            reader.set_reference(ref_path).unwrap();
            Ok(STRkitBAMReader { reader })
        } else {
            Err(PyErr::new::<PyValueError, _>(format!("Could not load BAM from path: {}", path)))
        }
    }

    #[getter]
    fn references(&self) -> Vec<String> {
        let names = self.reader.header().target_names();
        names.into_iter().map(|n| String::from_utf8_lossy(n).to_string()).collect()
    }

    fn get_overlapping_segments_and_related_data<'py>(
        &mut self,
        py: Python<'py>,
        contig: &str,
        left_coord: i64,
        right_coord: i64,
        max_reads: usize,
        logger: Bound<PyAny>,
        locus_log_str: &str,
    ) -> PyResult<(Bound<'py, PyArray1<PyObject>>, usize, Bound<'py, PyArray1<usize>>, HashMap<String, u8>, i64, i64)> {
        self.reader.fetch((contig, left_coord, right_coord)).unwrap();

        let mut left_most_coord = 999999999999i64;
        let mut right_most_coord = 0i64;

        let mut segments: Vec<PyObject> = Vec::new();
        let mut read_lengths: Vec<usize> = Vec::new();
        let mut chimeric_read_status: HashMap<String, u8> = HashMap::new();
        let mut seen_reads: HashSet<String> = HashSet::new();

        let mut record = Record::new();

        while let Some(r) = self.reader.read(&mut record) {
            match r {
                Ok(_) => {
                    let name = String::from_utf8_lossy(record.qname()).to_string();

                    let supp = record.is_supplementary();

                    // If we have two overlapping alignments for the same read, we have a chimeric read within the TR
                    // (so probably a large expansion...)

                    let crs = chimeric_read_status.entry(name.clone()).or_insert(0u8);
                    *crs |= if supp {2u8} else {1u8};

                    if supp {  // Skip supplemental alignments
                        logger.call_method1(
                            intern!(py, "debug"),
                            (intern!(py, "%s - skipping entry for read %s (supplemental)"), locus_log_str, name),
                        )?;
                        continue;
                    }

                    if seen_reads.contains(&name) {  // Skip already-seen reads
                        logger.call_method1(
                            intern!(py, "debug"),
                            (intern!(py, "%s - skipping entry for read %s (already seen)"), locus_log_str, name),
                        )?;
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
                        query_sequence,
                        _query_qualities: Array1::from_vec(record.qual().to_vec()),
                        _raw_cigar: Array1::from_vec(record.raw_cigar().to_vec()),
                        hp: _extract_i64_tag_value(record.aux(b"HP")),
                        ps: _extract_i64_tag_value(record.aux(b"PS")),
                    };

                    let nc = aligned_segment.name.clone();
                    seen_reads.insert(nc);

                    segments.push(Py::new(py, aligned_segment)?.into_py(py));
                    read_lengths.push(length);

                    left_most_coord = cmp::min(left_most_coord, start);
                    right_most_coord = cmp::max(right_most_coord, end);

                    if seen_reads.len() > max_reads {
                        break;
                    }
                }
                Err(_) => panic!("Error reading alignment record")
            }
        }

        Ok((
            segments.to_pyarray_bound(py),
            segments.len(),
            read_lengths.to_pyarray_bound(py),
            chimeric_read_status,
            left_most_coord,
            right_most_coord,
        ))
    }
}
