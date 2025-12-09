use numpy::{ndarray::Array1, PyArray1, PyArrayMethods};
use pyo3::{prelude::*, pybacked::PyBackedStr};
use std::iter::zip;

// TODO: use this
pub fn calculate_seq_with_wildcards(qs: &str, quals: Option<Array1<u8>>, base_wildcard_threshold: u8) -> String {
    if let Some(qls) = quals {
        zip(qs.chars(), qls).map(|(c, q)| {
            if q > base_wildcard_threshold { c } else { 'X' }
        }).collect::<String>()
    } else {
        qs.to_owned()
    }
}

#[pyfunction]
#[pyo3(name = "calculate_seq_with_wildcards")]
pub fn calculate_seq_with_wildcards_py(
    qs: PyBackedStr,
    quals: Option<&Bound<'_, PyArray1<u8>>>,
    base_wildcard_threshold: u8,
) -> String {
    if let Some(qls) = quals {
        let qls_ro = qls.readonly();
        zip(qs.chars(), qls_ro.as_array()).map(|(c, &q)| {
            if q > base_wildcard_threshold { c } else { 'X' }
        }).collect::<String>()
    } else {
        qs.to_string()
    }
}

pub fn starts_with_chr(contig: &str) -> bool {
    // Whether or not a contig name starts with the string "chr"
    contig.starts_with("chr")
}

pub fn normalize_contig(contig: &str, has_chr: bool) -> String {
    // Given a contig str and a boolean of whether we should have the UCSC-style "chr" prefix, normalize the contig
    // str to be either UCSC-style (e.g., chr1/chr2/.../chrX/chrY/chrM) or not (e.g., 1/2/.../X/Y/M).
    let mut res = String::new();
    if has_chr { res.push_str("chr"); }
    res.push_str(if starts_with_chr(contig) { &contig[3..] } else { contig });
    res
}

#[pyfunction]
#[pyo3(name = "normalize_contig")]
pub fn normalize_contig_py(contig: PyBackedStr, has_chr: bool) -> String {
    // Same as above, but for Python.
    // TODO: once every use moves to Rust, remove this.
    normalize_contig(&contig, has_chr)
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_starts_with_chr() {
        assert!(starts_with_chr("chr1"));
        assert!(!starts_with_chr("ch1"));
        assert!(!starts_with_chr("1chr"));
        assert!(!starts_with_chr("1"));
        assert!(!starts_with_chr("M"));
    }

    #[test]
    fn test_normalize_contig() {
        assert_eq!(normalize_contig("chr1", true), "chr1");
        assert_eq!(normalize_contig("chrM", true), "M");
        assert_eq!(normalize_contig("chr1", false), "1");
        assert_eq!(normalize_contig("1", true), "chr1");
        assert_eq!(normalize_contig("1", false), "1");
        assert_eq!(normalize_contig("X", true), "chrX");
        assert_eq!(normalize_contig("Y", false), "Y");
    }
}
