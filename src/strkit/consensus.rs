mod best_representatives;
mod poa;

use pyo3::intern;
use pyo3::prelude::*;
use pyo3::pybacked::PyBackedStr;
use pyo3::types::PyString;
use std::collections::HashSet;
// use std::io::BufWriter;
use std::ops::Deref;

use best_representatives::best_representative;
use poa::poa_consensus_seq;

static BLANK_STR: &str = "";


fn run_best_representatives<'py>(py: Python<'py>, seqs: &[&str], logger: Bound<'py, PyAny>) -> PyResult<Option<(String, &'py Bound<'py, PyString>)>> {
    match best_representative(seqs) {
        Some(best_rep) => Ok(Some((String::from(best_rep), intern!(py, "best_rep")))),
        None => {
            logger.call_method1(
                intern!(py, "debug"),
                (intern!(py, "Got no best representative from sequences"),),
            )?;
            Ok(None)
        }
    }
}


#[pyfunction]
pub fn consensus_seq<'py>(
    py: Python<'py>, seqs: Vec<PyBackedStr>, logger: Bound<'py, PyAny>, max_mdn_poa_length: usize
) -> PyResult<Option<(String, &'py Bound<'py, PyString>)>> {
    let mut n_seqs = seqs.len();

    if seqs.is_empty() {
        return Ok(None);
    }

    let n_blanks = seqs.iter().filter(|&s| s == BLANK_STR).count();

    if n_blanks as f64 > (n_seqs as f64) / 2.0 {
        // blanks make up majority, so blank is the consensus
        return Ok(
            Some((String::from(""), if n_blanks == n_seqs { intern!(py, "single") } else { intern!(py, "best_rep") }))
        );
    }

    // if blanks make up minority, filter them out for consensus
    n_seqs -= n_blanks;

    let mut seqs_no_blanks: Vec<PyBackedStr> = seqs.into_iter().filter(|x| !x.is_empty()).collect();

    let seqs_set = HashSet::<&PyBackedStr>::from_iter(seqs_no_blanks.iter());

    match seqs_set.len() {
        // With 1 sequence in the set, return it as the single value (if no blanks also present), or the majority value
        // otherwise.
        1 => Ok(Some((
            seqs_no_blanks[0].to_string(),
            if n_blanks == 0 { intern!(py, "single") } else { intern!(py, "best_rep") },
        ))),
        // With 2 sequences, return the majority representative, using the first item in the sorted deduplicated
        // sequence vector as a tiebreaker.
        2 => {
            let mut seqs_set_vec: Vec<&PyBackedStr> = seqs_set.into_iter().collect();
            seqs_set_vec.sort_unstable();
            let i0_count = seqs_no_blanks.iter().filter(|&s| s == seqs_set_vec[0]).count();
            // Shortcut: if index 0 count < n_seqs / 2, usize cast is 1 (so we return index 1); otherwise return index 0.
            Ok(Some((seqs_set_vec[(i0_count < n_seqs / 2) as usize].to_string(), intern!(py, "best_rep"))))
        },
        // Otherwise, return the POA alignment consensus or, if the sequences are too long to quickly run POA, or too
        // short to not crash the current version of rust-bio's POA function, the best representative strategy is used
        // instead.
        _ => {
            // sort the sequences by size so we can get the median size
            seqs_no_blanks.sort_unstable_by_key(|a| a.len());

            let sv: Vec<&str> = seqs_no_blanks.iter().map(|s| s.deref()).collect();

            // if the sequences are too long to quickly run POA --> we cannot quickly run POA
            // if the sequences are too short --> rust-bio's POA may give the wrong result or crash
            //   -> tracking: https://github.com/rust-bio/rust-bio/pull/605
            let mdn_seq_len = seqs_no_blanks[n_seqs / 2].len();
            if mdn_seq_len <= 1 || mdn_seq_len > max_mdn_poa_length {
                return run_best_representatives(py, &sv, logger);
            }

            match poa_consensus_seq(&sv) {
                Some(poa_c) => Ok(Some((poa_c, intern!(py, "poa")))),
                None => {
                    logger.call_method1(
                        intern!(py, "error"),
                        (intern!(py, "Got no POA consensus sequence from sequences %s; trying best representative strategy"), sv.clone()),
                    )?;
                    run_best_representatives(py, &sv, logger)
                }
            }
        }
    }
}
