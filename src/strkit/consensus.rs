use bio::alignment::pairwise::Scoring;
use bio::alignment::poa::*;
use pyo3::intern;
use pyo3::prelude::*;
use pyo3::pybacked::PyBackedStr;
use pyo3::types::PyString;
use std::cmp;
use std::collections::HashSet;
use std::ops::Deref;
use std::panic;
use strsim::normalized_levenshtein;

static GAP_CHAR_ORD: usize  = b'-' as usize;

static BLANK_STR: &str = "";


fn seqs_summed_levenshtein_similarity(seqs: &[&str]) -> Vec<f64> {
    let n_seqs = seqs.len();
    let mut similarity_memo: Vec<Option<f64>> = vec![None; n_seqs.pow(2)];
    let mut ds = vec![0f64; n_seqs];

    for i in 0..n_seqs {
        for j in 0..n_seqs {
            if i == j { continue; }
            // diagonal 2d matrix into 1d matrix (vector) index:
            let memo_idx = cmp::min(i, j) * n_seqs + cmp::max(i, j);
            if let Some(d) = similarity_memo[memo_idx] {
                // Distance already computed; use the cached value
                ds[i] += d;
            } else {
                // Need to compute the distance
                let d = normalized_levenshtein(seqs[i], seqs[j]);
                similarity_memo[memo_idx] = Some(d);
                ds[i] += d;
            }
        }
    }

    ds
}

fn best_representatives<'a>(seqs: &'a [&str]) -> HashSet<&'a str> {
    let ds = seqs_summed_levenshtein_similarity(seqs);
    let ms = ds.iter().max_by(|a, b| a.total_cmp(b));
    ms.map(|max_score| {
        ds.iter().enumerate().filter(|&(_, s)| s == max_score).map(|(i, _)| seqs[i]).collect()
    }).unwrap_or(HashSet::new())
}

fn best_representative<'a>(seqs: &'a [&str]) -> Option<&'a str> {
    /*
    Slightly different from a true consensus - returns the string with the minimum Levenshtein distance to all other
    strings for a particular allele. This roughly approximates a true consensus when |seqs| is large. If more than one
    best representative exist, the first one is returned. If |best| == |seqs| or |best| == 0, None is returned since there
    is effectively no true consensus.
    */
    best_representatives(seqs).into_iter().next()
}


fn majority_consensus_from_msa(aligned_seqs: &[&[u8]], aligned_len: usize) -> String {
    (0..aligned_len).filter_map(|i| {
        let mut counter = [0usize; 256];
        aligned_seqs.iter().for_each(|s| {
            counter[s[i] as usize] += 1;
        });
        let mode_idx = counter
            .iter()
            .enumerate()
            .max_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap_or(cmp::Ordering::Equal))
            .map(|(idx, _)| idx)
            .unwrap();

        (mode_idx != 255 && mode_idx != GAP_CHAR_ORD).then_some(mode_idx as u8 as char)
    }).collect::<String>()
}


fn poa_consensus_seq(seqs: &[&str]) -> Option<String> {
    let n_seqs = seqs.len();

    if n_seqs == 0 {
        return None;
    }

    panic::catch_unwind(|| {

        let scoring = Scoring::new(-1, 0, |a: u8, b: u8| if a == b { 1i32 } else { -1i32 });

        let _seqs: Vec<&[u8]> = seqs.iter().map(|s| s.as_bytes()).collect();

        let first_seq = _seqs[0];

        if n_seqs == 1 {
            return Some(first_seq.iter().map(|&b| b as char).collect::<String>());
        }

        let mut aligner = Aligner::new(scoring, first_seq);

        let mut max_len: usize = first_seq.len();
        _seqs[1..].iter().for_each(|y| {
            max_len = cmp::max(max_len, y.len());
            aligner.global(y).add_to_graph();
        });

        let pretty = aligner
            .alignment()
            .pretty(aligner.consensus().as_slice(), _seqs, aligner.graph(), max_len * 4);

        let pretty_split = pretty.split('\n').skip(1);

        let mut aligned_seqs: Vec<&[u8]> = Vec::with_capacity(n_seqs);
        let mut max_aligned_len: usize = 0;

        for s in pretty_split.into_iter().filter(|&z| !z.is_empty()) {
            let s_seq: &[u8] = s.split('\t').nth(1).unwrap().as_bytes();
            max_aligned_len = cmp::max(max_aligned_len, s_seq.len());
            aligned_seqs.push(s_seq);
        }

        Some(majority_consensus_from_msa(&aligned_seqs, max_aligned_len))
    }).ok()?
}


fn run_best_representatives<'py>(py: Python<'py>, seqs: &[&str], logger: Bound<'py, PyAny>) -> Option<(String, &'py Bound<'py, PyString>)> {
    best_representative(seqs).map_or_else(|| {
        logger
            .call_method1(
                intern!(py, "debug"),
                (intern!(py, "Got no best representative from sequences"),),
            )
            .unwrap();
        None
    }, |best_rep| Some((String::from(best_rep), intern!(py, "best_rep"))))
}


#[pyfunction]
pub fn consensus_seq<'py>(py: Python<'py>, seqs: Vec<PyBackedStr>, logger: Bound<'py, PyAny>, max_mdn_poa_length: usize) -> Option<(String, &'py Bound<'py, PyString>)> {
    let mut n_seqs = seqs.len();

    if seqs.is_empty() {
        return None;
    }

    let n_blanks = seqs.iter().filter(|&s| s == BLANK_STR).count();

    if n_blanks as f64 > (n_seqs as f64) / 2.0 {
        // blanks make up majority, so blank is the consensus
        return Some((String::from(""), if n_blanks == n_seqs { intern!(py, "single") } else { intern!(py, "best_rep") }));
    }

    // if blanks make up minority, filter them out for consensus
    n_seqs -= n_blanks;

    let mut seqs_no_blanks: Vec<PyBackedStr> = seqs.into_iter().filter(|x| !x.is_empty()).collect();

    let seqs_set = HashSet::<&PyBackedStr>::from_iter(seqs_no_blanks.iter());

    match seqs_set.len() {
        // With 1 sequence in the set, return it as the single value (if no blanks also present), or the majority value
        // otherwise.
        1 => Some((
            seqs_no_blanks[0].to_string(),
            if n_blanks == 0 { intern!(py, "single") } else { intern!(py, "best_rep") },
        )),
        // With 2 sequences, return the majority representative, using the first item in the sorted deduplicated
        // sequence vector as a tiebreaker.
        2 => {
            let mut seqs_set_vec: Vec<&PyBackedStr> = seqs_set.into_iter().collect();
            seqs_set_vec.sort_unstable();
            let i0_count = seqs_no_blanks.iter().filter(|&s| s == seqs_set_vec[0]).count();
            // Shortcut: if index 0 count < n_seqs / 2, usize cast is 1 (so we return index 1); otherwise return index 0.
            Some((seqs_set_vec[(i0_count < n_seqs / 2) as usize].to_string(), intern!(py, "best_rep")))
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

            poa_consensus_seq(&sv).map_or_else(|| {
                logger.call_method1(
                    intern!(py, "error"),
                    (intern!(py, "Got no POA consensus sequence from sequences %s; trying best representative strategy"), sv.clone()),
                ).unwrap();
                run_best_representatives(py, &sv, logger)
            }, |poa_c| Some((poa_c, intern!(py, "poa"))))
        }
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_seqs_summed_levenshtein_similarity() {
        assert_eq!(seqs_summed_levenshtein_similarity(&vec!["A", "A", "A"]), vec![2.0, 2.0, 2.0]);
        assert_eq!(seqs_summed_levenshtein_similarity(&vec!["AA", "AA", "AB"]), vec![1.5, 1.5, 1.0]);
    }

    #[test]
    fn test_best_representatives() {
        assert_eq!(best_representatives(&vec!["A", "A", "A"]), HashSet::from_iter(vec!["A"]));
        assert_eq!(best_representatives(&vec!["B", "A", "A"]), HashSet::from_iter(vec!["A"]));
        assert_eq!(best_representatives(&vec!["B", "A", "B"]), HashSet::from_iter(vec!["B"]));
        assert_eq!(best_representatives(&vec!["", "A", "A"]), HashSet::from_iter(vec!["A"]));
        assert_eq!(best_representatives(&vec!["A", "A", "B", "B"]), HashSet::from_iter(vec!["A", "B"]));
    }

    #[test]
    fn test_best_representative() {
        assert_eq!(best_representative(&vec![]), None);
        assert_eq!(best_representative(&vec!["A", "A", "A"]), Some("A"));
        assert_eq!(best_representative(&vec!["B", "A", "A"]), Some("A"));
        assert_eq!(best_representative(&vec!["B", "A", "B"]), Some("B"));
        assert_eq!(best_representative(&vec!["AA", "AB", "AA"]), Some("AA"));
        assert_eq!(best_representative(&vec!["AAACAAA", "AACAAA", "AAACAA"]), Some("AAACAAA"));
    }

    #[test]
    fn test_poa_consensus_seq() {
        assert_eq!(poa_consensus_seq(&vec![]), None);
        assert_eq!(poa_consensus_seq(&vec!["AA", "AB", "AA"]), Some(String::from("AA")));
    }
}
