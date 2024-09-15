use bio::alignment::pairwise::Scoring;
use bio::alignment::poa::*;
use pyo3::intern;
use pyo3::prelude::*;
use pyo3::types::PyString;
use std::cmp;
use std::collections::HashSet;
use std::panic;
use strsim::normalized_levenshtein;

static GAP_CHAR_ORD: usize  = b'-' as usize;

static BLANK_STR: &str = "";


fn best_representatives<'a>(seqs: &'a [&str]) -> HashSet<&'a str> {
    let mut ds = vec![0f64; seqs.len()];

    for i in 0..seqs.len() {
        for j in 0..seqs.len() {
            if i == j { continue; }
            ds[i] += normalized_levenshtein(seqs[i], seqs[j]);
        }
    }

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

    let mut res: Vec<&str> = Vec::from_iter(best_representatives(seqs));
    res.pop()
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
            let s_seq = s.split('\t').nth(1).unwrap().as_bytes();
            max_aligned_len = cmp::max(max_aligned_len, s_seq.len());
            aligned_seqs.push(s_seq);
        }

        Some(
            (0..max_aligned_len).filter_map(|i| {
                let mut counter = [0usize; 256];
                aligned_seqs.iter().for_each(|&s| {
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
        )
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
pub fn consensus_seq<'py>(py: Python<'py>, seqs: Vec<&str>, logger: Bound<'py, PyAny>, max_mdn_poa_length: usize) -> Option<(String, &'py Bound<'py, PyString>)> {
    let mut n_seqs = seqs.len();

    if seqs.is_empty() {
        return None;
    } else if HashSet::<&str>::from_iter(seqs.clone()).len() == 1 {
        return Some((seqs[0].to_owned(), intern!(py, "single")));
    }

    let n_blanks = seqs.iter().filter(|&&s| s == BLANK_STR).count();

    if n_blanks as f64 > (n_seqs as f64) / 2.0 {
        // blanks make up majority, so blank is the consensus
        return Some((String::from(""), intern!(py, "best_rep")));
    } else if n_blanks > 0 {
        // blanks make up minority, so filter them out for consensus
        // TODO: filter seqs_l
        n_seqs -= n_blanks;
    }

    let mut seqs_sorted: Vec<&str> = seqs.clone();
    seqs_sorted.sort();

    if seqs_sorted[n_seqs / 2].len() > max_mdn_poa_length {
        return run_best_representatives(py, &seqs_sorted, logger);
    }

    poa_consensus_seq(&seqs_sorted).map_or_else(|| {
        logger.call_method1(
            intern!(py, "error"),
            (intern!(py, "Got no consensus sequence from sequences %s; trying best representative strategy"), seqs_sorted.clone()),
        ).unwrap();
        run_best_representatives(py, &seqs_sorted, logger)
    }, |poa_c| Some((poa_c, intern!(py, "poa"))))
}

#[cfg(test)]
mod test {
    use super::*;

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
