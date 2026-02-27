use std::cmp;
use std::collections::HashSet;
use strsim::normalized_levenshtein;

fn seqs_summed_levenshtein_similarity(seqs: &[&str]) -> Vec<f64> {
    // TODO: is this right, and can this be accelerated?

    let n_seqs = seqs.len();
    // TODO: shrink similarity_memo
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

/// Best representative of a set of sequences:
/// Slightly different from a true consensus - returns the string with the minimum Levenshtein distance to all other
/// strings for a particular allele. This roughly approximates a true consensus when |seqs| is large. If more than one
/// best representative exist, the first one is returned. If |best| == |seqs| or |best| == 0, None is returned since
/// there is effectively no true consensus.
pub fn best_representative<'a>(seqs: &'a [&str]) -> Option<&'a str> {
    best_representatives(seqs).into_iter().next()
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
}
