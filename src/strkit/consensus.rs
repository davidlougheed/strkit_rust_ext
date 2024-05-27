use pyo3::prelude::pyfunction;
use bio::alignment::pairwise::Scoring;
use bio::alignment::poa::*;
use std::cmp;
use std::collections::HashSet;
use std::panic;
use strsim::normalized_levenshtein;

static GAP_CHAR_ORD: usize  = b'-' as usize;

#[pyfunction]
pub fn best_representatives(seqs: Vec<&str>) -> HashSet<&str> {
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

#[pyfunction]
pub fn consensus_seq(seqs: Vec<&str>) -> Option<String> {
    panic::catch_unwind(|| {

        let scoring = Scoring::new(-1, 0, |a: u8, b: u8| if a == b { 1i32 } else { -1i32 });

        let _seqs: Vec<&[u8]> = seqs.iter().map(|s| s.as_bytes()).collect();

        if _seqs.len() == 1 {
            return Some(_seqs[0].iter().map(|&b| b as char).collect::<String>());
        }

        let mut _seqs_iter: Vec<&[u8]> = _seqs.clone();

        let mut max_len: usize = 0;
        let mut n_seqs: usize = 0;

        _seqs_iter.pop().map(|x| {
            max_len = cmp::max(max_len, x.len());
            n_seqs += 1;

            let mut aligner = Aligner::new(scoring, x);

            for y in _seqs_iter {
                max_len = cmp::max(max_len, y.len());
                aligner.global(y).add_to_graph();
                n_seqs += 1;
            }

            let pretty = aligner
                .alignment()
                .pretty(aligner.consensus().as_slice(), _seqs, aligner.graph(), max_len * 4);

            let pretty_split = pretty.split('\n').skip(1);

            let mut aligned_seqs: Vec<&str> = Vec::with_capacity(n_seqs);
            let mut max_aligned_len: usize = 0;

            for s in pretty_split.filter(|&z| !z.is_empty()) {
                let s_seq = s.split('\t').nth(1).unwrap();
                max_aligned_len = cmp::max(max_aligned_len, s_seq.len());
                aligned_seqs.push(s_seq);
            }

            let mut consensus_chars: Vec<char> = Vec::with_capacity(max_aligned_len);

            for i in 0..max_aligned_len {
                let mut counter = [0usize; 256];
                for &s in aligned_seqs.iter() {
                    counter[s.as_bytes()[i] as usize] += 1;
                }
                let mode_idx = counter
                    .iter()
                    .enumerate()
                    .max_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap_or(cmp::Ordering::Equal))
                    .map(|(idx, _)| idx)
                    .unwrap();

                if mode_idx != 255 && mode_idx != GAP_CHAR_ORD {
                    // Will be 255 if everything remained 0 somehow
                    // We don't want gap characters in the consensus
                    consensus_chars.push(mode_idx as u8 as char);
                }
            }

            consensus_chars.iter().collect::<String>()
        })
    }).ok()?
}
