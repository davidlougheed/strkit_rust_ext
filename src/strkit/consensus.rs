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
pub fn best_representative(seqs: Vec<&str>) -> Option<&str> {
    let mut res: Vec<&str> = Vec::from_iter(best_representatives(seqs));
    res.pop()
}


#[pyfunction]
pub fn consensus_seq(seqs: Vec<&str>) -> Option<String> {
    let n_seqs = seqs.len();

    if n_seqs == 0 {
        return None;
    }

    panic::catch_unwind(|| {

        let scoring = Scoring::new(-1, 0, |a: u8, b: u8| if a == b { 1i32 } else { -1i32 });

        let _seqs: Vec<&[u8]> = seqs.into_iter().map(|s| s.as_bytes()).collect();

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
