use pyo3::prelude::pyfunction;
use bio::alignment::pairwise::Scoring;
use bio::alignment::poa::*;

#[pyfunction]
pub fn consensus_seq(seqs: Vec<&str>) -> Option<String> {
    let scoring = Scoring::new(-1, 0, |a: u8, b: u8| if a == b { 1i32 } else { -1i32 });
    let mut _seqs = seqs.clone();
    _seqs.pop().and_then(|x| {
        let mut aligner = Aligner::new(scoring, x.as_bytes());

        for y in _seqs {
            aligner.global(y.as_bytes()).add_to_graph();
        }

        Some(String::from_utf8(aligner.consensus()).unwrap())
    })
}
