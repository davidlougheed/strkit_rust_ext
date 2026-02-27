// use abpoa_rs;
use bio::alignment::pairwise::Scoring;
use bio::alignment::poa::*;
// use poasta::aligner::PoastaAligner;
// use poasta::aligner::config::AffineMinGapCost;
// use poasta::aligner::scoring::{AlignmentType, GapAffine};
// use poasta::graphs::poa::POAGraph;
// use poasta::io::fasta::poa_graph_to_fasta;
use std::cmp;
use std::panic;

static GAP_CHAR_ORD: usize  = b'-' as usize;

fn majority_consensus_from_msa(aligned_seqs: &[&[u8]], aligned_len: usize) -> String {
    (0..aligned_len).filter_map(|i| {
        let mut counter = [0usize; 256];
        aligned_seqs.iter().for_each(|s| {
            counter[s[i] as usize] += 1;
        });

        counter
            .iter()
            .enumerate()
            .max_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap_or(cmp::Ordering::Equal))
            .map(|(idx, _)| idx)
            .and_then(|mode_idx| (mode_idx != 255 && mode_idx != GAP_CHAR_ORD).then_some(mode_idx as u8 as char))
    }).collect::<String>()
}

pub fn poa_consensus_seq(seqs: &[&str]) -> Option<String> {
    let n_seqs = seqs.len();

    if n_seqs == 0 {
        return None;
    }

    let _seqs: Vec<&[u8]> = seqs.iter().map(|s| s.as_bytes()).collect();

    let first_seq = _seqs[0];

    if n_seqs == 1 {
        return Some(first_seq.iter().map(|&b| b as char).collect::<String>());
    }

    let scoring = Scoring::new(-1, 0, |a: u8, b: u8| if a == b { 1i32 } else { -1i32 });

    panic::catch_unwind(|| {
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

// fn poa_consensus_seq_2(seqs: &[&str]) -> Option<String> {
//     let aln_params = abpoa_rs::AlignmentParametersBuilder::new()
//         .alignment_mode(abpoa_rs::AlignmentMode::Global)
//         .gap_affine_penalties(1, 1, 1, 0)
//         .verbosity(abpoa_rs::Verbosity::Debug)
//         .build();

//     let mut graph = abpoa_rs::Graph::new(&aln_params);

//     for (i, &seq) in seqs.iter().enumerate() {
//         eprintln!("Sequence {}: {}", i + 1, seq);

//         let weights = vec![1; seq.len()];
//         let result = graph
//             .align_and_add_sequence(&aln_params, seq.as_bytes(), &weights, format!("seq{}", i).as_bytes())
//             .unwrap();

//         eprintln!("Sequence {}: score = {}", i + 1, result.get_best_score());
//     }

//     graph.generate_consensus(abpoa_rs::ConsensusAlgorithm::HeaviestBundle);

//     let cons = graph.get_consensus().unwrap();
//     let cons_string = String::from_utf8(aln_params.reverse_seq(cons.sequences().iter().next().unwrap())).unwrap();

//     eprintln!("Consensus: {}", cons_string);

//     Some(cons_string)
// }

// fn poa_consensus_seq_3(seqs: &[&str]) -> Option<String> {
//     let n_seqs = seqs.len();

//     if n_seqs == 0 {
//         return None;
//     }

//     let _seqs: Vec<&[u8]> = seqs.iter().map(|s| s.as_bytes()).collect();
//     let first_seq = _seqs[0];

//     if n_seqs == 1 {
//         return Some(first_seq.iter().map(|&b| b as char).collect::<String>());
//     }

//     panic::catch_unwind(|| {
//         let mut graph = POAGraph::<u32>::new();
//         let scoring = GapAffine::new(5, 1, 5);
//         let aligner = PoastaAligner::new(AffineMinGapCost(scoring), AlignmentType::Global);

//         let weights: Vec<usize> = vec![1; first_seq.len()];
//         graph.add_alignment_with_weights("s0", first_seq, None, &weights).unwrap();

//         for (i, &seq) in _seqs[1..].iter().enumerate() {
//             let result = aligner.align::<u32, POAGraph>(&graph, seq);
//             let w = vec![1; seq.len()];
//             graph.add_alignment_with_weights(
//                 format!("s{}", i).as_str(),
//                 seq,
//                 Some(&result.alignment),
//                  &w,
//             ).unwrap();
//         }

//         let mut buf = BufWriter::new(Vec::new());
//         poa_graph_to_fasta(&graph, &mut buf).unwrap();
//         let fasta_alignment = String::from_utf8(buf.into_inner().unwrap()).unwrap();

//         // eprintln!("{}", fasta_alignment);

//         let mut aligned_seqs: Vec<String> = Vec::with_capacity(n_seqs);
//         let mut max_aligned_len: usize = 0;

//         let mut iter = fasta_alignment.split("\n").into_iter();
//         let mut curr = iter.next();

//         // Hack: re-parse FASTA
//         while let Some(line) = curr {
//             if line.starts_with(">") {
//                 curr = iter.next();
//                 if let Some(last_seq) = aligned_seqs.last() {
//                     max_aligned_len = cmp::max(max_aligned_len, last_seq.len());
//                 }
//                 aligned_seqs.push(String::new());
//                 continue;
//             }

//             aligned_seqs.last_mut().unwrap().push_str(line);
//             curr = iter.next();
//         }

//         let aligned_seqs_bytes: Vec<Vec<u8>> = aligned_seqs.into_iter().map(|s| s.into_bytes()).collect();
//         Some(majority_consensus_from_msa(aligned_seqs_bytes, max_aligned_len))

//     }).ok()?
// }

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_poa_consensus_seq() {
        assert_eq!(poa_consensus_seq(&vec![]), None);
        assert_eq!(poa_consensus_seq(&vec!["AA", "AB", "AA"]), Some(String::from("AA")));
    }
}
