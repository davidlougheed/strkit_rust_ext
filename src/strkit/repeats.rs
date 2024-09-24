use once_cell::sync::Lazy;
use parasail_rs::{Aligner, Matrix, Profile};
use pyo3::pyfunction;
use std::{borrow::Borrow, cmp, collections::{HashMap, HashSet}};

const MATCH_SCORE: i32 = 2;
const MISMATCH_SCORE: i32 = -7;
const INDEL_PENALTY: i32 = 5;

const DNA_BASES: &[u8; 16] = b"ACGTRYSWKMBDHVNX";

static DNA_MATRIX: Lazy<Matrix> = Lazy::new(|| {
    let mut matrix = Matrix::create(DNA_BASES, MATCH_SCORE, MISMATCH_SCORE).unwrap();

    let dna_bases_map: HashMap<u8, i32> = DNA_BASES.iter().enumerate().map(|(i, &b)| (b, i as i32)).collect();

    // Build a hashmap of IUPAC codes and which bases they represent
    let mut dna_codes: HashMap<u8, Vec<u8>> = HashMap::new();
    dna_codes.insert(b'R', vec![b'A', b'G']);
    dna_codes.insert(b'Y', vec![b'C', b'T']);
    dna_codes.insert(b'S', vec![b'G', b'C']);
    dna_codes.insert(b'W', vec![b'A', b'T']);
    dna_codes.insert(b'K', vec![b'G', b'T']);
    dna_codes.insert(b'M', vec![b'A', b'C']);
    dna_codes.insert(b'B', vec![b'C', b'G', b'T']);
    dna_codes.insert(b'D', vec![b'A', b'C', b'T']);
    dna_codes.insert(b'H', vec![b'A', b'C', b'T']);
    dna_codes.insert(b'V', vec![b'A', b'C', b'G']);
    dna_codes.insert(b'N', vec![b'A', b'C', b'G', b'T']);

    // Special STRkit non-IUPAC character for matching low-quality bases
    dna_codes.insert(b'X', vec![b'A', b'C', b'G', b'T']);

    dna_codes.iter().for_each(|(&code, code_matches)| {
        code_matches.iter().for_each(|cm| {
            let value = if code != b'X' { 2 } else { 0 };
            matrix.set_value(dna_bases_map[&code], dna_bases_map[cm], value).unwrap();
            matrix.set_value(dna_bases_map[cm], dna_bases_map[&code], value).unwrap();
        });
    });

    matrix
});


fn score_candidate(aligner: &Aligner, motif: &str, motif_count: usize, flank_left_seq: &str, flank_right_seq: &str) -> i32 {
    let mut candidate = flank_left_seq.to_owned();
    let rep = motif.repeat(motif_count);
    candidate.push_str(rep.as_str());
    candidate.push_str(flank_right_seq);

    aligner.align(None, candidate.as_bytes()).unwrap().get_score()
}


#[pyfunction]
pub fn get_repeat_count(
    start_count: i32,
    tr_seq: &str,
    flank_left_seq: &str,
    flank_right_seq: &str,
    motif: &str,
    max_iters: usize,
    local_search_range: i32,
    step_size: i32,
) -> ((i32, i32), usize, i32) {
    let mut db_seq = flank_left_seq.to_owned();
    db_seq.push_str(tr_seq);
    db_seq.push_str(flank_right_seq);

    let db_seq_profile = Profile::new(
        db_seq.as_bytes(),
        false,
        DNA_MATRIX.borrow()
    ).unwrap();

    let aligner = Aligner::new()
        .gap_open(INDEL_PENALTY)
        .gap_extend(INDEL_PENALTY)
        .semi_global()
        .striped()
        .profile(db_seq_profile)
        .build();

    let max_init_score = (motif.len() as i32 * start_count + flank_left_seq.len() as i32 + flank_right_seq.len() as i32) * MATCH_SCORE;
    let start_score = score_candidate(&aligner, motif, start_count as usize, flank_left_seq, flank_right_seq);

    let score_diff: f64 = (start_score - max_init_score).abs() as f64 / max_init_score as f64;

    let mut lsr = local_search_range;
    let mut step = step_size;

    if score_diff < 0.05 {  // TODO: parametrize
        // If we're very close to the maximum, explore less.
        lsr = 1;
        step = 1;
    } else if score_diff < 0.1 && lsr > 2 {
        lsr = 2;
        step = 1;
    }

    let mut explored_sizes: HashSet<i32> = HashSet::new();
    explored_sizes.insert(start_count);

    let mut best_size: i32 = start_count;
    let mut best_score: i32 = start_score;
    let mut n_explored: usize = 1;
    let mut to_explore: Vec<(i32, bool)> = vec![(start_count - 1, false), (start_count + 1, true)];  // false = left, true = right

    while !to_explore.is_empty() && n_explored < max_iters {
        let (size_to_explore, going_right) = to_explore.pop().unwrap();

        if size_to_explore < 0 {
            continue;
        }

        let skip_search = step > lsr;  // whether we're skipping small areas for a faster search

        let mut best_size_this_round: Option<i32> = None;
        let mut best_score_this_round: i32 = -9999999;

        let start_size = cmp::max(size_to_explore - (if !going_right || skip_search {lsr} else {0}), 0);
        let end_size = size_to_explore + (if going_right || skip_search {lsr} else {0});

        (start_size..end_size+1).for_each(|i| {
            if !explored_sizes.contains(&i) {
                // Generate a candidate TR tract by copying the provided motif 'i' times & score it
                // Separate this from the .get() to postpone computation to until we need it

                explored_sizes.insert(i);
                let i_score = score_candidate(&aligner, motif, i as usize, flank_left_seq, flank_right_seq);

                if best_size_this_round.is_none() || i_score > best_score_this_round {
                    best_size_this_round = Some(i);
                    best_score_this_round = i_score;
                }

                n_explored += 1;
            }
        });

        if let Some(bstr) = best_size_this_round {
            // If this round is the best we've got so far, update the record size/score for the final return

            if best_score_this_round > best_score {
                best_size = bstr;
                best_score = best_score_this_round;

                if lsr > 1 && ((best_score - max_init_score).abs() as f64 / max_init_score as f64) < 0.05 {
                    // reduce search range as we approach an optimum
                    lsr = 1;
                }
            }

            let mr = bstr + step;
            let ml = bstr - step;

            if bstr > size_to_explore && !explored_sizes.contains(&mr) {
                if mr >= 0 {
                    to_explore.push((mr, true));
                }
            } else if bstr < size_to_explore && !explored_sizes.contains(&ml) && ml >= 0 {
                to_explore.push((ml, false));
            }
        }
    }

    ((best_size, best_score), n_explored, best_size - start_count)
}
