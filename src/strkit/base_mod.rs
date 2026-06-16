use rust_lapper::{Interval, Lapper};

use crate::coords::QueryCoord;

pub struct BaseModifications {
    pub tree: Lapper<QueryCoord, u8>, // For fast range lookups of base modifications
}

/// Expands MM/ML-tag-formatted modified bases TODO:
pub fn expand_modified_bases(seq: &str, unmod: &[usize], prob: &[u8]) -> BaseModifications {
    if prob.len() == 0 {
        return BaseModifications { tree: Lapper::new(Vec::new()) };
    }

    let mut queue_ptr = 0;
    let mut prob_ptr = 0;
    let mut counter = 0;

    let mut intervals: Vec<Interval<QueryCoord, u8>> = Vec::with_capacity(unmod.len());

    let mut current_unmod_count = unmod[queue_ptr];

    for (ch_idx, ch) in seq.chars().enumerate() {
        if !(ch == 'C' || ch == 'c') {
            continue;
        }

        if counter < current_unmod_count {
            counter += 1;
            continue;
        }

        let coord = ch_idx as QueryCoord;
        let prob_val = prob[prob_ptr];

        intervals.push(Interval { start: coord, stop: coord + 1, val: prob_val });

        // update unmodified base count slice pointer (move onto the next count of unmodified bases)
        queue_ptr += 1;
        prob_ptr += 1;
        current_unmod_count = unmod[queue_ptr];
        counter = 0;
    }

    BaseModifications { tree: Lapper::new(intervals) }
}

#[cfg(test)]
mod test {
    use super::*;
}
