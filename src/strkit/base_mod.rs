use crate::coords::QueryCoord;

pub struct BaseModifications {

}

/// Expands MM/ML-tag-formatted modified bases TODO:
pub fn expand_modified_bases(seq: &str, unmod: &[usize], prob: &[u8]) -> Vec<(QueryCoord, u8)> {
    if prob.len() == 0 {
        return Vec::new();
    }

    let mut queue_ptr = 0;
    let mut prob_ptr = 0;
    let mut counter = 0;

    let mut res = Vec::with_capacity(unmod.len());

    let mut current_unmod_count = unmod[queue_ptr];

    for (ch_idx, ch) in seq.chars().enumerate() {
        if !(ch == 'C' || ch == 'c') {
            continue;
        }

        if counter < current_unmod_count {
            counter += 1;
            continue;
        }

        res.push((ch_idx as QueryCoord, prob[prob_ptr]));

        // update unmodified base count slice pointer (move onto the next count of unmodified bases)
        queue_ptr += 1;
        prob_ptr += 1;
        current_unmod_count = unmod[queue_ptr];
        counter = 0;
    }

    // TODO: interval tree

    res
}

#[cfg(test)]
mod test {
    use super::*;
}
