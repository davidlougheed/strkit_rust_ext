use rust_htslib::bam::Record;
use rust_lapper::{Interval, Lapper};

use crate::coords::QueryCoord;

pub struct BaseModifications {
    pub tree: Lapper<QueryCoord, u8>, // For fast range lookups of base modifications
}

/// TODO:
pub fn build_modified_bases_tree(record: &Record) -> BaseModifications {
    let mut intervals: Vec<Interval<QueryCoord, u8>> = Vec::new();

    if let Ok(mods) = record.basemods_iter() {
        for res in mods {
            if let Ok((position, m)) = res {
                if m.modified_base as u8 as char == 'm' {
                    intervals.push(
                        Interval {
                            start: position as QueryCoord,
                            stop: (position + 1) as QueryCoord,
                            val: m.qual as u8,
                        }
                    );
                }
            }
        }
    }

    BaseModifications { tree: Lapper::new(intervals) }
}

#[cfg(test)]
mod test {
    use super::*;
}
