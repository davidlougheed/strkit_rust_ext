use pyo3::prelude::*;
use pyo3::types::{PyList, PyTuple};

#[pyfunction]
pub fn get_aligned_pair_matches<'py>(
    py: Python<'py>, 
    cigar: &PyList, 
    query_start: usize, 
    ref_start: usize,
) -> &'py PyTuple {
    let mut qi = query_start;
    let mut di = ref_start;

    // let mut res_vec: Vec<&PyTuple> = vec![];
    let mut qi_vec: Vec<usize> = Vec::new();
    let mut di_vec: Vec<usize> = Vec::new();

    for cigar_op in cigar.iter() {
        let dco0 = cigar_op.get_item(0).unwrap().extract::<usize>().unwrap();

        match dco0 {
            0 | 7 | 8 => {  // MATCH | SEQ_MATCH | SEQ_MISMATCH
                let dco1 = cigar_op.get_item(1).unwrap().extract::<usize>().unwrap();
                
                qi_vec.extend((0..dco1).map(|i| qi + i));
                di_vec.extend((0..dco1).map(|i| di + i));

                qi += dco1;
                di += dco1;
            },
            1 | 4 => {  // INSERTION | SOFT_CLIPPED
                let dco1 = cigar_op.get_item(1).unwrap().extract::<usize>().unwrap();
                qi += dco1;
            },
            2 | 3 => {  // DELETION | SKIPPED
                let dco1 = cigar_op.get_item(1).unwrap().extract::<usize>().unwrap();
                di += dco1;
            },
            _ => {  // HARD_CLIPPED | PADDING | (unknown cigar op)
                // Do nothing
            }
        }
    }

    PyTuple::new(py, [PyList::new(py, qi_vec), PyList::new(py, di_vec)])
}
