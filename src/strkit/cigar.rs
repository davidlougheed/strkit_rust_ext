use numpy::{ToPyArray, PyArray1};
use pyo3::prelude::*;
use pyo3::types::PyList;

#[pyfunction]
pub fn get_aligned_pair_matches<'py>(
    py: Python<'py>,
    cigar: &Bound<'py, PyList>, 
    query_start: u64, 
    ref_start: u64,
) -> (Bound<'py, PyArray1<u64>>, Bound<'py, PyArray1<u64>>) {
    let mut qi = query_start;
    let mut di = ref_start;

    let mut qi_vec: Vec<u64> = Vec::new();
    let mut di_vec: Vec<u64> = Vec::new();

    for cigar_op in cigar.iter() {
        let dco0 = cigar_op.get_item(0).unwrap().extract::<u64>().unwrap();

        match dco0 {
            0 | 7 | 8 => {  // MATCH | SEQ_MATCH | SEQ_MISMATCH
                let dco1 = cigar_op.get_item(1).unwrap().extract::<u64>().unwrap();
                
                qi_vec.extend(qi..qi+dco1);
                di_vec.extend(di..di+dco1);

                qi += dco1;
                di += dco1;
            },
            1 | 4 => {  // INSERTION | SOFT_CLIPPED
                let dco1 = cigar_op.get_item(1).unwrap().extract::<u64>().unwrap();
                qi += dco1;
            },
            2 | 3 => {  // DELETION | SKIPPED
                let dco1 = cigar_op.get_item(1).unwrap().extract::<u64>().unwrap();
                di += dco1;
            },
            _ => {  // HARD_CLIPPED | PADDING | (unknown cigar op)
                // Do nothing
            }
        }
    }

    (qi_vec.to_pyarray_bound(py), di_vec.to_pyarray_bound(py))
}
