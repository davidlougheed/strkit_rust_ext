use numpy::{PyArray1, PyArrayMethods};
use pyo3::prelude::*;

#[derive(Clone)]
#[pyclass(frozen)]
pub struct STRkitAlignedCoords {
    #[pyo3(get)]
    pub query_coords: Vec<u64>,
    #[pyo3(get)]
    pub ref_coords: Vec<u64>,
}

#[pymethods]
impl STRkitAlignedCoords {
    #[new]
    fn py_new<'py>(query_coords: &Bound<'py, PyArray1<u64>>, ref_coords: &Bound<'py, PyArray1<u64>>) -> PyResult<Self> {
        Ok(
            STRkitAlignedCoords {
                query_coords: query_coords.to_vec().unwrap(),
                ref_coords: ref_coords.to_vec().unwrap(),
            }
        )
    }

    fn query_coord_at_idx(&self, idx: usize) -> u64 {
        self.query_coords[idx]
    }
}
