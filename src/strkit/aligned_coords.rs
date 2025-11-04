use bincode;
use numpy::{PyArray1, PyArrayMethods};
use pyo3::prelude::*;
use pyo3::types::PyBytes;
use serde::{Deserialize, Serialize};

#[derive(Clone, Deserialize, Serialize)]
#[pyclass(module = "strkit_rust_ext")]
pub struct STRkitAlignedCoords {
    #[pyo3(get)]
    pub query_coords: Vec<u64>,
    #[pyo3(get)]
    pub ref_coords: Vec<u64>,
}

#[pymethods]
impl STRkitAlignedCoords {
    #[new]
    fn py_new<'py>(query_coords: Option<&Bound<'py, PyArray1<u64>>>, ref_coords: Option<&Bound<'py, PyArray1<u64>>>) -> PyResult<Self> {
        // Initializing with Nones MUST only happen during unpickling - this is ugly, sorry (see __getnewargs__ below).
        Ok(
            STRkitAlignedCoords {
                query_coords: if let Some(qc) = query_coords { qc.to_vec().unwrap() } else { Vec::new() },
                ref_coords: if let Some(rc) = ref_coords { rc.to_vec().unwrap() } else { Vec::new() },
            }
        )
    }

    fn query_coord_at_idx(&self, idx: usize) -> u64 {
        self.query_coords[idx]
    }

    // --- below are functions which make this class pickle-able ---

    pub fn __setstate__(&mut self, state: Bound<'_, PyBytes>) -> PyResult<()> {
        // TODO: replace unwrap with actual PyResult error
        (*self, _) = bincode::serde::decode_from_slice(state.as_bytes(), bincode::config::standard()).unwrap();
        Ok(())
    }

    pub fn __getstate__<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyBytes>> {
        // TODO: replace unwrap with actual PyResult error
        Ok(PyBytes::new(py, &bincode::serde::encode_to_vec(self, bincode::config::standard()).unwrap()))
    }

    // ugly hack with bad typing hygiene - initialize with Python None, then set state using bincode
    pub fn __getnewargs__(&self) -> PyResult<(Option<usize>, Option<usize>)> {
        Ok((None, None))
    }
}
