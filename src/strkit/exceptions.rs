use pyo3::exceptions::PyException;
use pyo3::prelude::*;

#[pyclass(extends=PyException)]
pub struct LowMeanBaseQual {
    #[pyo3(get)]
    pub mean_base_qual: f32,
}

impl<'py> IntoPyObject<'py> for LowMeanBaseQual {
    // https://stackoverflow.com/questions/79299476/return-a-custom-error-including-payload-in-pyo3

    type Target = PyAny;
    type Output = Bound<'py, Self::Target>;
    type Error = PyErr;

    fn into_pyobject(self, py: Python<'py>) -> Result<Self::Output, Self::Error> {
        py.get_type::<LowMeanBaseQual>().call1((self.mean_base_qual,))
    }
}

#[pymethods]
impl LowMeanBaseQual {
    #[new]
    fn new(mean_base_qual: f32) -> Self {
        Self { mean_base_qual }
    }

    fn __str__(&self) -> String {
        format!("LowMeanBaseQual ({})", self.mean_base_qual)
    }
}
