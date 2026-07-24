use std::collections::HashMap;

use ndarray::Axis;
use numpy::{PyArray, PyArray1, PyArray2, PyArrayMethods};
use pyo3::{exceptions::PyException, prelude::*, types::PyDict};
use smallvec::SmallVec;

pub type SeqAndConsensusMethod = (String, String); // TODO: enum for consensus method

/// Read-allele assignment method: how reads were assigned to particular alleles
#[pyclass(eq, eq_int, from_py_object)]
#[derive(Clone, PartialEq)]
pub enum AssignMethod {
    Dist = 0, // Slight misnomer, should be GMM
    SNV = 1, // Using flanking SNVs
    SNVAndDist = 2, // Using a distance metric combining copy number and SNVs
    Single = 3, // Haploid locus, only one allele
    HP = 4, // Haplotagged reads
}

#[pymethods]
impl AssignMethod {
    pub fn as_str(&self) -> &'static str {
        match self {
            AssignMethod::Dist => "dist",
            AssignMethod::SNV => "snv",
            AssignMethod::SNVAndDist => "snv+dist",
            AssignMethod::Single => "single",
            AssignMethod::HP => "hp",
        }
    }
}

pub struct CallPeaksData {
    pub means: SmallVec<[f64; 2]>,
    pub weights: SmallVec<[f64; 2]>,
    pub stdevs: SmallVec<[f64; 2]>,
    pub modal_n: u8,
    // Needs to be calculated at the end when we do read-peak assignment
    pub n_reads: Option<SmallVec<[u16; 2]>>,
    // if k-mer collection is enabled:
    pub kmers: Option<Vec<HashMap<String, u16>>>,
    // if consensus is enabled:
    pub seqs: Option<Vec<SeqAndConsensusMethod>>,
    pub start_anchor_seqs: Option<Vec<SeqAndConsensusMethod>>,
    // if 5mCpG methylation is enabled:
    pub am: Option<SmallVec<[f64; 2]>>,  // average methylation per peak
}

impl CallPeaksData {
    /// Reverses all peak data as part of a phasing fix-up process.
    pub fn reverse(&mut self) {
        self.means.reverse();
        self.weights.reverse();
        self.stdevs.reverse();
        if let Some(ref mut n_reads) = self.n_reads { n_reads.reverse(); }
        if let Some(ref mut kmers) = self.kmers { kmers.reverse(); }
        if let Some(ref mut seqs) = self.seqs { seqs.reverse(); }
        if let Some(ref mut start_anchor_seqs) = self.start_anchor_seqs { start_anchor_seqs.reverse(); }
        if let Some(ref mut am) = self.am { am.reverse(); }
    }

    pub fn to_dict<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        let res = PyDict::new(py);
        res.set_item("means", self.means.as_slice())?;
        res.set_item("weights", self.weights.as_slice())?;
        res.set_item("stdevs", self.stdevs.as_slice())?;
        res.set_item("modal_n", self.modal_n)?;
        res.set_item("n_reads", self.n_reads.as_ref().map(|nr| nr.as_slice()))?;
        match &self.kmers {
            Some(kmers) => res.set_item("kmers", kmers)?,
            None => {},
        };
        match &self.seqs {
            Some(seqs) => res.set_item("seqs", seqs)?,
            None => {},
        };
        match &self.start_anchor_seqs {
            Some(start_anchor_seqs) => res.set_item("start_anchor_seqs", start_anchor_seqs)?,
            None => {},
        };
        match &self.am {
            Some(am) => res.set_item("am", am.as_slice())?,
            None => {},
        };
        Ok(res)
    }
}

#[pyclass]
pub struct CallData {
    pub call: SmallVec<[i32; 2]>,
    pub call_95_cis: SmallVec<[(i32, i32); 2]>,
    pub call_99_cis: SmallVec<[(i32, i32); 2]>,
    pub peaks: CallPeaksData,
    pub assign_method: Option<AssignMethod>,
    #[pyo3(get)]
    pub ps: Option<i32>,  // phase set
}

fn bound_pyarray_ownedrepr<'py, T: Copy + numpy::Element, D: ndarray::Dimension>(
    x: Bound<'py, PyArray<T, D>>
) -> ndarray::ArrayBase<ndarray::OwnedRepr<T>, D, T> {
    x.readonly().as_array().as_standard_layout().into_owned()
}

fn bound_pyarray_to_smallvec2<'py, T: Copy + numpy::Element>(x: Bound<'py, PyArray1<T>>) -> SmallVec<[T; 2]> {
    bound_pyarray_ownedrepr(x).iter().map(|&i| i).collect()
}

#[pymethods]
impl CallData {
    #[new]
    fn py_new<'py>(
        call: Bound<'py, PyArray1<i32>>,
        call_95_cis: Bound<'py, PyArray2<i32>>,
        call_99_cis: Bound<'py, PyArray2<i32>>,
        means: Bound<'py, PyArray1<f64>>,
        weights: Bound<'py, PyArray1<f64>>,
        stdevs: Bound<'py, PyArray1<f64>>,
        modal_n: u8,
    ) -> PyResult<Self> {
        let c95cis = bound_pyarray_ownedrepr(call_95_cis);
        let c99cis = bound_pyarray_ownedrepr(call_99_cis);
        Ok(CallData {
            call: bound_pyarray_to_smallvec2(call),
            call_95_cis: c95cis.axis_iter(Axis(0)).map(|i| (i[0], i[1])).collect(),
            call_99_cis: c99cis.axis_iter(Axis(0)).map(|i| (i[0], i[1])).collect(),
            peaks: CallPeaksData {
                means: bound_pyarray_to_smallvec2(means),
                weights: bound_pyarray_to_smallvec2(weights),
                stdevs: bound_pyarray_to_smallvec2(stdevs),
                modal_n,
                n_reads: None,
                kmers: None,
                seqs: None,
                start_anchor_seqs: None,
                am: None,
            },
            assign_method: None,
            ps: None,
        })
    }

    #[getter]
    fn call<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<i32>> {
        PyArray1::from_iter(py, self.call.iter().map(|&i| i))
    }

    #[getter]
    fn call_95_cis<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyArray2<i32>>> {
        let vec2: Vec<Vec<i32>> = self.call_95_cis.iter().map(|i| vec![i.0, i.1]).collect();
        PyArray2::from_vec2(py, &vec2).map_err(|e| PyException::new_err(e.to_string()))
    }

    #[getter]
    fn peak_means<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        PyArray1::from_iter(py, self.peaks.means.iter().map(|&i| i))
    }

    #[getter]
    fn peak_weights<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        PyArray1::from_iter(py, self.peaks.weights.iter().map(|&i| i))
    }

    #[getter]
    fn peak_stdevs<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        PyArray1::from_iter(py, self.peaks.stdevs.iter().map(|&i| i))
    }

    #[getter]
    fn peak_modal_n(&self) -> u8 { self.peaks.modal_n }

    fn get_peaks_dict<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        self.peaks.to_dict(py)
    }

    fn get_assign_method_str(&self) -> &'static str {
        self.assign_method.as_ref().map_or("none", |am| am.as_str())
    }

    // --- Peak data setters ---

    fn set_n_reads<'py>(&mut self, n_reads: Bound<'py, PyArray1<u16>>) {
        self.peaks.n_reads = Some(bound_pyarray_to_smallvec2(n_reads));
    }

    fn set_kmers(&mut self, kmers: Option<Vec<HashMap<String, u16>>>) {
        self.peaks.kmers = kmers;
    }

    fn set_seqs(&mut self, seqs: Vec<SeqAndConsensusMethod>, start_anchor_seqs: Vec<SeqAndConsensusMethod>) {
        self.peaks.seqs = Some(seqs);
        self.peaks.start_anchor_seqs = Some(start_anchor_seqs);
    }

    fn set_am<'py>(&mut self, am: Bound<'py, PyArray1<f64>>) {
        self.peaks.am = Some(bound_pyarray_to_smallvec2(am));
    }

    // ---

    fn set_assign_method(&mut self, assign_method: AssignMethod) {
        self.assign_method = Some(assign_method);
    }

    fn set_ps(&mut self, ps: i32) {
        self.ps = Some(ps);
    }

    /// Reverses all call and peak data as part of a phasing fix-up process.
    fn reverse(&mut self) {
        // call
        self.call.reverse();
        self.call_95_cis.reverse();
        self.call_99_cis.reverse();
        // peaks
        self.peaks.reverse();
    }

    pub fn to_dict<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        let res = PyDict::new(py);
        res.set_item("call", self.call.as_slice())?;
        res.set_item("call_95_cis", self.call_95_cis.as_slice())?;
        res.set_item("call_99_cis", self.call_99_cis.as_slice())?;
        res.set_item("peaks", self.peaks.to_dict(py)?)?;
        res.set_item("assign_method", self.assign_method.as_ref().map(|am| am.as_str()))?;
        res.set_item("ps", self.ps)?;
        Ok(res)
    }
}

/// Combines multiple CallData structs together
#[pyfunction]
pub fn combine_call_data<'py>(calls: Vec<Bound<'py, CallData>>) -> PyResult<CallData> {
    let mut call: SmallVec<[i32; 2]> = SmallVec::new();
    let mut call_95_cis: SmallVec<[(i32, i32); 2]> = SmallVec::new();
    let mut call_99_cis: SmallVec<[(i32, i32); 2]> = SmallVec::new();
    let mut assign_method: Option<AssignMethod> = None;

    let mut means: SmallVec<[f64; 2]> = SmallVec::new();
    let mut weights: SmallVec<[f64; 2]> = SmallVec::new();
    let mut stdevs: SmallVec<[f64; 2]> = SmallVec::new();
    let mut n_reads: Option<SmallVec<[u16; 2]>> = None;
    let mut kmers: Option<Vec<HashMap<String, u16>>> = None;
    let mut seqs: Option<Vec<SeqAndConsensusMethod>> = None;
    let mut start_anchor_seqs: Option<Vec<SeqAndConsensusMethod>> = None;
    let mut am: Option<SmallVec<[f64; 2]>> = None;

    for cd in &calls {
        let call_data = cd.borrow();

        call.extend_from_slice(&call_data.call);
        call_95_cis.extend_from_slice(&call_data.call_95_cis);
        call_99_cis.extend_from_slice(&call_data.call_99_cis);

        if let Some(am) = &call_data.assign_method {
            if assign_method.is_none() {
                assign_method = Some(am.clone());
            } else if let Some(am2) = &assign_method && am != am2 {
                return Err(PyException::new_err("mismatch in assignment methods"));
            }
        }

        means.extend_from_slice(&call_data.peaks.means);
        weights.extend_from_slice(&call_data.peaks.weights);
        stdevs.extend_from_slice(&call_data.peaks.stdevs);
        if let Some(nr) = &call_data.peaks.n_reads {
            match &mut n_reads {
                Some(n_reads_inner) => n_reads_inner.extend_from_slice(nr),
                None => n_reads = Some(nr.clone()),
            }
        }
        if let Some(km) = &call_data.peaks.kmers {
            match &mut kmers {
                Some(kmers_inner) => kmers_inner.extend_from_slice(km),
                None => kmers = Some(km.clone()),
            }
        }
        if let Some(sq) = &call_data.peaks.seqs {
            match &mut seqs {
                Some(seqs_inner) => seqs_inner.extend_from_slice(sq),
                None => seqs = Some(sq.clone()),
            }
        }
        if let Some(sas) = &call_data.peaks.start_anchor_seqs {
            match &mut start_anchor_seqs {
                Some(sas_inner) => sas_inner.extend_from_slice(sas),
                None => start_anchor_seqs = Some(sas.clone()),
            }
        }
        if let Some(am_) = &call_data.peaks.am {
            match &mut am {
                Some(am_inner) => am_inner.extend_from_slice(am_),
                None => am = Some(am_.clone()),
            }
        }
    }

    let wsum = weights.iter().sum::<f64>();

    for i in 0..weights.len() {
        weights[i] /= wsum;
    }

    Ok(CallData {
        call,
        call_95_cis,
        call_99_cis,
        peaks: CallPeaksData {
            means,
            weights,
            stdevs,
            modal_n: calls.len() as u8,
            n_reads,
            kmers,
            seqs,
            start_anchor_seqs,
            am,
        },
        assign_method,
        ps: None,
    })
}
