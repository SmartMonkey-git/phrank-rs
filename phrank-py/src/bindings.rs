use numpy::ToPyArray;
use ontolius::io::OntologyLoaderBuilder;
use ontolius::ontology::csr::FullCsrOntology;
use phrank::Phrank;
use phrank::ontology::ontolius_adapter::CachedOntologyAdapter;
use pyo3::prelude::*;
use pyo3::types::PyDict;
use std::collections::HashMap;
use std::fs::File;

/// A high-performance, phenotype-driven similarity engine.
///
/// This engine calculates the structural similarity between patient cohorts
/// using Information Content (IC) derived from an underlying ontology. It leverages
/// Rust's parallelization and zero-copy memory transfers to SciPy.
#[pyclass(name = "PyPhrank")]
pub struct PyPhrank {
    inner: Phrank<CachedOntologyAdapter<FullCsrOntology>>,
}

#[pymethods]
impl PyPhrank {
    /// Initialize the Phrank Engine with a specific ontology.
    ///
    /// Args:
    ///     ontology_path (str): The absolute path to the ontology JSON file.
    ///
    /// Returns:
    ///     PhrankEngine: A ready-to-use similarity engine.
    #[new]
    pub fn new(ontology_path: &str, cache_size: u64) -> PyResult<Self> {
        let loader = OntologyLoaderBuilder::new().obographs_parser().build();

        let ontology_file = File::open(ontology_path)?;
        let ontology = loader
            .load_from_read(ontology_file)
            .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(e.to_string()))?;
        let adapter = CachedOntologyAdapter::new(ontology, cache_size);

        let inner = Phrank::new(adapter);

        Ok(Self { inner })
    }

    /// Calculate the pairwise similarity matrix for a patient cohort.
    ///
    /// This method computes the similarity using parallelized Rust operations
    /// and returns a zero-copy SciPy CSR matrix directly to Python memory.
    ///
    /// Args:
    ///     cohort (dict[str, list[str]]): A dictionary mapping Patient IDs
    ///         to a list of phenotype/HPO terms.
    ///
    /// Returns:
    ///     tuple[scipy.sparse.csr_matrix, dict[int, str]]:
    ///         - The N x N similarity matrix.
    ///         - A dictionary mapping matrix row/col indices to the original Patient IDs.
    pub fn calculate_similarity<'py>(
        &self,
        py: Python<'py>, // Inject the Python GIL token
        cohort: HashMap<String, Vec<String>>,
    ) -> PyResult<(Bound<'py, PyAny>, HashMap<usize, String>)> {
        let num_patients = cohort.len();
        let (matrix, bimap) = self
            .inner
            .calculate_similarity(&cohort)
            .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(e.to_string()))?;

        let rows = matrix.row_inds();
        let cols = matrix.col_inds();
        let data = matrix.data();

        let py_rows = rows.to_pyarray(py);
        let py_cols = cols.to_pyarray(py);
        let py_data = data.to_pyarray(py);

        let scipy_sparse = PyModule::import(py, "scipy.sparse")?;

        let kwargs = PyDict::new(py);
        kwargs.set_item("shape", (num_patients, num_patients))?;

        let coo_args = (py_data, (py_rows, py_cols));
        let coo_matrix = scipy_sparse
            .getattr("coo_matrix")?
            .call((coo_args,), Some(&kwargs))?;

        let csr_matrix = coo_matrix.call_method0("tocsr")?;

        let mut id_map = HashMap::new();
        for (idx, patient_id) in bimap.iter() {
            id_map.insert(*idx, patient_id.clone());
        }

        Ok((csr_matrix, id_map))
    }
}

/// Calculate the pairwise similarity matrix for a patient cohort.
///
/// This method computes the similarity using parallelized Rust operations
/// and returns a zero-copy SciPy CSR matrix directly to Python memory.
///
/// Args:
///     cohort (dict[str, list[str]]): A dictionary mapping Patient IDs
///         to a list of phenotype/HPO terms.
///
/// Returns:
///     tuple[scipy.sparse.csr_matrix, dict[int, str]]:
///         - The N x N similarity matrix.
///         - A dictionary mapping matrix row/col indices to the original Patient IDs.
#[pymodule]
fn phrank_py(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyPhrank>()?;
    Ok(())
}
