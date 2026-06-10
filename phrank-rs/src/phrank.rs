use crate::cohort_entity::CohortEntity;
use crate::error::PhrankError;
use crate::traits::OntologyTraversal;
use bimap::BiBTreeMap;
use itertools::Itertools;
use rayon::iter::IntoParallelRefIterator;
use rayon::iter::ParallelIterator;
use sprs::TriMat;
use std::collections::{HashMap, HashSet};
use std::marker::PhantomData;

mod sealed {
    pub trait Sealed {}
}

pub trait ICMode: sealed::Sealed {}

pub struct PrecomputedIC;
pub struct RuntimeIC;

impl sealed::Sealed for PrecomputedIC {}
impl sealed::Sealed for RuntimeIC {}

impl ICMode for PrecomputedIC {}
impl ICMode for RuntimeIC {}

/// A structural comparison tool for calculating phenotype-driven similarity
/// across a cohort of patients.
///
/// `Phrank` uses the Information Content (IC) of ontological features
/// to weight the rarity and significance of shared phenotypes.
pub struct Phrank<O, ModeIC: ICMode> {
    ontology: O,
    ic_strategy: PhantomData<ModeIC>,
    ic: HashMap<String, f64>,
}

impl<O> Phrank<O, RuntimeIC>
where
    O: OntologyTraversal,
{
    pub fn new(ontology: O) -> Self {
        Self {
            ontology,
            ic_strategy: PhantomData,
            ic: HashMap::new(),
        }
    }

    pub fn into_precomputed(self, ic: HashMap<String, f64>) -> Phrank<O, PrecomputedIC> {
        Phrank::with_precomputed(self.ontology, ic)
    }

    pub fn calculate_similarity(
        &self,
        cohort: &[CohortEntity],
    ) -> Result<(TriMat<f64>, BiBTreeMap<usize, String>), PhrankError> {
        if cohort.len() <= 2 {
            return Err(PhrankError::CohortTooSmall(2, cohort.len()));
        }
        let ic = Phrank::<O, RuntimeIC>::calculate_ic(cohort, &self.ontology)?;
        self.inner_calculate_similarity(cohort, &ic)
    }
}

impl<O> Phrank<O, PrecomputedIC>
where
    O: OntologyTraversal,
{
    pub fn with_precomputed(ontology: O, ic: HashMap<String, f64>) -> Self {
        Self {
            ontology,
            ic_strategy: PhantomData,
            ic,
        }
    }

    pub fn ic(&self) -> &HashMap<String, f64> {
        &self.ic
    }

    pub fn into_runtime(self) -> Phrank<O, RuntimeIC> {
        Phrank::new(self.ontology)
    }

    pub fn calculate_similarity(
        &self,
        cohort: &[CohortEntity],
    ) -> Result<(TriMat<f64>, BiBTreeMap<usize, String>), PhrankError> {
        if cohort.len() <= 1 {
            return Err(PhrankError::CohortTooSmall(1, cohort.len()));
        }
        self.inner_calculate_similarity(cohort, &self.ic)
    }
}

impl<O, IC> Phrank<O, IC>
where
    O: OntologyTraversal,
    IC: ICMode,
{
    /// Calculates the Information Content (IC) for all terms present in the cohort.
    ///
    /// This method automatically propagates annotations up the ontology tree.
    /// If a patient is annotated with a specific term, they are also considered
    /// annotated with all of its ancestors.
    ///
    /// # Arguments
    /// * `cohort` - A reference to a `HashMap` where the key is a Patient ID
    ///   and the value is a `Vec` of feature/phenotype IDs.
    ///
    /// # Returns
    /// A `Result` containing a `HashMap` mapping phenotype IDs to their respective
    /// IC score as an `f32`, or a `PhrankError` if ontology traversal fails.
    pub fn calculate_ic(
        cohort: &[CohortEntity],
        ontology: &O,
    ) -> Result<HashMap<String, f64>, PhrankError> {
        let cohort_size = cohort.len() as f64;

        let mut direct_associations: HashMap<&str, HashSet<&str>> = HashMap::new();

        for entity in cohort.iter() {
            for feature_id in entity.features() {
                direct_associations
                    .entry(*feature_id)
                    .or_default()
                    .insert(entity.id());
            }
        }

        let mut phenotype_patient_association: HashMap<String, HashSet<&str>> = HashMap::new();

        for (pt_id, patients) in direct_associations {
            let mut ancestors = ontology.get_ancestor_ids(pt_id)?;
            ancestors.push(pt_id.to_owned());

            for ancestor_id in ancestors {
                phenotype_patient_association
                    .entry(ancestor_id)
                    .or_default()
                    .extend(patients.iter().copied());
            }
        }

        Ok(phenotype_patient_association
            .into_iter()
            .map(|(pt_id, patients)| {
                let information_content = -(patients.len() as f64 / cohort_size).log2();
                (pt_id, information_content)
            })
            .collect())
    }

    /// Computes the pairwise similarity matrix for an entire cohort.
    ///
    /// This function performs a parallelized Cartesian product over the cohort.
    /// The similarity between two patients is calculated by summing the
    /// Information Content (IC) of the intersection of their direct feature sets.
    ///
    /// # Arguments
    /// * `cohort` - A reference to a `HashMap` mapping Patient IDs to a list of feature IDs.
    ///
    /// # Returns
    /// A `Result` containing a tuple:
    /// 1. A `TriMat<f32>`: A sparse coordinate matrix containing the similarity scores.
    /// 2. A `BiMap<usize, String>`: A bidirectional map linking the numeric indices of
    ///    the sparse matrix to the original string-based Patient IDs.
    fn inner_calculate_similarity(
        &self,
        cohort: &[CohortEntity],
        ic: &HashMap<String, f64>,
    ) -> Result<(TriMat<f64>, BiBTreeMap<usize, String>), PhrankError> {
        if cohort
            .iter()
            .map(|entity| entity.id())
            .collect::<HashSet<&str>>()
            .len()
            < cohort.len()
        {
            return Err(PhrankError::DuplicateIDs);
        }

        let mut matrix = TriMat::<f64>::new((cohort.len(), cohort.len()));
        let pp_to_matrix_id: BiBTreeMap<usize, String> = cohort
            .iter()
            .enumerate()
            .map(|(idx, entity)| (idx, entity.id().to_owned()))
            .collect();

        let mut closures: Vec<HashSet<String>> = Vec::with_capacity(cohort.len());

        for entity in cohort.iter() {
            let mut closure = HashSet::new();
            for &feature_id in entity.features() {
                closure.insert(feature_id.to_string());

                let ancestors = self.ontology.get_ancestor_ids(feature_id)?;
                for ancestor in ancestors {
                    closure.insert(ancestor);
                }
            }
            closures.push(closure);
        }

        let entities_with_closures: Vec<(&CohortEntity, &HashSet<String>)> =
            cohort.iter().zip(closures.iter()).collect();

        let product: Vec<_> = entities_with_closures
            .iter()
            .cartesian_product(entities_with_closures.iter())
            .collect();

        let results: Vec<(usize, usize, f64)> = product
            .par_iter()
            .map(|&(&(entity_1, closure_1), &(entity_2, closure_2))| {
                let mut similarity = 0.0_f64;

                for key in closure_1.intersection(closure_2) {
                    if let Some(&information_content) = ic.get(key) {
                        similarity += information_content;
                    } else {
                        panic!("Missing Information Content for input {key}.");
                    }
                }

                let row = *pp_to_matrix_id.get_by_right(entity_1.id()).unwrap();
                let col = *pp_to_matrix_id.get_by_right(entity_2.id()).unwrap();
                (row, col, similarity)
            })
            .collect();

        for (row, col, sim) in results {
            matrix.add_triplet(row, col, sim);
        }

        Ok((matrix, pp_to_matrix_id))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
    use sprs::CsMat;
    use std::collections::HashMap;

    #[allow(unused)]
    fn print_patient_grid(trimat: &TriMat<f64>, patient_map: &BiBTreeMap<usize, String>) {
        let (rows, cols) = trimat.shape();
        let csr: CsMat<f64> = trimat.to_csr();

        let get_id = |index: &usize| -> &str {
            patient_map
                .get_by_left(index)
                .map(|s| s.as_str())
                .unwrap_or("???")
        };

        // 1. Print Column Headers
        print!("{:>12} ", ""); // Empty space for the top-left corner
        for c in 0..cols {
            // Truncate column headers to 8 characters to fit the numbers
            print!("{:>8.8} ", get_id(&c));
        }
        println!();

        for r in 0..rows {
            print!("{:>12.12} [", get_id(&r));

            for c in 0..cols {
                if let Some(val) = csr.get(r, c) {
                    print!("{:>8.2} ", val);
                } else {
                    print!("{:>8} ", "·");
                }
            }
            println!("]");
        }
    }

    struct MockOntology {
        ancestor_map: HashMap<String, Vec<String>>,
    }

    impl OntologyTraversal for MockOntology {
        fn get_ancestor_ids(&self, id: &str) -> Result<Vec<String>, PhrankError> {
            Ok(self.ancestor_map.get(id).cloned().unwrap_or_default())
        }
    }

    fn mock_ontology() -> MockOntology {
        let mut ancestor_map = HashMap::new();
        ancestor_map.insert("HP:001".to_string(), vec!["HP:000".to_string()]);
        ancestor_map.insert("HP:002".to_string(), vec!["HP:000".to_string()]);
        ancestor_map.insert(
            "HP:003".to_string(),
            vec!["HP:002".to_string(), "HP:000".to_string()],
        );

        MockOntology { ancestor_map }
    }
    fn setup_mock_phrank() -> Phrank<MockOntology, RuntimeIC> {
        Phrank::<MockOntology, RuntimeIC>::new(mock_ontology())
    }

    #[test]
    fn test_calculate_ic() {
        let cohort = vec![
            CohortEntity::new("P1", vec!["HP:001"]),
            CohortEntity::new("P2", vec!["HP:002"]),
        ];

        let ic_map = Phrank::<MockOntology, RuntimeIC>::calculate_ic(&cohort, &mock_ontology())
            .expect("Failed to calculate IC");

        assert_eq!(ic_map.get("HP:000").copied().unwrap_or(f64::NAN), 0.0);
        assert_eq!(ic_map.get("HP:001").copied().unwrap_or(f64::NAN), 1.0);
        assert_eq!(ic_map.get("HP:002").copied().unwrap_or(f64::NAN), 1.0);
    }

    #[test]
    fn test_calculate_ic_same() {
        let cohort = vec![
            CohortEntity::new("P1", vec!["HP:001"]),
            CohortEntity::new("P2", vec!["HP:001"]),
        ];

        let ic_map = Phrank::<MockOntology, RuntimeIC>::calculate_ic(&cohort, &mock_ontology())
            .expect("Failed to calculate IC");

        assert_eq!(ic_map.get("HP:000").copied().unwrap_or(f64::NAN), 0.0);
        assert_eq!(ic_map.get("HP:001").copied().unwrap_or(f64::NAN), 0.0);
    }

    #[test]
    fn test_no_similarity() {
        let phrank = setup_mock_phrank();

        let cohort = vec![
            CohortEntity::new("P1", vec!["HP:999"]),
            CohortEntity::new("P2", vec!["HP:002"]),
            CohortEntity::new("P3", vec!["HP:003"]),
        ];

        let (matrix, bimap) = phrank
            .calculate_similarity(&cohort)
            .expect("Failed to calculate similarity");

        assert!(bimap.contains_right("P1"));
        assert!(bimap.contains_right("P2"));

        let id_a = *bimap.get_by_right("P1").unwrap();
        let id_b = *bimap.get_by_right("P2").unwrap();

        let csr_matrix: sprs::CsMat<f64> = matrix.to_csr();

        let sim_score = csr_matrix.get(id_a, id_b).copied().unwrap_or(0.0); // Safely unwrap to 0.0 if the sparse matrix naturally omitted it

        assert_eq!(sim_score, 0.0);
    }
    #[test]
    fn test_similarity() {
        let phrank = setup_mock_phrank();

        let cohort = vec![
            CohortEntity::new("P1", vec!["HP:001"]),
            CohortEntity::new("P2", vec!["HP:001"]),
            CohortEntity::new("P3", vec!["HP:002"]),
        ];

        let (matrix, bimap) = phrank
            .calculate_similarity(&cohort)
            .expect("Failed to calculate similarity");

        for (idx, (_matrix_id, pp_id)) in bimap.iter().enumerate() {
            assert_eq!(pp_id, cohort[idx].id())
        }

        assert!(bimap.contains_right("P1"));
        assert!(bimap.contains_right("P2"));

        let p_id_1 = *bimap.get_by_right("P1").unwrap();
        let p_id_2 = *bimap.get_by_right("P2").unwrap();
        let p_id_3 = *bimap.get_by_right("P3").unwrap();

        let csr_matrix: CsMat<f64> = matrix.to_csr();

        let sim_score_p1_p2 = csr_matrix
            .get(p_id_1, p_id_2)
            .copied()
            .expect("Failed to get sim score");

        assert_relative_eq!(sim_score_p1_p2, 0.5849625, epsilon = 1e-5);

        let sim_score_p3_p3 = csr_matrix
            .get(p_id_3, p_id_3)
            .copied()
            .expect("Failed to get sim score");

        assert!(sim_score_p3_p3 > sim_score_p1_p2);
    }

    #[test]
    fn test_matrix_shape_three() {
        let phrank = setup_mock_phrank();

        let cohort = vec![
            CohortEntity::new("P1", vec!["HP:001"]),
            CohortEntity::new("P2", vec!["HP:002"]),
            CohortEntity::new("P3", vec!["HP:002"]),
        ];

        let (matrix, _bimap) = phrank
            .calculate_similarity(&cohort)
            .expect("Failed to calculate similarity");

        assert_eq!(matrix.rows(), 3);
        assert_eq!(matrix.cols(), 3);

        assert_eq!(matrix.col_inds().len(), 9);
        assert_eq!(matrix.row_inds().len(), 9);
    }

    #[test]
    fn test_cohort_too_small() {
        let phrank = setup_mock_phrank();

        let cohort = vec![CohortEntity::new("P1", vec!["HP:001"])];

        let res = phrank.calculate_similarity(&cohort);

        match res.err().unwrap() {
            PhrankError::CohortTooSmall(2, cohort_len) => {
                assert_eq!(cohort_len, 1);
            }
            _ => panic!("Wrong error"),
        }

        let res = phrank.calculate_similarity(vec![].as_slice());

        match res.err().unwrap() {
            PhrankError::CohortTooSmall(2, cohort_len) => {
                assert_eq!(cohort_len, 0);
            }
            _ => panic!("Wrong error"),
        }
    }

    #[test]
    fn test_cohort_duplicate_ids() {
        let phrank = setup_mock_phrank();

        let cohort = vec![
            CohortEntity::new("P1", vec!["HP:001"]),
            CohortEntity::new("P1", vec!["HP:001"]),
            CohortEntity::new("P1", vec!["HP:001"]),
        ];

        let res = phrank.calculate_similarity(&cohort);

        match res.err().unwrap() {
            PhrankError::DuplicateIDs => {}
            _ => panic!("Wrong error"),
        }

        let res = phrank.calculate_similarity(vec![].as_slice());

        match res.err().unwrap() {
            PhrankError::CohortTooSmall(2, cohort_len) => {
                assert_eq!(cohort_len, 0);
            }
            _ => panic!("Wrong error"),
        }
    }
}
