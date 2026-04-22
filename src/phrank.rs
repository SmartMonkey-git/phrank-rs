use crate::error::PhrankError;
use crate::traits::OntologyTraversal;
use bimap::BiMap;
use itertools::Itertools;
use rayon::iter::IntoParallelRefIterator;
use rayon::iter::ParallelIterator;
use sprs::TriMat;
use std::collections::{HashMap, HashSet};

/// A structural comparison tool for calculating phenotype-driven similarity
/// across a cohort of patients.
///
/// `Phrank` uses the Information Content (IC) of ontological features
/// to weight the rarity and significance of shared phenotypes.
pub struct Phrank<O> {
    ontology: O,
}

impl<O> Phrank<O>
where
    O: OntologyTraversal,
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
    fn _calculate_ic(
        &self,
        cohort: &HashMap<String, Vec<String>>,
    ) -> Result<HashMap<String, f32>, PhrankError> {
        let n_patients = cohort.len() as f32;

        let mut direct_associations: HashMap<&String, HashSet<&str>> = HashMap::new();

        for (id, features) in cohort.iter() {
            for feature_id in features.iter() {
                direct_associations
                    .entry(&feature_id)
                    .or_default()
                    .insert(id);
            }
        }

        let mut phenotype_patient_association: HashMap<String, HashSet<&str>> = HashMap::new();

        for (pt_id, patients) in direct_associations {
            let mut ancestors = self.ontology.get_ancestor_ids(pt_id)?;
            ancestors.push(pt_id.clone());

            for ancestor_id in ancestors {
                phenotype_patient_association
                    .entry(ancestor_id)
                    .or_default()
                    .extend(patients.iter().copied()); // Fast union of HashSets
            }
        }

        Ok(phenotype_patient_association
            .into_iter()
            .map(|(pt_id, patients)| {
                let information_content = -1.0 * (patients.len() as f32 / n_patients).log2();
                (pt_id, information_content)
            })
            .collect())
    }

    /// Computes the pairwise similarity matrix for an entire cohort.
    ///
    /// This function performs a parallelized Cartesian product over the cohort.
    /// The similarity between two patients is calculated by summing the
    /// Information Content (IC) of the union of their direct feature sets.
    ///
    /// # Arguments
    /// * `cohort` - A reference to a `HashMap` mapping Patient IDs to a list of feature IDs.
    ///
    /// # Returns
    /// A `Result` containing a tuple:
    /// 1. A `TriMat<f32>`: A sparse coordinate matrix containing the similarity scores.
    /// 2. A `BiMap<usize, String>`: A bidirectional map linking the numeric indices of
    ///    the sparse matrix to the original string-based Patient IDs.
    pub fn calculate_similarity(
        &self,
        cohort: &HashMap<String, Vec<String>>,
    ) -> Result<(TriMat<f32>, BiMap<usize, String>), PhrankError> {
        let ic = self._calculate_ic(cohort)?;

        let mut matrix = TriMat::new((cohort.len(), cohort.len()));
        let pp_to_matrix_id: BiMap<usize, String> = cohort
            .iter()
            .enumerate()
            .map(|(idx, (id, _))| (idx, id.clone()))
            .collect();

        let product: Vec<_> = cohort
            .iter()
            .cartesian_product(cohort.iter())
            .filter(|(a, b)| !(*a == *b))
            .collect();

        let results: Vec<(usize, usize, f32)> = product
            .par_iter()
            .map(
                |((id_1, features_1), (id_2, features_2)): &(
                    (&String, &Vec<String>),
                    (&String, &Vec<String>),
                )| {
                    let mut similarity = 0.0;
                    for key in
                        HashSet::<&String>::from_iter(features_1.iter().chain(features_2.iter()))
                    {
                        similarity += ic.get(key).unwrap_or(&0.0);
                    }

                    let row = *pp_to_matrix_id.get_by_right(id_1.as_str()).unwrap();
                    let col = *pp_to_matrix_id.get_by_right(id_2.as_str()).unwrap();
                    (row, col, similarity)
                },
            )
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
    use std::collections::HashMap;

    struct MockOntology {
        ancestor_map: HashMap<String, Vec<String>>,
    }

    impl OntologyTraversal for MockOntology {
        fn get_ancestor_ids(&self, id: &str) -> Result<Vec<String>, PhrankError> {
            Ok(self.ancestor_map.get(id).cloned().unwrap_or_default())
        }
    }

    fn setup_mock_phrank() -> Phrank<MockOntology> {
        let mut ancestor_map = HashMap::new();
        ancestor_map.insert("HP:001".to_string(), vec!["HP:000".to_string()]);
        ancestor_map.insert("HP:002".to_string(), vec!["HP:000".to_string()]);

        let ontology = MockOntology { ancestor_map };
        Phrank { ontology }
    }

    #[test]
    fn test_calculate_ic() {
        let phrank = setup_mock_phrank();
        let mut cohort = HashMap::new();
        cohort.insert("Patient_A".to_string(), vec!["HP:001".to_string()]);
        cohort.insert("Patient_B".to_string(), vec!["HP:002".to_string()]);

        let ic_map = phrank
            ._calculate_ic(&cohort)
            .expect("Failed to calculate IC");

        assert_eq!(ic_map.get("HP:000").copied().unwrap_or(f32::NAN), 0.0);
        assert_eq!(ic_map.get("HP:001").copied().unwrap_or(f32::NAN), 1.0);
        assert_eq!(ic_map.get("HP:002").copied().unwrap_or(f32::NAN), 1.0);
    }

    #[test]
    fn test_calculate_similarity() {
        let phrank = setup_mock_phrank();
        let mut cohort = HashMap::new();

        cohort.insert("Patient_A".to_string(), vec!["HP:001".to_string()]);
        cohort.insert("Patient_B".to_string(), vec!["HP:002".to_string()]);

        let (matrix, bimap) = phrank
            .calculate_similarity(&cohort)
            .expect("Failed to calculate similarity");

        assert!(bimap.contains_right("Patient_A"));
        assert!(bimap.contains_right("Patient_B"));

        let id_a = *bimap.get_by_right("Patient_A").unwrap();
        let id_b = *bimap.get_by_right("Patient_B").unwrap();

        assert_eq!(matrix.rows(), 2);
        assert_eq!(matrix.cols(), 2);

        let csr_matrix: sprs::CsMat<f32> = matrix.to_csr();

        let sim_score = csr_matrix.get(id_a, id_b).copied().unwrap_or(0.0);
        assert_eq!(sim_score, 2.0);
    }
}
