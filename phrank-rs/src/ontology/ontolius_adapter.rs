use crate::error::PhrankError;
use crate::traits::OntologyTraversal;
use moka::sync::Cache;
use ontolius::TermId;
use ontolius::ontology::HierarchyWalks;
use ontolius::ontology::csr::FullCsrOntology;
use std::str::FromStr;

/// An adapter that bridges a concrete ontology implementation with the
/// `OntologyTraversal` trait required by the Phrank algorithm.
///
/// This struct includes an embedded, thread-safe cache to store the results
/// of expensive ancestor lookups, significantly speeding up repetitive queries
/// across large patient cohorts.
pub struct OntologyAdapter<Ontology> {
    ontology: Ontology,
    query_cache: Cache<String, Vec<String>>,
}

impl<Ontology> OntologyAdapter<Ontology> {
    pub fn new(ontology: Ontology) -> Self {
        Self {
            ontology,
            query_cache: Cache::new(250),
        }
    }
}

impl OntologyTraversal for OntologyAdapter<FullCsrOntology> {
    /// Retrieves all ancestor IDs for a given child term ID.
    ///
    /// This method first checks the internal cache. If the ancestors for the requested
    /// ID have been queried recently, it returns the cached result immediately.
    /// Otherwise, it parses the term ID, traverses the `FullCsrOntology` graph to
    /// gather all ancestors, caches the result, and then returns it.
    ///
    /// # Arguments
    /// * `child_id` - A string slice representing the phenotype ID to look up.
    ///
    /// # Returns
    /// A `Result` containing a `Vec` of ancestor IDs as `String`s, or a `PhrankError`
    /// if the `child_id` cannot be parsed into a valid `TermId`.
    fn get_ancestor_ids(&self, child_id: &str) -> Result<Vec<String>, PhrankError> {
        if let Some(ancestors) = self.query_cache.get(child_id) {
            return Ok(ancestors);
        }

        let term_id = TermId::from_str(child_id);

        match term_id {
            Ok(term) => {
                let ancestors: Vec<String> = self
                    .ontology
                    .iter_ancestor_ids(&term)
                    .map(|t_id| t_id.to_string())
                    .collect();

                self.query_cache
                    .insert(child_id.to_string(), ancestors.clone());
                Ok(ancestors)
            }
            Err(_) => Err(PhrankError::TermIdNotFound(child_id.to_string())),
        }
    }
}
