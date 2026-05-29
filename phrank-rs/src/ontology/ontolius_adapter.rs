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
pub struct CachedOntologyAdapter<Ontology> {
    ontology: Ontology,
    cache: Cache<String, Vec<String>>,
}

impl<Ontology> CachedOntologyAdapter<Ontology> {
    pub fn new(ontology: Ontology, cache_size: u64) -> Self
    where
        Ontology: OntologyTraversal,
    {
        Self {
            ontology,
            cache: Cache::new(cache_size),
        }
    }

    pub fn ontology(&self) -> &Ontology
    where
        Ontology: OntologyTraversal,
    {
        &self.ontology
    }
}

impl<Ontology> OntologyTraversal for CachedOntologyAdapter<Ontology>
where
    Ontology: OntologyTraversal,
{
    fn get_ancestor_ids(&self, child: &str) -> Result<Vec<String>, PhrankError> {
        if let Some(ancestors) = self.cache.get(child) {
            return Ok(ancestors);
        }
        let ancestors = self.ontology.get_ancestor_ids(child)?;
        self.cache.insert(child.to_string(), ancestors.clone());
        Ok(ancestors)
    }
}

impl OntologyTraversal for FullCsrOntology {
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
        let term_id = TermId::from_str(child_id);

        match term_id {
            Ok(term) => {
                let ancestors: Vec<String> = self
                    .iter_ancestor_ids(&term)
                    .map(|t_id| t_id.to_string())
                    .collect();

                Ok(ancestors)
            }
            Err(_) => Err(PhrankError::TermIdNotFound(child_id.to_string())),
        }
    }
}

#[cfg(feature = "obo")]
pub mod obo {
    use crate::error::PhrankError;
    use crate::traits::OntologyTraversal;
    use fastobo::ast::{EntityFrame, OboDoc, TermClause};
    use std::collections::{HashMap, HashSet, VecDeque};

    impl OntologyTraversal for OboDoc {
        fn get_ancestor_ids(&self, child_id: &str) -> Result<Vec<String>, PhrankError> {
            let mut parent_map: HashMap<String, Vec<String>> = HashMap::new();
            let mut term_exists = false;

            let entities: &[EntityFrame] = self.as_ref();
            for entity in entities {
                if let EntityFrame::Term(term_frame) = entity {
                    let id = term_frame.id().as_inner().to_string();

                    if id == child_id {
                        term_exists = true;
                    }

                    let mut parents = Vec::new();

                    for line in term_frame.as_ref() {
                        if let TermClause::IsA(parent_ident) = line.as_inner() {
                            parents.push(parent_ident.to_string());
                        }
                    }

                    parent_map.insert(id, parents);
                }
            }

            if !term_exists {
                return Err(PhrankError::TermIdNotFound(child_id.to_string()));
            }

            let mut ancestors = HashSet::new();
            let mut queue = VecDeque::new();

            if let Some(direct_parents) = parent_map.get(child_id) {
                for parent in direct_parents {
                    ancestors.insert(parent.clone());
                    queue.push_back(parent.clone());
                }
            }

            while let Some(current_id) = queue.pop_front() {
                if let Some(parents) = parent_map.get(&current_id) {
                    for parent_id in parents {
                        if ancestors.insert(parent_id.clone()) {
                            queue.push_back(parent_id.clone());
                        }
                    }
                }
            }

            let ancestor_vec: Vec<String> = ancestors.into_iter().collect();
            Ok(ancestor_vec)
        }
    }

    #[cfg(test)]
    mod tests {
        use crate::traits::OntologyTraversal;
        use ontolius::io::OntologyLoaderBuilder;
        use ontolius::ontology::csr::FullCsrOntology;
        use ontology_registry::{
            BioRegistryMetadataProvider, FileSystemOntologyRegistry, FileType, OboLibraryProvider,
            OntologyRegistration, RegistryKey, SupportedOntology, Version,
        };
        use std::env::temp_dir;
        use std::io::BufReader;

        #[test]
        #[ignore]
        fn test_consistency() {
            let tepmdir = temp_dir();

            let registry = FileSystemOntologyRegistry::new(
                tepmdir,
                BioRegistryMetadataProvider::default(),
                OboLibraryProvider::default(),
            );

            let json_hpo_key =
                RegistryKey::new(SupportedOntology::HP, Version::Latest, FileType::Json);

            let json_hpo = registry.register(json_hpo_key).unwrap();
            let loader = OntologyLoaderBuilder::new().obographs_parser().build();
            let ontolius: FullCsrOntology = loader.load_from_read(json_hpo).unwrap();
            let mut json_ancestors = ontolius.get_ancestor_ids("HP:0006803").unwrap();
            json_ancestors.sort();

            let obo_hpo_key =
                RegistryKey::new(SupportedOntology::HP, Version::Latest, FileType::Obo);
            let ontology_path = registry.register(obo_hpo_key).unwrap();
            let mut reader = BufReader::new(ontology_path);
            let obo_doc = fastobo::from_reader(&mut reader).unwrap();

            let mut obo_ancestors = obo_doc.get_ancestor_ids("HP:0006803").unwrap();
            obo_ancestors.sort();

            assert_eq!(json_ancestors, obo_ancestors);
        }
    }
}
