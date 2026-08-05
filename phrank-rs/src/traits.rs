use crate::error::PhrankError;
use std::collections::HashSet;

pub trait OntologyTraversal {
    fn get_ancestor_ids(&self, child: &str) -> Result<HashSet<String>, PhrankError>;
}
