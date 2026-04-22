use crate::error::PhrankError;

pub trait OntologyTraversal {
    fn get_ancestor_ids(&self, child: &str) -> Result<Vec<String>, PhrankError>;
}
