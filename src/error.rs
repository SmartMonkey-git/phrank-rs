use thiserror::Error;


#[derive(Error, Debug)]
pub enum PhrankError{
    #[error("Could not load ontology: {0}")]
    ReadError(String),
    #[error("TermID not found: {0}")]
    TermIdNotFound(String),
}

