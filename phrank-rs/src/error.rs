use thiserror::Error;

#[derive(Error, Debug)]
pub enum PhrankError {
    #[error("Could not load ontology: {0}.")]
    ReadError(String),
    #[error("TermID not found: {0}.")]
    TermIdNotFound(String),
    #[error("Cohort needs to be >= 2. Got {0}.")]
    CohortTooSmall(usize),
    #[error("Found duplicate ID's in cohort.")]
    DuplicateIDs,
}
