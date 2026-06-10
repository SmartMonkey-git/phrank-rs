use thiserror::Error;

#[derive(Error, Debug)]
pub enum PhrankError {
    #[error("TermID not found: {0}.")]
    TermIdNotFound(String),
    #[error("Cohort needs to be > {0}. Got {1}.")]
    CohortTooSmall(usize, usize),
    #[error("Found duplicate ID's in cohort.")]
    DuplicateIDs,
}
