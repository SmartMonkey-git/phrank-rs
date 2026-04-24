#[derive(Debug, Hash, PartialEq, Eq, PartialOrd, Ord, Clone)]
#[cfg_attr(feature = "python", derive(pyo3::FromPyObject))]
pub struct CohortEntity {
    id: String,
    features: Vec<String>,
}

impl CohortEntity {
    pub fn new(id: impl Into<String>, features: Vec<String>) -> Self {
        Self {
            id: id.into(),
            features,
        }
    }
    pub fn id(&self) -> &str {
        &self.id
    }
    pub fn features(&self) -> &[String] {
        &self.features
    }
}
