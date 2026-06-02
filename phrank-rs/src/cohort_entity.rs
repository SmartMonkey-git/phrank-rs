#[derive(Debug, Hash, PartialEq, Eq, PartialOrd, Ord, Clone)]
pub struct CohortEntity<'a> {
    id: &'a str,
    features: Vec<&'a str>,
}

impl<'a> CohortEntity<'a> {
    pub fn new(id: &'a str, features: Vec<&'a str>) -> Self {
        Self {
            id: id.into(),
            features,
        }
    }
    pub fn id(&self) -> &str {
        self.id
    }
    pub fn features(&self) -> &[&'a str] {
        &self.features
    }
}
