//! # Phrank Similarity Engine
//!
//! A high-performance, phenotype-driven similarity engine for comparing patient cohorts.
//!
//! This crate implements the Phrank algorithm to calculate the similarity between
//! sets of phenotypic features (such as Human Phenotype Ontology (HPO) terms).
//! It leverages information theory, where the rarity of a shared phenotype—determined
//! by its Information Content (IC) across the cohort—dictates the weight of the similarity.
//!
//! By utilizing sparse matrices (`sprs`), parallel processing (`rayon`), and
//! memory-efficient caching (`moka`), this crate is optimized to compute pairwise
//! similarity matrices for large patient cohorts rapidly.
//!
//! ## Core Architecture
//!
//! The crate is structurally divided into four main modules:
//!
//! * [`phrank`] - The core algorithmic implementation. Contains the `Phrank` struct
//!   responsible for calculating Information Content (IC) and generating the pairwise
//!   similarity matrix.
//! * [`ontology`] - Contains implementations and adapters for traversing ontology graphs.
//!   Includes `OntologyAdapter` for bridging external libraries (like `ontolius`) with
//!   internal caching mechanisms to speed up ancestor lookups.
//! * [`traits`] - Defines the shared abstractions used throughout the crate, such as
//!   `OntologyTraversal`, allowing the core algorithm to remain agnostic to the
//!   underlying ontology data structure.
//! * [`error`] - Centralized error handling types, including `PhrankError`, representing
//!   domain-specific failures like invalid Term IDs or traversal errors.
//!
//! ## Quick Start
//!
//! ```rust,no_run
//! use std::collections::HashMap;
//! use my_crate::phrank::Phrank;
//! // Assuming you have initialized your chosen ontology:
//! use my_crate::ontology::OntologyAdapter;
//! let adapter = OntologyAdapter::new(my_ontology)?;
//! let phrank = Phrank { ontology: adapter };
//!
//! // 1. Define your cohort (Patient ID -> Vec<Phenotype IDs>)
//! let mut cohort = HashMap::new();
//! cohort.insert("Patient_1".to_string(), vec!["HP:0001250".to_string()]);
//! cohort.insert("Patient_2".to_string(), vec!["HP:0001250".to_string(), "HP:0000001".to_string()]);
//!
//! // 2. Calculate the similarity matrix
//! let (matrix, id_map) = phrank.calculate_similarity(&cohort)?;
//! ```

pub mod error;
pub mod ontology;
pub use ontology::ontolius_adapter::OntologyAdapter;
pub mod phrank;
pub use phrank::Phrank;
pub mod traits;
pub(crate) mod types;
