# Phrank Similarity Engine [![License: MIT][license-badge]][license] [![Static Badge](https://img.shields.io/badge/RobinsonGroup-8A2BE2?style=flat&color=Green)](https://robinsongroup.github.io/)

[license]: https://opensource.org/licenses/MIT

[license-badge]: https://img.shields.io/badge/License-MIT-blue.svg
Phrank is a high-performance, phenotype-driven similarity engine designed to calculate the similarity between patient
cohorts. By leveraging information theory, Phrank quantifies the significance of shared phenotypic features (such as
Human Phenotype Ontology terms) based on their Information Content (IC). This allows for the rarity of a shared
phenotype to dictate the weight of the similarity across the cohort.
The algorithm is available in Rust and Python. Check these distribution platforms for examples

- [Pypi](https://pypi.org/project/phrank-py/)
- [Crate.io](https://crates.io/crates/phrank/0.2.1)

# Credit

Original Publication by Karthik A. Jagadeesh et al. [here](https://www.nature.com/articles/s41436-018-0072-y)