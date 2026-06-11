# Phrank Similarity Engine

Phrank is a high-performance, phenotype-driven similarity engine designed to calculate the similarity between patient
cohorts. By leveraging information theory, Phrank quantifies the significance of shared phenotypic features (such as
Human Phenotype Ontology terms) based on their Information Content (IC). This allows for the rarity of a shared
phenotype to dictate the weight of the similarity across the cohort.

Built with speed and scalability in mind, this crate utilizes parallel processing, sparse matrix representation, and
efficient caching to compute pairwise similarity matrices for large patient cohorts rapidly.

# 🚀 Key Features

Information-Theory Driven: Uses Information Content (IC) to weight rare phenotypes higher than common ones. The
algorithm automatically propagates annotations up the ontology tree.

**High Performance**:
**Parallelism**: Utilizes rayon to perform parallelized Cartesian product calculations across the cohort.

**Efficient Memory Usage**: Employs sprs for sparse matrix storage to generate the coordinate matrix of similarity
scores.

**Smart Caching**: Uses moka to cache expensive ancestor lookups in a thread-safe manner, significantly reducing
redundant ontology traversals.

**Extensible Architecture**: Designed around the OntologyTraversal trait, allowing you to plug in different ontology
backends seamlessly. It natively includes an adapter for ontolius.

# Credit

Original Publication by Karthik A. Jagadeesh et al. [here](https://www.nature.com/articles/s41436-018-0072-y)