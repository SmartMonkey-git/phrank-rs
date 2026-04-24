# Phrank Python Bindings (phrank_py)

phrank_py provides Python bindings for the Phrank similarity engine. It is a high-performance, phenotype-driven
similarity engine that calculates the structural similarity between patient cohorts using Information Content (IC)
derived from an underlying ontology.

By wrapping the core Rust implementation, phrank_py delivers the performance of Rust's parallelization combined with
zero-copy memory transfers directly into Python's SciPy ecosystem.

# 🚀 Key Features

**High Performance**: Leverages Rust's multithreading to compute pairwise similarity matrices rapidly.

**Zero-Copy SciPy Integration**: Returns the similarity matrix directly as a scipy.sparse.csr_matrix without duplicating
large arrays in memory.

**Phenopacket Compatible**: Easily parses standard Phenopacket JSON files to build feature dictionaries.

# 🛠 Getting Started

Prerequisites
You will need the scipy package installed to handle the sparse matrix output. If you are parsing Phenopackets, ensure
you have the appropriate protobuf/JSON parsers installed.

Example Usage
The following example demonstrates how to load an ontology, parse a directory of Phenopacket JSON files, and compute a
similarity matrix for the entire cohort.

```python
import os
from pathlib import Path
from google.protobuf.json_format import Parse
from phenopackets import Phenopacket
from phrank_py import PyPhrank

# 1. Initialize the Phrank Engine with your ontology JSON
phrank = PyPhrank("./hp.json", cache_size=1500, normalize=false)

# 2. Load your patient cohort (e.g., from a directory of Phenopackets)
pp_dir = Path(os.path.expanduser("./phenopackets"))
cohort: list[Phenopacket] = [
Parse(json_file.read_text(encoding="utf-8"), Phenopacket())
for json_file in pp_dir.glob("*.json")
]

# 3. Map Patient IDs to a list of their phenotypic feature IDs
id_by_feature_id = {
pp.id: list({pt.type.id for pt in pp.phenotypic_features})
for pp in cohort
}

# 4. Calculate the similarity matrix
# Returns a SciPy CSR matrix and a mapping of matrix indices to Patient IDs
matrix, mapping = phrank.calculate_similarity(id_by_feature_id)

print(f"Generated sparse matrix of shape: {matrix.shape}")
```

# Credit

Original Publication by Karthik A. Jagadeesh et al. [here](https://www.nature.com/articles/s41436-018-0072-y)