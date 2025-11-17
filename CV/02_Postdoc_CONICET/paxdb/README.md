# Amino Acid Usage under Metabolic Constraints (PaxDB project)

This repository contains a cleaned–up, modular reimplementation of my 2013–2015
PaxDB analysis used in:

> Krick T., Verstraete N., Alonso L.G., Shub D.A., Ferreiro D.U., Sanchez I.E.  
> **Amino Acid Metabolism Conflicts with Protein Diversity.**  
> *Molecular Biology and Evolution*, 2014. DOI:10.1093/molbev/msu228

The goal is to reproduce, in a modern and readable way, the main steps of my
original scripts:

- parsing FASTA proteomes for multiple organisms
- mapping PaxDB protein abundance data to sequences
- computing amino acid usage, with and without abundance weighting
- building species × amino acid tables
- computing correlations between amino acid usage and metabolic cost
- providing a clean entry–point script (`scripts/analyze_amino_acids.py`)

The structure is intentionally simple and pipeline-oriented, but close in spirit
to the original files (`amino_acid_count*.py`, `newdef_protein*.py`,
`aa_relationship.py`, etc.).

---

## Directory layout

```text
paxdb/
├── data/
│   ├── fasta/          # input proteome FASTA files (one per species)
│   ├── abundances/     # PaxDB abundance tables (one per species)
│   └── metadata/       # species metadata, amino acid costs, etc.
├── src/
│   ├── __init__.py
│   ├── utils.py
│   ├── fasta_parser.py
│   ├── abundance_loader.py
│   ├── protein.py
│   ├── aa_metrics.py
│   └── relationships.py
├── scripts/
│   └── analyze_amino_acids.py
└── notebooks/
    └── 01_visualize_amino_acid_bias.ipynb   # optional, not required by pipeline
