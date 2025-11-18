# Evolutionary Constraints on Amino Acid Usage in Proteomes

**Affiliation:** INQUIMAE – CONICET, University of Buenos Aires  
**Period:** 2013–2015  
**Publication:** [Molecular Biology and Evolution, 2014](https://doi.org/10.1093/molbev/msu228)  

---

## Context
Protein composition reflects a balance between evolutionary pressures, biochemical stability, and energetic cost.  
This project aimed to understand how **metabolic constraints** influence amino acid usage across proteomes, and how this trade-off affects protein diversity across species.

---

## Objectives
- Quantify amino acid usage across species in relation to biochemical properties and synthesis cost.  
- Evaluate evolutionary trade-offs between energetic efficiency, chemical stability, and proteome diversity.  
- Develop a quantitative framework describing these multi-objective constraints.  

---

## Methods
- **Data acquisition:** Protein abundance data retrieved from **PaxDB** for 17 model organisms.  
- **Processing:** Cleaning, structuring, and cross-validation of abundance-weighted proteome datasets.  
- **Analysis:** Calculation of amino acid usage frequencies weighted by abundance.  
- **Comparative modeling:** Statistical evaluation of correlations between amino acid cost, stability, and frequency.  
- **Visualization:** Multidimensional scaling and regression analyses representing trade-offs between energy and diversity.  

---

## Contributions
- Retrieved and processed proteome-wide abundance data from PaxDB.  
- Quantified amino acid usage patterns and computed cost–diversity correlations.  
- Contributed to statistical modeling and graphical representation of trade-offs.  
- Participated in manuscript preparation and data validation.  

---

## Reference
Krick T., *Verstraete N.*, Alonso L.G., Shub D.A., Ferreiro D.U., Sanchez I.E.  
*Amino Acid Metabolism Conflicts with Protein Diversity.*  
*Molecular Biology and Evolution*, 2014. [DOI:10.1093/molbev/msu228](https://doi.org/10.1093/molbev/msu228)


---



# Amino Acid Usage under Metabolic Constraints (PaxDB project)

This repository contains a cleaned–up, modular reimplementation of my 2013–2015
PaxDB analysis used in:

> Krick T., Verstraete N., Alonso L.G., Shub D.A., Ferreiro D.U., Sanchez I.E.  
> **Amino Acid Metabolism Conflicts with Protein Diversity.**  
> *Molecular Biology and Evolution*, 2014. DOI:10.1093/molbev/msu228

The goal is to reproduce the main steps of the project :

- parsing FASTA proteomes for multiple organisms
- mapping PaxDB protein abundance data to sequences
- computing amino acid usage, with and without abundance weighting
- building species × amino acid tables
- computing correlations between amino acid usage and metabolic cost
- providing a clean entry–point script (`scripts/analyze_amino_acids.py`)

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
