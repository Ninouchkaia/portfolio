# Systemic Effects of SARS-CoV-2 on Host Cellular Functions

**Affiliation:** INSERM U1037 – Centre de Recherches en Cancérologie de Toulouse (CRCT)  
**Period:** 2020–2021  
**Publication:** [Communications Physics, Nature, 2021](https://www.nature.com/articles/s42005-021-00582-8)  

---

## Context
SARS-CoV-2 infection perturbs multiple cellular processes beyond its primary respiratory targets. This collaborative project investigated the **system-wide impact** of viral proteins on host molecular networks to understand how infection reshapes fundamental biological functions. The study combined **virus–host interactomics**, **transcriptomics**, and **functional enrichment** analyses to delineate the molecular footprint of COVID-19 at the cellular level.

---

## Objectives
- Identify cellular pathways and biological processes disrupted by SARS-CoV-2 proteins.  
- Integrate host interaction data with transcriptomic responses to infection.  
- Quantify the extent of systemic perturbation across tissues and molecular functions.  

---

## Methods
- **Data integration:** Protein–protein interaction data from mass spectrometry and curated databases.  
- **Functional enrichment:** Gene Ontology and Reactome analyses to map perturbed biological processes.  
- **Network analysis:** Graph-based propagation of perturbation scores across molecular interaction networks.  
- **Validation:** Comparison with transcriptomic responses from infected cell lines and tissues.  

---

## Contributions
- Performed large-scale functional enrichment analyses of viral–host interaction datasets.  
- Quantified pathway-level perturbations to identify globally affected cellular functions.  
- Contributed to the visualization of systemic functional disruption and cross-tissue comparisons.  
- Participated in the interpretation of network propagation results and manuscript preparation.  

---

## Key Skills
Functional genomics · Network propagation · Enrichment analysis · Multi-omics integration · Systems biology  

---

## Reference
*Ghavasieh A.*, et al. *Multiscale statistical physics of the pan-viral interactome unravels the systemic nature of SARS-CoV-2 infections.*  
*Communications Physics (Nature)*, 2021. [DOI:10.1038/s42005-021-00582-8](https://www.nature.com/articles/s42005-021-00582-8)

---

# SARS-CoV-2 Host Interactome – Reproducible Analysis Pipeline

This repository provides a strict and fully reproducible pipeline for computing **viral–host interactors**, **multi-order network propagation**, and **functional enrichment (Reactome / GO)** for SARS-CoV-2 and a panel of human viruses. The workflow combines **Python** (NetworkX, Pandas) and **R** (clusterProfiler, ReactomePA), and is designed for publication-grade systems biology analysis.

---


## 🔧 Project Structure

```

covid_networks/
│
├── data/
│   ├── raw/                        # raw interactomes (nodes.csv + edges.csv per virus)
│   ├── intermediate/
│   │   ├── interactors/            # 1st–5th order interactors per virus
│   │   ├── gene_lists/             # gene lists for enrichment
│   │   └── enrichment/             # clusterProfiler outputs
│   └── results/
│       ├── tables/                 # gene × virus matrices
│       └── figures/                # plots (optional)
│
├── src/
│   ├── config.py                   # paths
│   ├── io_utils.py                 # loaders for nodes/edges
│   ├── network.py                  # multi-order BFS interactors
│   ├── interactors.py              # per-virus interactor extraction
│   ├── tables.py                   # gene lists + matrices
│   ├── enrichment_wrapper.py       # Python → Rscript interface
│   └── logging_utils.py
│
├── r/
│   ├── enrich_reactome_compareCluster.R
│   └── (optional) enrich_go_bp_compareCluster.R
│
└── analysis.py                     # main command-line orchestrator

```

---

## Getting Started

Place each virus’ interactome under:

```

data/raw/Virus_host_interactomes_thresh25/thresh_0.25/<virus_name>/
├── nodes.csv
└── edges.csv

```

`nodes.csv` format:  
```

node_id, gene_symbol, node_type

```
Where:
- `node_type = 0` → viral protein  
- `node_type = 1` → human protein  

---

## Pipeline Usage

All steps rely on the unified CLI:

```

python analysis.py <command> [options]

````

### Compute multi-order interactors (NetworkX)

```bash
python analysis.py compute_interactors --max-order 4
````

Generates:

```
data/intermediate/interactors/<virus_name>/
    direct_interactors.txt
    only_second_range_interactors.txt
    only_Third_range_interactors.txt
    only_Fourth_range_interactors.txt
    only_Fifth_range_interactors.txt
```

---

### Build gene lists for enrichment (R-compatible)

```bash
python analysis.py gene_lists --range-mode direct_and_second
```

Available modes:

* `direct`
* `only_second`
* `only_third`
* `direct_and_second` (recommended)
* `orders_1_to_3`

Output:

```
data/intermediate/gene_lists/genes_list_<range-mode>.txt
```

---

### Functional enrichment (Reactome via Rscript)

```bash
python analysis.py enrich_reactome --range-mode direct_and_second
```

Calls:

```
r/enrich_reactome_compareCluster.R
```

Output:

```
data/intermediate/enrichment/enrichPathway.tsv
```

---

### Build gene × virus binary matrices

```bash
python analysis.py gene_virus_table --range-mode direct_and_second
```

Output:

```
data/results/tables/gene_virus_table_direct_and_second.tsv
```

---

## Example Full Reproducible Workflow

```bash
python analysis.py compute_interactors --max-order 4
python analysis.py gene_lists --range-mode direct_and_second
python analysis.py enrich_reactome --range-mode direct_and_second
python analysis.py gene_virus_table --range-mode direct_and_second
```

---

## Dependencies

**Python**

* networkx
* pandas

**R**

* clusterProfiler
* ReactomePA
* reactome.db
* readr

---

## Notes

* All results are automatically organized under `data/intermediate/` and `data/results/`.
* The pipeline is deterministic, modular, and suitable for publication workflows.
* Additional enrichment scripts (GO BP/CC/MF) can be added under `r/`.

---
