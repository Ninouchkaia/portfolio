# SARS-CoV-2 Host Interactome – Reproducible Analysis Pipeline

This repository provides a pipeline for computing  
**viral–host interactors**, **multi-order propagation on the human interactome**, and  
**functional enrichment analysis** (Reactome / GO) for SARS-CoV-2 and a panel of human viruses.  

The workflow is implemented in **Python (NetworkX + Pandas)** and **R (clusterProfiler + ReactomePA)**,  
and is suitable for reproducible systems-biology analyses and publication-grade outputs.

---

## 🔧 Project Structure

```

covid_networks/
│
├── data/
│   ├── raw/                     # nodes.csv + edges.csv for each virus
│   ├── intermediate/
│   │   ├── interactors/         # direct + 2nd/3rd/4th/5th order interactors
│   │   ├── gene_lists/          # gene lists used for enrichment (txt)
│   │   └── enrichment/          # raw enrichment TSV
│   └── results/
│       ├── tables/              # binary matrices genes × viruses
│       └── figures/             # heatmaps, plots (optional)
│
├── src/
│   ├── config.py                # project-wide paths
│   ├── io_utils.py              # loading nodes/edges
│   ├── network.py               # BFS-based multi-order interactors
│   ├── interactors.py           # per-virus interactors extraction
│   ├── tables.py                # gene lists / matrices
│   ├── enrichment_wrapper.py    # Python → Rscript interface
│   └── logging_utils.py
│
├── r/
│   └── enrich_reactome_compareCluster.R
│
└── analysis.py                  # main CLI orchestrator

```

---

## 🚀 Getting Started

Before running the pipeline, place your interactome files under:

```

data/raw/Virus_host_interactomes_thresh25/thresh_0.25/<virus_name>/
├── nodes.csv
└── edges.csv

```

Each `nodes.csv` must contain:  
```

node_id, gene_symbol, node_type

```
Where:
- `node_type = 0` → viral protein  
- `node_type = 1` → human protein  

---

## 🧩 Pipeline Steps

All steps are run through the unified orchestrator:

```

python analysis.py <command> [options]

````

---

## 1️⃣ Compute multi-order interactors (NetworkX)

Compute direct interactors, second order, third order, etc., for *all* viruses:

```bash
python analysis.py compute_interactors --max-order 4
````

Results are written to:

```
data/intermediate/interactors/<virus_name>/
    direct_interactors.txt
    only_second_range_interactors.txt
    only_Third_range_interactors.txt
    only_Fourth_range_interactors.txt
```

---

## 2️⃣ Build gene lists for enrichment

These gene lists are formatted exactly for use with `compareCluster()` in R:

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

## 3️⃣ Functional enrichment (Reactome via Rscript)

Runs the clusterProfiler workflow in R automatically:

```bash
python analysis.py enrich_reactome --range-mode direct_and_second
```

This internally generates the gene list (same `range-mode`) then calls:

```
r/enrich_reactome_compareCluster.R
```

Output:

```
data/intermediate/enrichment/enrichPathway.tsv
```

Contains enrichment for all viruses-side-by-side, as in the publication.

---

## 4️⃣ Build gene–virus matrices

Create binary matrices summarizing which human proteins interact with which viruses:

```bash
python analysis.py gene_virus_table --range-mode direct_and_second
```

Output:

```
data/results/tables/gene_virus_table_direct_and_second.tsv
```

Rows = human genes
Columns = viruses
Entries = 1 if gene interacts (directly or at 2nd order), else 0.

---

## 📦 Example Full Workflow (all steps)

```bash
python analysis.py compute_interactors --max-order 4
python analysis.py gene_lists --range-mode direct_and_second
python analysis.py enrich_reactome --range-mode direct_and_second
python analysis.py gene_virus_table --range-mode direct_and_second
```

This produces interactors, gene lists, enrichment tables, and gene–virus matrices.

---

## 🧬 Dependencies

**Python**

* networkx
* pandas

**R**

* clusterProfiler
* ReactomePA
* reactome.db
* readr

---






```python
readme_md = r'''
# SARS-CoV-2 – Host Interactome Pipeline (Strict Reproducible Version)

This repository contains a fully structured, reproducible pipeline for computing
viral–host interactors and pathway enrichments using Python + R.

## Steps

### 1. Compute interactors (NetworkX)
    python analysis.py compute_interactors --max-order 4

### 2. Build gene lists for enrichment
    python analysis.py gene_lists --range-mode direct_and_second

### 3. Enrichment (Reactome via Rscript)
    python analysis.py enrich_reactome --range-mode direct_and_second

### 4. Build gene–virus matrices
    python analysis.py gene_virus_table --range-mode direct_and_second

All results are written to:
- data/intermediate/*
- data/results/*
'''


```
covid_networks/
│
├── data/
│   ├── raw/
│   │   └── Virus_host_interactomes_thresh25/
│   │       └── thresh_0.25/
│   │           └── <virus_name>/
│   │               ├── nodes.csv
│   │               └── edges.csv
│   ├── intermediate/
│   │   ├── interactors/          # fichiers direct/2nd/3rd... par virus
│   │   ├── gene_lists/           # genes_list_*.txt pour R
│   │   └── enrichment/           # TSV issus de compareCluster()
│   └── results/
│       ├── tables/               # matrices virus×gènes, virus×pathways
│       └── figures/              # heatmaps, etc.
│
├── r/
│   ├── enrich_reactome_compareCluster.R
│   └── enrich_go_bp_compareCluster.R
│
├── src/
│   ├── config.py
│   ├── io_utils.py
│   ├── network.py
│   ├── interactors.py
│   ├── tables.py
│   ├── enrichment_wrapper.py
│   └── logging_utils.py
│
├── analysis.py
└── README.md
```