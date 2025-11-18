
# SARS-CoV-2 – Host Interactome Pipeline (Strict Reproducible Version)

This repository contains a fully structured, reproducible pipeline for computing
viral–host interactors and pathway enrichments using Python + R.

## Steps

### 1. Compute interactors (NetworkX)
```bash
python analysis.py compute_interactors --max-order 4
````

### 2. Build gene lists for enrichment

```bash
python analysis.py gene_lists --range-mode direct_and_second
```

### 3. Enrichment (Reactome via Rscript)

```bash
python analysis.py enrich_reactome --range-mode direct_and_second
```

### 4. Build gene–virus matrices

```bash
python analysis.py gene_virus_table --range-mode direct_and_second
```

All results are written to:

* `data/intermediate/*`
* `data/results/*`
  