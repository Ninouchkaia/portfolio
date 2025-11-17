# Nina Verstraete – Bioinformatics & Computational Biology

Computational biologist with experience in:
- End-to-end pipelines (Snakemake, Python, R)
- Bulk RNA-seq deconvolution and immune profiling
- High-throughput barcode data, QC and differential abundance
- Network medicine and multilayer graphs
- Agent-based modelling integrated with reproducible workflows

Experience in molecular and cell biology : here LINK

---

## Quick links for bioinformatics projects

Analysis and pipeline modules, with input/output examples and documentation.

- **Bulk RNA-seq deconvolution pipeline (lung cancer immunotherapy)**  
  Tumour microenvironment deconvolution using multiple published signatures, QC, and prediction models.  
  → `Snakemake + R + Python`, modular `pipeline/` folder.  
  👉 [CV/04_Postdoc_INSERM/03_Deconvolution](04_Postdoc_INSERM/03_Deconvolution)

- **Barcode-based clonal dynamics (high-throughput counts → QC → DESeq2 → networks)**  
  Data cleaning, QC of technical controls, construction of DESeq2 inputs, and drug–drug correlation networks.  
  → `Python + R (DESeq2)`, structured `data/`, `scripts/`, `results/`.  
  👉 [CV/04_Postdoc_INSERM/07_Barcodes](04_Postdoc_INSERM/07_Barcodes)

- **COVID-19 multilayer network medicine pipeline**  
  Construction of multilayer networks (protein–GO, protein–drug, etc.), bootstrap random networks, Z-scores and ranking of targets.  
  → `Python`, `01_multilayer_pipeline/` + `02_bootstrap_pipeline/` clean scripts.  
  👉 [CV/04_Postdoc_INSERM/05_covid_network_medicine](04_Postdoc_INSERM/05_covid_network_medicine)

- **Patient-specific agent-based tumour model: analysis pipeline**  
  Python toolkit around a NetLogo/OpenMOLE ABM, for parameter exploration (NSGA-II), model validation, sensitivity analysis and advanced statistics.  
  → `Python package structure`, `abm_pipeline/` with `cli.py`.  
  👉 [CV/04_Postdoc_INSERM/01/abm_pipeline](04_Postdoc_INSERM/01/abm_pipeline)

- **Amino-acid usage evolution with PaxDB**  
  Modular Python library to parse FASTA, load PaxDB abundances, compute amino-acid metrics and explore correlations.  
  → `src/` package (`fasta_parser.py`, `abundance_loader.py`, `aa_metrics.py`, `relationships.py`).  
  👉 [CV/02_Postdoc_CONICET/paxdb](02_Postdoc_CONICET/paxdb)

- **Ankyrin modularity: PFAM and ELM enrichment**  
  Modern refactor of legacy scripts into a clear pipeline to study modular organisation of ankyrin repeats and associated motifs.  
  → `io.py`, `enrichment.py`, `conservation.py`, `elm.py`, `pipeline.py` + `README`.  
  👉 [CV/02_Postdoc_CONICET/ANKYRIN_MODULARITY](02_Postdoc_CONICET/ANKYRIN_MODULARITY)

- - **Figures and visuals**
  👉 [figures_visuals/](figures_visuals)

---

## CV

- **PTEFb_Regulation/PTEFb_PhD/ PhD work
  - `PTEFb_HEXIM_project/` – data, analysis scripts and figures for key mechanistic studies  
  - `PTEFb_Regulation/`, `PTEFb_PhD/` – thesis material, methods and results

- **02_Postdoc_CONICET/** – Structural bioinformatics / ankyrin modularity  
  - `paxdb/` – amino-acid usage evolution pipeline (see quick links above)  
  - `ANKYRIN_MODULARITY/` – ankyrin modularity and ELM/PFAM enrichment (see quick links above)

- **03_Industry_AIRBUS/** – critical systems development (Java)  
  - `README.md`, `critical_systems_development.md`, `soutenance_nina.pdf`

- **04_Postdoc_INSERM/** – CRCT / tumour ecosystem and network medicine  
  - `01/` – agent-based tumour ecosystem modelling + Python analysis pipeline  
  - `03_Deconvolution/` – bulk RNA-seq deconvolution pipeline  
  - `04_LungPredict/` – TF activity analysis in lung cancer  
  - `05_covid_network_medicine/` – multilayer network medicine  
  - `07_Barcodes/` – clonal dynamics from barcode counts

---
