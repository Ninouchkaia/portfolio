# **Portfolio - Research, Bioinformatics & Computational Biology**

This repository gathers a selection of scientific, technical and software-engineering projects developed during my PhD, postdocs (CNRS, CONICET, INSERM), and industry experience (Airbus/Capgemini). It showcases **biological modeling**, **multi-omics analysis**, **network science**, **data collection and analysis pipelines** and **reproducible workflows** in Python, R, Snakemake, NetLogo and OpenMOLE.

<div style="margin-top:10px"></div>

---

# **Structure Overview**

```
.
├── 01_PhD_CNRS/            ← Functional genomics & transcription regulation
├── 02_Postdoc_CONICET/     ← Evolution, structural bioinformatics, PaxDB, ankyrin modularity
├── 03_Industry_AIRBUS/     ← Software engineering in critical systems
├── 04_Postdoc_INSERM/      ← Large-scale bioinformatics pipelines (7 major projects)
└── figures_visuals/        ← Visual summaries for portfolio and presentations
```

---
# **2020-2023 Tumor heterogeneity in immuno-oncology (INSERM)**
[Postdoc, CRCT, Toulouse](04_Postdoc_INSERM)

### **1. Agent-Based Modeling of Tumor Ecosystems (INSERM)**

`04_Postdoc_INSERM/01_AgentBasedModel/`
**NetLogo → OpenMOLE (NSGA-II) → Python validation**
A complete multi-objective calibration pipeline for a cancer–immune ABM, including BehaviorSpace generation, Pareto/knee-point selection, patient-specific fitting and advanced statistical analysis.

### **2. Boolean Modeling of Immunotherapy Responses**

`04_Postdoc_INSERM/02_BooleanModel/`
Pipeline for computing transcription factor signatures (DoRothEA-like), scaling, normalisation and prediction modules.

### **3. RNA-Seq Deconvolution & Immunotherapy Prediction**

`04_Postdoc_INSERM/03_RNASeqDeconvolution/`
Full Snakemake-style workflow in R/Python for TPM conversion, signature scoring, deconvolution (EPIC, MCPCounter, quanTIseq), and predictive modelling.

### **4. LungPredict – Transcription Factor Network Analysis**

`04_Postdoc_INSERM/04_LungPredict/`
Large R pipeline for TF activity inference, heatmaps, multi-omics integrative annotation and patient stratification.

### **5. Network Medicine – Multilayer Graph Analysis**

`04_Postdoc_INSERM/05_NetworkMedicine/`
Python pipeline for multilayer protein–GO networks, bootstrap null models, z-scores and ranking analyses.

### **6. SARS-CoV-2 Interactome & Systemic Effects**

`04_Postdoc_INSERM/06_Sarscov2Interactome/`
Data integration pipeline for viral–host interactions, enrichment (Reactome), network propagation, and ranking.

### **7. Drug Screening with DNA Barcodes**

`04_Postdoc_INSERM/07_BarcodesDrugScreening/`
Preprocessing, QC, DESeq2 inputs, fold-change networks, drug–drug correlation matrices, and figure generator for publication.

---

# **2017-2020 Industry – Software Engineering for Critical Systems (Airbus)**
(in construction)

`03_Industry_AIRBUS/`
Development, integration and maintenance of PLM/PDM systems (APS, CASPARE, PASS V3 migration, Windchill v6→v11).
(in construction)

---

# **2013-2015 Evolution & Structural Bioinformatics (CONICET)**
[Postdoc, UBA, Buenos Aires](02_Postdoc_CONICET)

### **Ankyrin Repeat Modularity Pipeline**

`02_Postdoc_CONICET/ANKYRIN_MODULARITY/`
Conservation, co-occurrence statistics, Pfam/ELM enrichment, protein family clustering and structural analysis.

### **Amino Acid Usage & Metabolic Cost (PaxDB)**

`02_Postdoc_CONICET/paxdb/`
Weighted frequencies, proteome cost metrics, correlations, PCA and domain-specific notebooks.

---

# **2008-2012 PhD – Regulation of Transcription Elongation (CNRS)**
[PhD, ENS, Paris](01_PhD_CNRS)

`01_PhD_CNRS/`
Work on P-TEFb, HEXIM1, Cyclin T1 mapping, structural hotspots and HIV Tat interference.

---

# **Visual Summaries**

`figures_visuals/`
Contains PNG visual overviews of the main projects (ABM, Barcodes, Network Medicine, GEMDECAN, LungPredict, Systemic COVID).

---

# **How to Navigate**

Every sub-project contains its own **README**, with:

* scientific context
* pipeline description
* file tree
* commands
* figures and results
* reproducibility instructions

The top-level repository serves as a **CV companion**, connecting all research and engineering work.

---

# **Contact**

For collaborations or technical discussions, contact me at verstraete.nina[at]gmail.com.

---






