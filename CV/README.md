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
### 🎯 **INSERM CRCT – Postdoctoral Research (2020–2023)**
*Systems oncology, tumor–immune modeling, transcriptomics, and drug response prediction*

- [Tumor Ecosystem Modeling](postdoc_CRCT_2020-2023/tumor_ecosystem_modeling.md): A complete multi-objective calibration pipeline for a cancer–immune ABM, including BehaviorSpace generation, Pareto/knee-point selection, patient-specific fitting and advanced statistical analysis.  
- [Macrophage Polarization in CLL](postdoc_CRCT_2020-2023/macrophage_polarization.md): Pipeline for computing transcription factor signatures (DoRothEA-like), scaling, normalisation and prediction modules.  
- [Predicting Response to Immunotherapy (GEMDECAN)](postdoc_CRCT_2020-2023/immunotherapy_prediction.md)  
- [Tumor Microenvironment Analysis (LungPredict)](postdoc_CRCT_2020-2023/tumor_microenvironment_LungPredict.md)  
- [COVID-19 Drug Repurposing via Network Medicine](postdoc_CRCT_2020-2023/covid_network_medicine.md)  
- [Systemic Effects of SARS-CoV-2](postdoc_CRCT_2020-2023/sarscov2_systemic_effects.md)  
- [Clonal Dynamics and scRNA-seq under Treatment](postdoc_CRCT_2020-2023/clonal_dynamics_scRNAseq.md)

---

### 🧩 **INQUIMAE – CONICET (Buenos Aires, 2013–2015)**
*Structural bioinformatics, protein evolution, and interaction modularity*

- [Amino Acid Usage under Evolutionary Constraints](postdoc_CONICET_2013-2015/aa_usage_evolution.md)  
- [Structure and Dynamics of Ankyrin Repeats](postdoc_CONICET_2013-2015/ankyrin_structure_dynamics.md)  
- [Functional Modularity of Ankyrin Proteins and Partners](postdoc_CONICET_2013-2015/ankyrin_modularity.md)

---

### 🧫 **PhD – ENS / CNRS IBENS (Paris, 2008–2012)**
*Regulation of transcriptional elongation and structure–function analysis*

- [P-TEFb Regulation by HEXIM1 and HIV-1 Tat](phd_ENS_2008-2012/transcription_regulation_hexim_tat.md)

---

### **2020-2023 Tumor heterogeneity in immuno-oncology (INSERM)**

#### [**Tumor Ecosystem Modeling**](04_Postdoc_INSERM/01_AgentBasedModel/)
A complete multi-objective calibration pipeline for a cancer–immune ABM, including BehaviorSpace generation, Pareto/knee-point selection, patient-specific fitting and advanced statistical analysis.

#### [**2. Macrophage Polarization in CLL**](04_Postdoc_INSERM/02_BooleanModel)
Pipeline for computing transcription factor signatures (DoRothEA-like), scaling, normalisation and prediction modules.

#### [**3. Predicting Response to Immunotherapy (GEMDECAN)**](04_Postdoc_INSERM/03_RNASeqDeconvolution/)
Snakemake workflow in R/Python for TPM conversion, signature scoring, deconvolution (EPIC, MCPCounter, quanTIseq), and predictive modelling.

#### [**4. LungPredict - Transcription Factor Network Analysis**](04_Postdoc_INSERM/04_LungPredict)
R pipeline for TF activity inference, heatmaps, multi-omics integrative annotation and patient stratification.

#### [**5. COVID-19 Drug Repurposing via Network Medicine**](04_Postdoc_INSERM/05_NetworkMedicine/)
Python pipeline for multilayer protein–GO networks, bootstrap null models, z-scores and ranking analyses.

#### (**6. Systemic Effects of SARS-CoV-2**)(04_Postdoc_INSERM/06_Sarscov2Interactome/)
Data integration pipeline for viral–host interactions, enrichment (Reactome), network propagation, and ranking.

#### **[7. Clonal Dynamics and scRNA-seq under Treatment**](04_Postdoc_INSERM/07_BarcodesDrugScreening/)
Preprocessing, QC, DESeq2 inputs, fold-change networks, drug–drug correlation matrices, and figure generator for publication.

---

# **2017-2020 Industry – Software Engineering for Critical Systems (Airbus)**
(in construction)

Development, integration and maintenance of PLM/PDM systems (APS, CASPARE, PASS V3 migration, Windchill v6→v11).
(in construction)

---

# **2013-2015 Evolution & Structural Bioinformatics (CONICET)**
[Postdoc, UBA, Buenos Aires](02_Postdoc_CONICET)

### **1. Ankyrin Repeat Modularity Pipeline**
[1. AnkyrinModularity](02_Postdoc_CONICET/ANKYRIN_MODULARITY/)

Conservation, co-occurrence statistics, Pfam/ELM enrichment, protein family clustering and structural analysis.

### **2. Amino Acid Usage & Metabolic Cost (PaxDB)**
[2. AminoAcidsAbundance](02_Postdoc_CONICET/paxdb/)

Weighted frequencies, proteome cost metrics, correlations, PCA and domain-specific notebooks.

---

# **2008-2012 PhD – Regulation of Transcription Elongation (CNRS)**
[PhD, ENS, Paris](01_PhD_CNRS)

Work on P-TEFb, HEXIM1, Cyclin T1 mapping, structural hotspots and HIV Tat interference.

---

# **Visual Summaries**
[Visuals](figures_visuals/)

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















