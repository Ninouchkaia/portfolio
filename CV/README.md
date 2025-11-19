# **Portfolio - Research, Bioinformatics & Computational Biology**

## Overview

This GitHub repository documents the **methods**, **projects**, and **conceptual frameworks** I developed or contributed to across academia and research institutions.  It serves as a structured overview of my approach to **data analysis**, **modeling**, and **scientific reproducibility**. Each folder corresponds to a specific research project and contains Markdown files summarizing the objectives, methods, and results of individual projects.

---

## Research Focus

- **Bioinformatics & Data Analysis**
  - RNA-seq (bulk & single-cell), deconvolution, and gene regulatory inference  
  - Network medicine and multi-omics data integration  
  - Statistical modeling and visualization of biological data  

- **Modeling & Systems Biology**
  - Agent-based modeling of tumor-immune ecosystems  
  - Multi-scale simulations of cellular dynamics (NetLogo, OpenMOLE)  
  - Network perturbation and system-level functional analysis  

- **Molecular & Structural Biology**
  - Protein evolution and structural modularity  
  - Regulatory mechanisms in transcription and post-transcriptional control  

- **Scientific Communication & Pedagogy**
  - Teaching, mentoring, and outreach in bioinformatics and digital creativity  
  - Integration of artistic coding practices (music & generative systems) for public engagement  

---

## Research Projects

### **INSERM CRCT - Postdoctoral Research (2020-2023)**
*Systems oncology, tumor-immune modeling, transcriptomics, and drug response prediction*

#### [1. Tumor Ecosystem Modeling](04_Postdoc_INSERM/01_AgentBasedModel/)
A complete multi-objective calibration pipeline for a cancer-immune ABM, including BehaviorSpace generation, Pareto/knee-point selection, patient-specific fitting and advanced statistical analysis.

#### [2. Macrophage Polarization in CLL](04_Postdoc_INSERM/02_BooleanModel)
Pipeline for computing transcription factor signatures (DoRothEA-like), scaling, normalisation and prediction modules.

#### [3. Predicting Response to Immunotherapy (GEMDECAN)](04_Postdoc_INSERM/03_RNASeqDeconvolution/)
Snakemake workflow in R/Python for TPM conversion, signature scoring, deconvolution (EPIC, MCPCounter, quanTIseq), and predictive modelling.

#### [4. LungPredict - Transcription Factor Network Analysis](04_Postdoc_INSERM/04_LungPredict)
R pipeline for TF activity inference, heatmaps, multi-omics integrative annotation and patient stratification.

#### [5. COVID-19 Drug Repurposing via Network Medicine](04_Postdoc_INSERM/05_NetworkMedicine/)
Python pipeline for multilayer protein-GO networks, bootstrap null models, z-scores and ranking analyses.

#### [6. Systemic Effects of SARS-CoV-2](04_Postdoc_INSERM/06_Sarscov2Interactome/)
Data integration pipeline for viral-host interactions, enrichment (Reactome), network propagation, and ranking.

#### [7. Clonal Dynamics and scRNA-seq under Treatment](04_Postdoc_INSERM/07_BarcodesDrugScreening/)
Preprocessing, QC, DESeq2 inputs, fold-change networks, drug-drug correlation matrices, and figure generator for publication.

---

### **Software Engineering for Airbus Critical Systems (Toulouse, 2017-2020)**
*Development, integration and maintenance of PLM/PDM systems (APS, CASPARE, PASS V3 migration, Windchill v6→v11).*

---

### **CONICET INQUIMAE - Postdoctoral Research (Buenos Aires, 2013-2015)**
*Structural bioinformatics, protein evolution, and interaction modularity*

#### [Amino Acid Usage under Evolutionary Constraints](02_Postdoc_CONICET/paxdb/)
Weighted frequencies, proteome cost metrics, correlations, PCA and domain-specific notebooks.

#### [Structure and Dynamics of Ankyrin Repeats](02_Postdoc_CONICET/AnkyrinStructure/)
Dataset construction, repeat detection, structural analysis, visualization.

#### [Functional Modularity of Ankyrin Proteins and Partners](02_Postdoc_CONICET/ANKYRIN_MODULARITY/)
Conservation, co-occurrence statistics, Pfam/ELM enrichment, protein family clustering and structural analysis.

---

### **CNRS IBENS - PhD (Paris, 2008-2012)**
*Regulation of transcriptional elongation and structure-function analysis*

#### [P-TEFb Regulation by HEXIM1 and HIV-1 Tat](01_PhD_CNRS)
Work on P-TEFb, HEXIM1, Cyclin T1 mapping, structural hotspots and HIV Tat interference.

---

## **Visual Summaries**
[Visuals](figures_visuals/)
Contains PNG visual overviews of the main projects (ABM, Barcodes, Network Medicine, GEMDECAN, LungPredict, Systemic COVID).

---

## **How to Navigate**

Every sub-project contains its own **README**, with:

* scientific context
* pipeline description
* commands
* figures and results

---

## **Contact**

For collaborations or technical discussions, contact me at verstraete.nina[at]gmail.com.

---























