# **Portfolio - Research, Bioinformatics & Computational Biology**

## Overview

This GitHub repository documents the **methods**, **projects**, and **conceptual frameworks** I developed or contributed to across research institutions.  It serves as a structured overview of my approach to **data analysis**, **modeling**, and **scientific reproducibility**. 

Each folder corresponds to a specific research project and contains Markdown files summarizing the project structure, objectives, methods and eventual details to run the scripts.

---

## Research Focus

- **Bioinformatics & Data Analysis**
  - RNA-seq (bulk & single-cell), deconvolution, and gene regulatory inference  
  - Network medicine and multi-omics data integration  
  - Statistical modeling and visualization of biological data  

- **Modeling & Systems Biology**
  - Agent-based modeling of tumor-immune ecosystems  
  - Multi-scale simulations of cellular dynamics 
  - Network perturbation and system-level functional analysis  

- **Molecular & Structural Biology**
  - Protein evolution and structural modularity  
  - Regulatory mechanisms in transcription and post-transcriptional control  

- **Scientific Communication & Pedagogy**
  - Teaching, mentoring, and outreach in coding
  - Integration of artistic coding practices for public engagement  

---

## Research Projects

### **INSERM/CRCT - Postdoc (2020-2023)**

#### [1. Tumor Ecosystem Modeling](04_Postdoc_INSERM/01_AgentBasedModel/)

This project models the spatial and temporal dynamics of CLL tumor cells and nurse-like cells to understand how microenvironmental interactions shape patient-specific growth trajectories. The goal was to calibrate an agent-based model against experimental co-culture data using multi-objective optimization. The study identifies parameter regimes that reproduce observed CLL-NLC coexistence patterns.

#### [2. Macrophage Polarization in CLL](04_Postdoc_INSERM/02_BooleanModel)

This project investigates how transcription factors drive monocyte-to-macrophage polarization in the CLL microenvironment. A Boolean regulatory model was used to map distinct polarization trajectories and link them to experimental signatures. The analysis identifies regulatory determinants underlying M1/M2/NLC-like phenotypes.

#### [3. Predicting Response to Immunotherapy (GEMDECAN)](04_Postdoc_INSERM/03_RNASeqDeconvolution/)

The project explores how immune composition and gene expression programs predict checkpoint inhibitor response in solid tumors. A multi-cohort RNA-seq workflow quantifies immune infiltration and computes signature-based predictors of clinical outcome. The analysis aims to identify biomarkers associated with therapeutic benefit.

#### [4. LungPredict - Transcription Factor Network Analysis](04_Postdoc_INSERM/04_LungPredict)

This project examines how transcription factor activity patterns stratify lung cancer patients into biologically meaningful subgroups. Multi-omic annotations and PCA reveal TF-driven regulatory programs associated with clinical phenotypes. The analysis highlights regulatory axes shaping tumor heterogeneity.

#### [5. COVID-19 Drug Repurposing via Network Medicine](04_Postdoc_INSERM/05_NetworkMedicine/)

This work models how pharmacological perturbations propagate through SARS-CoV-2 host interaction layers. By comparing observed multilayer network structure to bootstrapped null models, it identifies drugs targeting essential viral-host modules. The objective is to highlight compounds with mechanistically grounded repurposing potential.

#### [6. Systemic Effects of SARS-CoV-2](04_Postdoc_INSERM/06_Sarscov2Interactome)

This project characterizes the systemic impact of SARS-CoV-2 by integrating viral-host interaction datasets and pathway annotations. Enrichment analysis highlights the cellular pathways and tissues most affected by viral proteins. The goal is to map functional vulnerabilities and host responses.

#### [7. Clonal Dynamics and scRNA-seq under Treatment](04_Postdoc_INSERM/07_BarcodesDrugScreening/)

This project investigates how clonal populations respond to targeted and chemotherapeutic drugs by quantifying barcode representation under treatment. Differential abundance analysis reveals drug-specific clonal signatures and functional clusters. Correlation networks map convergent and divergent therapeutic effects.

---

### **Software Engineering for Airbus Systems (Toulouse, 2017-2020)**
*Areas of interest: Development, integration and maintenance of PLM/PDM systems (APS, CASPARE, PASS V3 migration, Windchill v6→v11).*

---

### **CONICET/UBA - Postdoc (2013-2015)**
*Structural bioinformatics, protein evolution, and interaction modularity*

#### [1. Amino Acid Usage under Evolutionary Constraints](02_Postdoc_CONICET/paxdb/)

This work explores how metabolic cost and biosynthetic limitations influence amino acid composition across proteomes. By combining abundance-weighted PaxDB datasets with biochemical cost models, it tests evolutionary trade-offs between energy efficiency and proteome diversity. The results provide a quantitative framework linking metabolism to sequence composition.

#### [2. Structure and Dynamics of Ankyrin Repeats](02_Postdoc_CONICET/AnkyrinStructure/)

This project examines how structural variability and conserved features define the stability of ankyrin repeat proteins. A curated structural dataset enables analysis of contacts, repeat architecture, and positional variability. The study provides insight into how repeat proteins achieve robustness despite modular design.

#### [3. Functional Modularity of Ankyrin Proteins and Partners](02_Postdoc_CONICET/ANKYRIN_MODULARITY/)

The goal of this project is to understand how ankyrin domains combine with additional motifs to encode specific interaction functions. Conservation, Pfam/ELM enrichment and co-occurrence analyses reveal modular rules governing partner specificity. The work outlines how repeat and non-repeat elements assemble into functional architectures.

---

### **CNRS/IBENS - PhD (2008-2012)**
*Areas of interest: Regulation of transcriptional elongation and structure-function analysis*

#### [1. P-TEFb Regulation by HEXIM1 and HIV-1 Tat](01_PhD_CNRS)

This project identifies the Cyclin T1 structural determinants underlying binding to HEXIM1 and HIV-1 Tat. Mutagenesis and biochemical assays map critical hotspots involved in transcriptional elongation regulation. The findings clarify how viral proteins hijack the host machinery.

#### [2. P-TEFb mobility](01_PhD_CNRS)
This study examines how complex formation modulates P-TEFb nuclear mobility using FRAP/FLIP imaging. Truncation mutants and chemical perturbations reveal diffusion-binding equilibria within the nucleus. The work connects molecular assembly to dynamic regulation of transcription.

#### [3. P-TEFb conservation across Metazoa](01_PhD_CNRS)

This project investigates whether CyclinT-HEXIM interactions are conserved across metazoans. Comparative analysis and *C. elegans* validation show that critical interface residues remain functionally preserved. The study highlights ancestral mechanisms controlling transcriptional elongation.


### Research Projects Summary

| Project                                     | Summary                                                                 | Biological Question                                                    | Contribution                                                                                                 | Core Methods                                   |
|---------------------------------------------|-------------------------------------------------------------------------|------------------------------------------------------------------------|-----------------------------------------------------------------------------------------------------------------|------------------------------------------------|
| **Tumor Ecosystem Modeling**                | Calibration of a cancer-immune ABM using NSGA-II, knee-point selection, patient-specific fitting and statistical validation. | How do CLL-NLC interactions reproduce patient-specific tumor dynamics? | Full ABM pipeline (NetLogo → OpenMOLE → Python), NSGA-II calibration, patient fitting, validation, advanced stats | ABM, NSGA-II, RMSE scoring, PCA, sensitivity   |
| **Macrophage Polarization (Boolean Model)** | TF activity scoring (DoRothEA/VIPER), scaling and signature-based phenotype prediction. | Which TFs drive monocyte → macrophage polarization in CLL?             | TF activity scoring (VIPER/DoRothEA), scaling/normalisation                                                       | Regulatory inference, TF activity scoring      |
| **GEMDECAN Immunotherapy Response**         | Snakemake pipeline for RNA-seq TPM conversion, immune deconvolution and predictive modelling. | Can RNA-seq deconvolution predict responders?                          | Snakemake workflow: TPM, immune deconvolution, modelling                                                          | EPIC, quanTIseq, MCPCounter, signature scoring |
| **LungPredict (TF Networks)**               | TF activity inference, multi-omic annotation and PCA-based patient stratification. | How do TF programs stratify lung cancer patients?                      | TF activity inference, heatmaps, multi-omic annotation, PCA                                                       | VIPER, heatmap generation, clinical annotation |
| **COVID-19 Drug Repurposing**               | Multilayer network integration, bootstrap null models, drug z-score ranking. | Which drugs perturb SARS-CoV-2 host networks?                          | Multilayer networks, bootstrap null models, z-score ranking                                                       | Network integration, bootstrap, z-scores       |
| **SARS-CoV-2 Interactome**                  | Viral-host integration, Reactome enrichment, annotation harmonisation and figure generation. | Which pathways and tissues are most impacted?                          | Viral-host merge, Reactome enrichment, annotation, figures                                                        | Reactome, enrichment, annotation               |
| **Barcode Drug Screening**                  | QC, DESeq2 inputs, logFC matrices, drug-drug correlation networks and clustering. | How do clonal populations respond to drugs?                            | QC, DESeq2 input generation, logFC matrices, correlation networks                                                 | DESeq2, clustering, correlation networks       |
| **Amino Acid Usage (PaxDB)**                | Abundance-weighted amino acid frequencies, metabolic cost models and evolutionary trade-offs. | How do metabolic costs shape amino acid usage?                         | PaxDB parsing, weighted AA frequencies, cost modelling                                                            | Proteome parsing, evolutionary constraints     |
| **Ankyrin Structure**                       | Structural dataset construction, repeat detection, contact maps and variability profiles. | What features shape ankyrin repeat stability?                          | Dataset creation, repeat detection, structural analysis                                                           | PDB parsing, contact maps                      |
| **ANKYRIN_MODULARITY**                      | Pfam/ELM enrichment, domain co-occurrence statistics, conservation analysis and clustering. | How do ankyrin modules encode function?                                | Conservation, Pfam/ELM enrichment, co-occurrence networks                                                         | API retrieval, motif stats, clustering         |
| **P-TEFb Regulation**                       | Mutagenesis and structure-function mapping of CyclinT1. | Which residues control HEXIM1/Tat interactions?                        | Mutagenesis, interface mapping, assays                                                                            | Binding assays, structure-function             |
| **P-TEFb Mobility**                         | FRAP/FLIP diffusion modelling of P-TEFb complexes. | How does complex assembly affect mobility?                             | FRAP/FLIP quantification, diffusion modelling                                                                     | Imaging, biophysics                            |
| **P-TEFb Evolution**                        | Conservation analysis and experimental validation in *C. elegans*. | Is CyclinT-HEXIM interaction conserved?                                | Conservation analysis + *C. elegans* validation                                                                   | Evolutionary analysis, transgenics             |


---

## **Contact**

For collaborations or technical discussions, contact me at verstraete.nina[at]gmail.com.

---







































