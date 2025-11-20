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

### **1. INSERM CRCT - Postdoctoral Research (2020-2023)**
*Areas of interest: Systems oncology, tumor-immune modeling, transcriptomics, and drug response prediction*

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

#### [6. Systemic Effects of SARS-CoV-2](04_Postdoc_INSERM/06_Sarscov2Interactome)
Data integration pipeline for viral-host interactions, enrichment (Reactome), network propagation, and ranking.

#### [7. Clonal Dynamics and scRNA-seq under Treatment](04_Postdoc_INSERM/07_BarcodesDrugScreening/)
Preprocessing, QC, DESeq2 inputs, fold-change networks, drug-drug correlation matrices, and figure generator for publication.

---

### **2. Software Engineering for Airbus Critical Systems (Toulouse, 2017-2020)**
*Areas of interest: Development, integration and maintenance of PLM/PDM systems (APS, CASPARE, PASS V3 migration, Windchill v6→v11).*

---

### **3. CONICET INQUIMAE - Postdoctoral Research (Buenos Aires, 2013-2015)**
*Structural bioinformatics, protein evolution, and interaction modularity*

#### [1. Amino Acid Usage under Evolutionary Constraints](02_Postdoc_CONICET/paxdb/)
Weighted frequencies, proteome cost metrics, correlations, PCA and domain-specific notebooks.

#### [2. Structure and Dynamics of Ankyrin Repeats](02_Postdoc_CONICET/AnkyrinStructure/)
Dataset construction, repeat detection, structural analysis, visualization.

#### [3. Functional Modularity of Ankyrin Proteins and Partners](02_Postdoc_CONICET/ANKYRIN_MODULARITY/)
Conservation, co-occurrence statistics, Pfam/ELM enrichment, protein family clustering and structural analysis.

---

### **4. CNRS IBENS - PhD (Paris, 2008-2012)**
*Areas of interest: Regulation of transcriptional elongation and structure-function analysis*

#### [1. P-TEFb Regulation by HEXIM1 and HIV-1 Tat](01_PhD_CNRS)
Work on P-TEFb, HEXIM1, Cyclin T1 mapping, structural hotspots and HIV Tat interference.

#### [2. P-TEFb mobility](01_PhD_CNRS)
Work on nuclear mobility of Cyclin T1 and P-TEFb complexes, using FRAP/FLIP to assess how truncation mutants and transcriptional inhibition modulate diffusion and partner association.

#### [3. P-TEFb conservation across Metazoa](01_PhD_CNRS)
Work on evolutionary conservation of the HEXIM-Cyclin T complex, identifying C. elegans homologs and showing that key interaction residues are functionally preserved.

---

## **Research Projects Summary**

### Research Projects Summary

| Project                                     | Summary                                                                 | Biological Question                                                    | Contribution                                                                                                 | Core Methods                                   |
|---------------------------------------------|-------------------------------------------------------------------------|------------------------------------------------------------------------|-----------------------------------------------------------------------------------------------------------------|------------------------------------------------|
| **Tumor Ecosystem Modeling**                | Calibration of a cancer–immune ABM using NSGA-II, knee-point selection, patient-specific fitting and statistical validation. | How do CLL–NLC interactions reproduce patient-specific tumor dynamics? | Full ABM pipeline (NetLogo → OpenMOLE → Python), NSGA-II calibration, patient fitting, validation, advanced stats | ABM, NSGA-II, RMSE scoring, PCA, sensitivity   |
| **Macrophage Polarization (Boolean Model)** | TF activity scoring (DoRothEA/VIPER), scaling and signature-based phenotype prediction. | Which TFs drive monocyte → macrophage polarization in CLL?             | TF activity scoring (VIPER/DoRothEA), scaling/normalisation                                                       | Regulatory inference, TF activity scoring      |
| **GEMDECAN Immunotherapy Response**         | Snakemake pipeline for RNA-seq TPM conversion, immune deconvolution and predictive modelling. | Can RNA-seq deconvolution predict responders?                          | Snakemake workflow: TPM, immune deconvolution, modelling                                                          | EPIC, quanTIseq, MCPCounter, signature scoring |
| **LungPredict (TF Networks)**               | TF activity inference, multi-omic annotation and PCA-based patient stratification. | How do TF programs stratify lung cancer patients?                      | TF activity inference, heatmaps, multi-omic annotation, PCA                                                       | VIPER, heatmap generation, clinical annotation |
| **COVID-19 Drug Repurposing**               | Multilayer network integration, bootstrap null models, drug z-score ranking. | Which drugs perturb SARS-CoV-2 host networks?                          | Multilayer networks, bootstrap null models, z-score ranking                                                       | Network integration, bootstrap, z-scores       |
| **SARS-CoV-2 Interactome**                  | Viral–host integration, Reactome enrichment, annotation harmonisation and figure generation. | Which pathways and tissues are most impacted?                          | Viral-host merge, Reactome enrichment, annotation, figures                                                        | Reactome, enrichment, annotation               |
| **Barcode Drug Screening**                  | QC, DESeq2 inputs, logFC matrices, drug–drug correlation networks and clustering. | How do clonal populations respond to drugs?                            | QC, DESeq2 input generation, logFC matrices, correlation networks                                                 | DESeq2, clustering, correlation networks       |
| **Amino Acid Usage (PaxDB)**                | Abundance-weighted amino acid frequencies, metabolic cost models and evolutionary trade-offs. | How do metabolic costs shape amino acid usage?                         | PaxDB parsing, weighted AA frequencies, cost modelling                                                            | Proteome parsing, evolutionary constraints     |
| **Ankyrin Structure**                       | Structural dataset construction, repeat detection, contact maps and variability profiles. | What features shape ankyrin repeat stability?                          | Dataset creation, repeat detection, structural analysis                                                           | PDB parsing, contact maps                      |
| **ANKYRIN_MODULARITY**                      | Pfam/ELM enrichment, domain co-occurrence statistics, conservation analysis and clustering. | How do ankyrin modules encode function?                                | Conservation, Pfam/ELM enrichment, co-occurrence networks                                                         | API retrieval, motif stats, clustering         |
| **P-TEFb Regulation**                       | Mutagenesis and structure–function mapping of CyclinT1. | Which residues control HEXIM1/Tat interactions?                        | Mutagenesis, interface mapping, assays                                                                            | Binding assays, structure-function             |
| **P-TEFb Mobility**                         | FRAP/FLIP diffusion modelling of P-TEFb complexes. | How does complex assembly affect mobility?                             | FRAP/FLIP quantification, diffusion modelling                                                                     | Imaging, biophysics                            |
| **P-TEFb Evolution**                        | Conservation analysis and experimental validation in *C. elegans*. | Is CyclinT–HEXIM interaction conserved?                                | Conservation analysis + *C. elegans* validation                                                                   | Evolutionary analysis, transgenics             |


---

## **Contact**

For collaborations or technical discussions, contact me at verstraete.nina[at]gmail.com.

---






































