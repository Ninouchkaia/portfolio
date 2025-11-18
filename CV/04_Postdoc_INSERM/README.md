# Postdoctoral Research – INSERM CRCT (2020–2023)
This section summarizes my work at the Cancer Research Center of Toulouse (CRCT - Team NetBio²), focused on modeling tumor–immune interactions, transcriptomics, and integrative bioinformatics.


---
## **1. Agent-Based Modeling of the Tumor Ecosystem (iScience 2023)**
<p align="center">
  <img src="figures/tumor_ecosystem_modeling.png" width="520px">
</p>

### **Scientific Objective**
Simulate the spatio-temporal evolution of solid tumors by integrating immune dynamics, diffusive gradients, and intercellular interactions.

### **Contributions**
* Design of the agent-based model (NetLogo).
* Automation of hundreds of simulations using OpenMOLE.
* Sensitivity analysis, parameter exploration, metric extraction.
* Generation of figures used in the iScience publication.

  
---
## **2. Macrophage Polarization in Chronic Lymphocytic Leukemia**
<p align="center">
  <img src="figures/macrophage_polarization.png" width="520px">
</p>

### **Scientific Objective**
Understand the transition of macrophages toward a pro-tumoral state (NLC) and identify the underlying regulatory programs.

### **Contributions**
* Estimation of transcription factor activity.
* Contribution to the dynamic model of macrophage polarization.
* Multi-dataset analysis and cross-validation.


---
## **3. Prediction of Immunotherapy Response — GEMDECAN**
<p align="center">
  <img src="figures/immunotherapy_prediction.png" width="520px">
</p>

### **Scientific Objective**
Identify robust transcriptomic signatures associated with immunotherapy response (NK cells).

### **Contributions**
* Construction of parts of the bulk RNA-seq pipeline.
* Quantification, normalization, differential expression (DESeq2).
* Inference of TF activity (DoRothEA, VIPER).
* Modeling of expression–response associations.


---
## **4. Tumor Microenvironment Analysis — LungPredict**
<p align="center">
  <img src="figures/lungpredict_deconvolution.png" width="520px">
</p>

### **Scientific Objective**
Characterize cell-type composition and transcriptional programs in lung tumors to identify signals associated with clinical outcomes.

### **Contributions**
* Full RNA-seq preprocessing: QC → trimming → alignment → quantification.
* Immune deconvolution (EPIC, CIBERSORTx).
* Regulatory profiling.
* Contribution to the integrated microenvironment analysis.


---
## **5. Drug Repurposing for COVID-19 through Network Medicine**
[Access the project folder](./05_covid_network_medicine/)

### **Scientific Objective**
* Integrate multi-omics and molecular interaction data to construct a virus–host–drug network.
* Identify biologically plausible drug candidates based on topological proximity and pathway enrichment.
* Assess robustness of network-based predictions using simulated perturbations.

### **Scientific Objective**
* Implemented random network simulations to evaluate robustness of predicted drug–disease associations.
* Automated analysis of node connectivity and topological metrics for drug ranking.
* Contributed to visualization and reporting of systemic network perturbations.
* Participated in manuscript review and interpretation of results.

## **Reference**

*Verstraete N.*, et al. *CovMulNet19, Integrating Proteins, Diseases, Drugs, and Symptoms: A Network Medicine Approach to COVID-19.*
*Network and Systems Medicine*, 2020. DOI:10.1089/nsm.2020.0011


---
## **6. Systemic Effects of SARS-CoV-2**
<p align="center">
  <img src="figures/sarscov2_systemic.png" width="520px">
</p>

### **Scientific Objective**
Describe how viral proteins perturb cellular functions across multiple tissues and characterize systemic effects.

### **Scientific Objective**
* GO/Reactome/WikiPathways analyses.
* Identification of perturbed processes.
* Contribution to mechanistic figures.


---
## **7. Clonal Dynamics and Single-Cell RNA-seq**
<p align="center">
  <img src="figures/scRNA_clonality.png" width="520px">
</p>

### **Scientific Objective**
Evaluate how tumor clones diversify under treatment and connect transcriptomic trajectories to emerging resistance.

### **Contributions**
* scRNA-seq pipeline: filtering, normalization, clustering, UMAP.
* Integration of barcodes → clones → transcriptional programs.
* Analysis of clonal diversity and state trajectories.

---
