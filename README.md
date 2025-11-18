# Nina Verstraete – Bioinformatics & Computational Biology

Computational biologist with experience in:
- End-to-end pipelines (Snakemake, Python, R)
- Bulk RNA-seq deconvolution and immune profiling
- High-throughput barcode data, QC and differential abundance
- Network medicine and multilayer graphs
- Agent-based modelling integrated with reproducible workflows

For experience in molecular and cell biology, see here: [CV/01_PhD_CNRS](CV/01_PhD_CNRS)

---

## **1. Agent-Based Modeling of the Tumor Ecosystem (iScience 2023)**

<p align="center"><img src="CV/figures_visuals/abm.png" width="520px"></p>

### **Scientific Objective**

Simulate the spatio-temporal evolution of CLL cells and monocyte-derived myeloid cells, including NLC differentiation and microenvironmental support.

### **Contributions**

* NetLogo ABM design
* OpenMOLE NSGA-II exploration
* Sensitivity analysis, patient-specific calibration
* Python validation and figure-generation pipeline

### **Reference**

*Verstraete N.*, Marku M., Domagala M., Arduin H., Bordenave J., Fournié J.-J., Ysebaert L., Poupot M., Pancaldi V.
*An agent-based model of monocyte differentiation into tumour-associated macrophages in chronic lymphocytic leukemia.*
*iScience*, 2023. [https://doi.org/10.1016/j.isci.2023.106897](https://doi.org/10.1016/j.isci.2023.106897)

### **Project folder**

[Access the folder 01_AgentBasedModel](../CV/04_Postdoc_INSERM/01_AgentBasedModel)

---

## **2. Macrophage Polarization in CLL (Cancers 2020)**

### **Scientific Objective**

Study monocyte → macrophage → NLC polarization through transcriptomics and Boolean network modeling.

### **Contributions**

* TF activity inference with DoRothEA / VIPER
* Boolean model analysis and perturbations
* Integration of multiple datasets

### **Reference**

Marku, M., *Verstraete, N.*, Raynal, F., Madrid-Mencía, M., Domagala, M., Fournié, J.-J., Ysebaert, L., Poupot, M., & Pancaldi, V. 
*Insights on TAM formation from a Boolean model of macrophage polarization based on in vitro studies.* 
*Cancers*, 2020. [https://doi.org/10.3390/cancers12123664](https://doi.org/10.3390/cancers12123664)

### **Project Folder**

[Access the project folder 02_macrophage_polarization](CV/04_Postdoc_INSERM/02_macrophage_polarization/)

---


## **3. GEM-DeCan: Deconvolution Pipeline and Prediction of Immunotherapy Response **

<p align="center"><img src="CV/figures_visuals/gemdecan.png" width="320px"></p>

### **Scientific Objective**

Identify microenvironmental features predictive of immunotherapy response.

### **Contributions**

* Bulk RNA-seq preprocessing (QC → trimming → alignment → quantification)
* Deconvolution (DNA methylation + expression)
* Predictive modeling

### **Project folder**

[Access the folder 03_Deconvolution](CV/04_Postdoc_INSERM/03_Deconvolution)

---


## **4. Tumor Microenvironment Analysis — LungPredict**

<p align="center"><img src="CV/figures_visuals/LP.png" width="520px"></p>

### **Scientific Objective**

Analyze bulk RNA-seq from lung tumors to characterize immune composition and regulatory programs.

### **Contributions**

* QC, trimming, alignment, quantification
* Immune deconvolution
* TF activity inference and regulatory analysis


### **Project folder**

[Access the folder 04_LungPredict](CV/04_Postdoc_INSERM/04_LungPredict)

---

## **5. Drug Repurposing for COVID-19 — Network Medicine (Network & Systems Medicine 2020)**
<p align="center"><img src="CV/figures_visuals/multilayer.png" width="520px"></p>

### **Scientific Objective**

Integrate proteins, drugs, diseases, and symptoms in a multilayer graph to identify repurposing candidates.

### **Contributions**

* Random network simulations
* Drug ranking via proximity metrics
* Robustness analysis (bootstrap)
* Visualisation of multilayer interactions

### **Reference**

*Verstraete N.*, Jurman G., Bertagnolli G., Ghavasieh A., Pancaldi V., De Domenico M.
*CovMulNet19: Integrating proteins, diseases, drugs, and symptoms: A network medicine approach to COVID-19.*
*Network and Systems Medicine*, 2020.
[https://doi.org/10.1089/nsm.2020.0011](https://doi.org/10.1089/nsm.2020.0011)

### **Project folder**

[Access the folder 05_covid_network_medicine](CV/04_Postdoc_INSERM/05_covid_network_medicine)

---

## **6. Systemic Effects of SARS-CoV-2 (Communications Physics 2021)**

<p align="center"><img src="CV/figures_visuals/systemic.png" width="520px"></p>

### **Scientific Objective**

Describe how viral proteins perturb host cell systems at multiple scales.

### **Contributions**

* GO/Reactome/WikiPathways enrichment
* Viral module clustering
* Interpretation of systemic perturbations

### **Reference**

Ghavasieh A., Bontorin S., Artime O., *Verstraete N.*, De Domenico M.
*Multiscale statistical physics of the pan-viral interactome unravels the systemic nature of SARS-CoV-2 infections.*
*Communications Physics*, 2021.
[https://doi.org/10.1038/s42005-021-00582-8](https://doi.org/10.1038/s42005-021-00582-8)

### **Project folder**

[Access the folder 06_SystemicEffects](CV/04_Postdoc_INSERM/06_SystemicEffects)

---

## **7. Clonal Dynamics & scRNA-seq under Treatment**

<p align="center"><img src="CV/figures_visuals/barcodes.png" width="520px"></p>

### **Scientific Objective**

Link barcode-based clonal expansion to transcriptomic reprogramming under drug treatment.

### **Contributions**

* Processing of barcode count matrices
* Differential clonal abundance
* scRNA-seq pipeline (QC, clustering, UMAP)
* Integration clonal identity ↔ gene expression

### **Project folder**

[Access the folder 07_Barcodes](CV/04_Postdoc_INSERM/07_Barcodes)
