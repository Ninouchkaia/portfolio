# Nina Verstraete – Bioinformatics & Computational Biology

Computational biologist with experience in:
- End-to-end pipelines (Snakemake, Python, R)
- Bulk RNA-seq deconvolution and immune profiling
- High-throughput barcode data, QC and differential abundance
- Network medicine and multilayer graphs
- Agent-based modelling integrated with reproducible workflows

Experience in molecular and cell biology : here [CV/01_PhD_CNRS](CV/01_PhD_CNRS)

---




## **1. Agent-Based Modeling of the Tumor Ecosystem (iScience 2023)**

### **Reference**

*Verstraete N.*, Marku M., Domagala M., Arduin H., Bordenave J., Fournié J.-J., Ysebaert L., Poupot M., Pancaldi V.
*An agent-based model of monocyte differentiation into tumour-associated macrophages in chronic lymphocytic leukemia.*
*iScience*, 2023. [https://doi.org/10.1016/j.isci.2023.106897](https://doi.org/10.1016/j.isci.2023.106897)

### **Project folder**

[Access the folder 01_AgentBasedModel](CV/04_Postdoc_INSERM/01_AgentBasedModel)

---

## **2. Macrophage Polarization in Chronic Lymphocytic Leukemia**

### **Reference**

*Verstraete N.*, et al. *Regulatory Programs Driving Macrophage Polarization in Chronic Lymphocytic Leukemia.*
*(unpublished / internal project)*

### **Project Folder**

[Access the project folder 02_macrophage_polarization](CV/04_Postdoc_INSERM/02_macrophage_polarization/)

---

## **3. Prediction of Immunotherapy Response — GEMDECAN**

### **Reference**

Xie T., Pernet J., *Verstraete N.*, Madrid-Mencía M., Kuo M., Hucteau A., Coullomb A., Solórzano J., Delfour O., Cruzalegui F., Pancaldi V.
*GEM-DeCan: Improving tumor immune microenvironment profiling by the integration of novel gene expression and DNA methylation deconvolution signatures.*
*bioRxiv*, 2021. [https://doi.org/10.1101/2021.04.09.439207](https://doi.org/10.1101/2021.04.09.439207)

### **Project folder**

[Access the folder](CV/04_Postdoc_INSERM/03_Deconvolution)

---

## **4. Tumor Microenvironment Analysis — LungPredict**

### **Project folder**

[Access the folder](CV/04_Postdoc_INSERM/04_LungPredict)

---

## **5. Drug Repurposing for COVID-19 through Network Medicine**

### **Reference**

*Verstraete N.*, Jurman G., Bertagnolli G., Ghavasieh A., Pancaldi V., De Domenico M.
*CovMulNet19: Integrating proteins, diseases, drugs, and symptoms: A network medicine approach to COVID-19.*
*Network and Systems Medicine*, 2020.
[https://doi.org/10.1089/nsm.2020.0011](https://doi.org/10.1089/nsm.2020.0011)

### **Project folder**

[Access the folder](CV/04_Postdoc_INSERM/05_covid_network_medicine)

---

## **6. Systemic Effects of SARS-CoV-2**

### **Reference**

Ghavasieh A., Bontorin S., Artime O., *Verstraete N.*, De Domenico M.
*Multiscale statistical physics of the pan-viral interactome unravels the systemic nature of SARS-CoV-2 infections.*
*Communications Physics*, 2021.
[https://doi.org/10.1038/s42005-021-00582-8](https://doi.org/10.1038/s42005-021-00582-8)

### **Project folder**

[Access the folder](CV/04_Postdoc_INSERM/06_SystemicEffects)

---

## **7. Clonal Dynamics and Single-Cell RNA-seq**

### **Reference**

Chhouri H., *Verstraete N.*, Lee J., Alexandre D., Pancaldi V., Grumolato L.
*Inferring the mechanism of action of new drugs through the analysis of predetermined heterogeneous response in cancer cell subpopulations.*
NetBioMed 2022 (conference). *(No DOI)*

### **Project folder**

[Access the folder](CV/04_Postdoc_INSERM/07_Barcodes)

---















# **1. Agent-Based Modeling of the Tumor Ecosystem (iScience 2023)**

<p align="center"><img src="figures/tumor_ecosystem_modeling.png" width="520px"></p>

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

👉 [Access the folder 01_AgentBasedModel](CV/04_Postdoc_INSERM/01_AgentBasedModel)

---

# **2. Macrophage Polarization in CLL (Cancers 2020)**

<p align="center"><img src="figures/macrophage_polarization.png" width="520px"></p>

### **Scientific Objective**

Study monocyte → macrophage → NLC polarization through transcriptomics and Boolean network modeling.

### **Contributions**

* TF activity inference with DoRothEA / VIPER
* Boolean model analysis and perturbations
* Integration of multiple datasets

### **Reference**

*Verstraete N.*, et al. *Regulatory Programs Driving Macrophage Polarization in Chronic Lymphocytic Leukemia.*
*(unpublished / internal project)*

### **Project Folder**

[Access the project folder 02_macrophage_polarization](CV/04_Postdoc_INSERM/02_macrophage_polarization/)

---


# **3. Prediction of Immunotherapy Response — GEM-DeCan (bioRxiv 2021)**

<p align="center"><img src="figures/immunotherapy_prediction.png" width="520px"></p>

### **Scientific Objective**

Identify microenvironmental features predictive of immunotherapy response.

### **Contributions**

* Bulk RNA-seq preprocessing (QC → trimming → alignment → quantification)
* Deconvolution (DNA methylation + expression)
* Predictive modeling

### **Project folder**

[Access the folder](CV/04_Postdoc_INSERM/04_LungPredict)

---


# **4. Tumor Microenvironment Analysis — LungPredict**

<p align="center"><img src="figures/lungpredict_deconvolution.png" width="520px"></p>

### **Scientific Objective**

Analyze bulk RNA-seq from lung tumors to characterize immune composition and regulatory programs.

### **Contributions**

* QC, trimming, alignment, quantification
* Immune deconvolution
* TF activity inference and regulatory analysis

*(no publication yet)*

### **Project folder**

👉 [Access the folder](CV/04_Postdoc_INSERM/04_LungPredict)

---

# **5. Drug Repurposing for COVID-19 — Network Medicine (Network & Systems Medicine 2020)**

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

👉 [Access the folder](CV/04_Postdoc_INSERM/05_covid_network_medicine)

---

# **6. Systemic Effects of SARS-CoV-2 (Communications Physics 2021)**

<p align="center"><img src="figures/sarscov2_systemic.png" width="520px"></p>

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

👉 [Access the folder](CV/04_Postdoc_INSERM/06_SystemicEffects)

---

# **7. Clonal Dynamics & scRNA-seq under Treatment (NetBioMed 2022)**

<p align="center"><img src="figures/scRNA_clonality.png" width="520px"></p>

### **Scientific Objective**

Link barcode-based clonal expansion to transcriptomic reprogramming under drug treatment.

### **Contributions**

* Processing of barcode count matrices
* Differential clonal abundance
* scRNA-seq pipeline (QC, clustering, UMAP)
* Integration clonal identity ↔ gene expression

### **Reference**

Chhouri H., *Verstraete N.*, Lee J., Alexandre D., Pancaldi V., Grumolato L.
*Inferring the mechanism of action of new drugs through the analysis of predetermined heterogeneous response in cancer cell subpopulations.*
NetBioMed 2022 (conference). *(No DOI)*

### **Project folder**

👉 [Access the folder](CV/04_Postdoc_INSERM/07_Barcodes)













---
to remove !
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
