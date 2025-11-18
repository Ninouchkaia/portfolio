# Postdoctoral Research – INSERM CRCT (2020–2023)
This section summarizes my work at the Cancer Research Center of Toulouse (CRCT),
focused on modeling tumor–immune interactions, transcriptomics, and integrative bioinformatics.

## **1. Modélisation multi-agents de l’écosystème tumoral (iScience 2023)**

<p align="center">
  <img src="figures/tumor_ecosystem_modeling.png" width="520px">
</p>

### **Objectif scientifique**  
Simuler l’évolution spatio-temporelle de tumeurs solides en incluant dynamique immunitaire, gradients diffusifs et interactions intercellulaires.

### **Contributions**
- Conception du modèle multi-agents (NetLogo).  
- Automatisation de centaines de simulations via OpenMOLE.  
- Analyse de sensibilité, exploration paramétrique, extraction de métriques.  
- Production de figures utilisées dans l’article iScience.

---

## **2. Polarisation macrophagique dans la leucémie lymphoïde chronique**

<p align="center">
  <img src="figures/macrophage_polarization.png" width="520px">
</p>

### **Objectif scientifique**  
Comprendre la transition des macrophages vers un état pro-tumoral (NLC) et identifier les programmes régulationnels sous-jacents.

### **Contributions**
- Estimation de l’activité de régulateurs transcriptionnels.  
- Participation au modèle dynamique de polarisation.  
- Analyse multi-datasets et validation croisée.  

---

## **3. Prédiction de la réponse à l’immunothérapie — GEMDECAN**

<p align="center">
  <img src="figures/immunotherapy_prediction.png" width="520px">
</p>

### **Objectif scientifique**  
Identifier des signatures transcriptomiques robustes associées à la réponse à l’immunothérapie (NK cells).

### **Contributions**
- Construction de parties de la pipeline bulk RNAseq.  
- Quantification, normalisation, tests différentiels (DESeq2).  
- Inférence TF activity (DoRothEA, VIPER).  
- Modélisation expression ↔ réponse thérapeutique.

---

## **4. Analyse du microenvironnement tumoral — LungPredict**

<p align="center">
  <img src="figures/lungpredict_deconvolution.png" width="520px">
</p>

### **Objectif scientifique**  
Caractériser la composition cellulaire et les programmes transcriptionnels des tumeurs pulmonaires afin d’identifier les signaux associés à l’évolution clinique.

### **Contributions**
- Prétraitement complet RNAseq : QC → trimming → alignement → quantification.  
- Déconvolution immunitaire (EPIC, CIBERSORTx).  
- Profilage régulationnel.  
- Contribution à l’analyse intégrée du microenvironnement.

---
# **5. Drug Repurposing for COVID-19 through Network Medicine**

[Accéder au dossier du projet](./05_covid_network_medicine/)

**Affiliation:** INSERM U1037 – Centre de Recherches en Cancérologie de Toulouse (CRCT)  
**Period:** 2020–2021  
**Publication:** [Network and Systems Medicine, 2020](https://www.liebertpub.com/doi/10.1089/nsm.2020.0011)


## 🧭 Context
At the onset of the COVID-19 pandemic, identifying potential therapeutic candidates required integrative strategies beyond single-target screening.  
This project used **network medicine** approaches to explore interactions between SARS-CoV-2 proteins, host cellular pathways, and drug targets, with the goal of repositioning existing compounds.  


## 🎯 Objectives
- Integrate multi-omics and molecular interaction data to construct a **virus–host–drug network**.  
- Identify biologically plausible drug candidates through **topological proximity** and **pathway enrichment**.  
- Test robustness of network-based predictions using simulated perturbations.  


## 🧪 Methods
- **Data integration:** Host–virus interactome from public datasets (BioGRID, IntAct), drug–target relationships from DrugBank and ChEMBL.  
- **Network modeling:** Weighted graph representation of molecular associations.  
- **Simulation:** Random rewiring and node removal to assess prediction stability.  
- **Analysis:** Centrality and community detection to highlight key druggable modules.  
- **Validation:** Cross-checking candidate lists with published clinical data and ongoing trials.  

## 💡 Contributions
- Implemented random network simulations to evaluate robustness of predicted drug–disease associations.  
- Automated analysis of node connectivity and topological metrics for ranking candidate drugs.  
- Contributed to visualization and reporting of systemic network perturbations.  
- Participated in manuscript review and interpretation of results.  

## 📘 Key Skills
Network biology · Data integration · Graph theory · Drug repurposing · Python/R scripting  

## 🔗 Reference
*Verstraete N.*, et al. *CovMulNet19, Integrating Proteins, Diseases, Drugs, and Symptoms: A Network Medicine Approach to COVID-19.*  
*Network and Systems Medicine*, 2020. [DOI:10.1089/nsm.2020.0011](https://www.liebertpub.com/doi/10.1089/nsm.2020.0011)

---

## **6. Effets systémiques du SARS-CoV-2**

<p align="center">
  <img src="figures/sarscov2_systemic.png" width="520px">
</p>

### **Objectif scientifique**  
Décrire comment les protéines virales perturbent les fonctions cellulaires dans différents tissus, et caractériser les effets systémiques.

### **Contributions**
- Analyses GO/Reactome/WikiPathways.  
- Identification de processus perturbés.  
- Contribution aux figures mécanistiques.

---

## **7. Dynamique clonale & single-cell RNAseq**

<p align="center">
  <img src="figures/scRNA_clonality.png" width="520px">
</p>

### **Objectif scientifique**  
Évaluer comment des clones tumoraux se diversifient sous traitement, et relier trajectoires transcriptomiques et résistance émergente.

### **Contributions**
- Pipeline scRNA-seq : filtrage, normalisation, clustering, UMAP.  
- Intégration barcodes → clones → programmes transcriptionnels.  
- Analyse de la diversité clonale et trajectoires d’états.  

---
