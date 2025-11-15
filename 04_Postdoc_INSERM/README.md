# Postdoctoral Research – INSERM CRCT (2020–2023)
This section summarizes my work at the Cancer Research Center of Toulouse (CRCT),
focused on modeling tumor–immune interactions, transcriptomics, and integrative bioinformatics.
Travaux menés au Centre de Recherches en Cancérologie de Toulouse : modélisation multi-agents, analyses RNAseq, et bioinformatique translationnelle.

# **4. Postdoc Inserm — CRCT Toulouse (2020–2023)**

## **4.1. Modélisation multi-agents de l’écosystème tumoral (iScience 2023)**

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

### **Scripts associés**
- `scripts/ABM/run_batch.sh`  
- `scripts/ABM/aggregate_results.py`

---

## **4.2. Polarisation macrophagique dans la leucémie lymphoïde chronique**

<p align="center">
  <img src="figures/macrophage_polarization.png" width="520px">
</p>

### **Objectif scientifique**  
Comprendre la transition des macrophages vers un état pro-tumoral (NLC) et identifier les programmes régulationnels sous-jacents.

### **Contributions**
- Estimation de l’activité de régulateurs transcriptionnels.  
- Participation au modèle dynamique de polarisation.  
- Analyse multi-datasets et validation croisée.  

### **Scripts associés**
- `scripts/NLC/regulon_activity.R`

---

## **4.3. Prédiction de la réponse à l’immunothérapie — GEMDECAN**

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

### **Scripts associés**
- `scripts/bulk_RNAseq/Snakefile`  
- `scripts/bulk_RNAseq/dorothea_activity.R`

---

## **4.4. Analyse du microenvironnement tumoral — LungPredict**

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

### **Scripts associés**
- `scripts/LungPredict/deconvolution.R`

---

## **4.5. Repositionnement thérapeutique COVID-19 — Network Medicine**

<p align="center">
  <img src="figures/covid_network.png" width="520px">
</p>

### **Objectif scientifique**  
Identifier des candidats thérapeutiques via analyse d'interactions virus–hôte et signatures pharmacologiques, puis tester la robustesse par simulations sur réseaux.

### **Contributions**
- Randomisation réseau sur interactomes et graphes fonctionnels.  
- Analyse multi-omics intégrée.  
- Figures réseau et visualisations de modules impactés.

### **Scripts associés**
- `scripts/COVID/network_randomization.py`

---

## **4.6. Effets systémiques du SARS-CoV-2**

<p align="center">
  <img src="figures/sarscov2_systemic.png" width="520px">
</p>

### **Objectif scientifique**  
Décrire comment les protéines virales perturbent les fonctions cellulaires dans différents tissus, et caractériser les effets systémiques.

### **Contributions**
- Analyses GO/Reactome/WikiPathways.  
- Identification de processus perturbés.  
- Contribution aux figures mécanistiques.

### **Scripts associés**
- `scripts/SARSCoV2/enrichment_analysis.R`

---

## **4.7. Dynamique clonale & single-cell RNAseq**

<p align="center">
  <img src="figures/scRNA_clonality.png" width="520px">
</p>

### **Objectif scientifique**  
Évaluer comment des clones tumoraux se diversifient sous traitement, et relier trajectoires transcriptomiques et résistance émergente.

### **Contributions**
- Pipeline scRNA-seq : filtrage, normalisation, clustering, UMAP.  
- Intégration barcodes → clones → programmes transcriptionnels.  
- Analyse de la diversité clonale et trajectoires d’états.  

### **Scripts associés**
- `scripts/singlecell/preprocess.py`  
- `scripts/singlecell/clonal_integration.R`

---
