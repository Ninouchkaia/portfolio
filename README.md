# **Portfolio scientifique — Nina Verstraete**

Ce dépôt rassemble une sélection détaillée de mes travaux en bioinformatique, modélisation, analyse de données multi-omics, biologie moléculaire et ingénierie logicielle.  
Il reflète un parcours au croisement de la biologie, de la computation et des systèmes complexes, couvrant près de quinze années d’expériences en recherche académique, R&D et analyse de données.

L’objectif est de présenter clairement les approches méthodologiques que j’ai développées, les raisons de ces choix, et les contributions scientifiques associées — dans une perspective reproductible et structurée.

---

# **📌 Structure générale du portfolio**

- **Thèse — CNRS / ENS Paris (2008–2012)**  
  *(Régulation transcriptionnelle, microscopie avancée, mutagenèse, biologie cellulaire)*

- **Postdoc — CONICET Buenos Aires (2013–2015)**  
  *(Bioinformatique structurale, évolution moléculaire, modularité protéique)*

- **Ingénierie logicielle — Airbus Defence & Space (2017–2019)**  
  *(Développement logiciel, systèmes critiques, automatisation, architecture)*

- **Postdoc — CRCT Inserm Toulouse (2020–2023)**  
  *(Modélisation multi-agents, RNAseq, single-cell, barcoding, immunologie computationnelle)*

Chaque projet inclut :
- **Objectif scientifique**  
- **Contributions principales**  
- **Figures illustratives**  
- **Chemins génériques vers scripts/notebooks**

---

# **1. Thèse — CNRS / ENS Paris (2008–2012)**  
## **Régulation transcriptionnelle de P-TEFb par HEXIM1 et TAT (HIV-1)**

<p align="center">
  <img src="figures/transcription_hexintat.png" width="520px">
</p>

### **Objectif scientifique**  
Comprendre comment P-TEFb (CDK9/Cyclin T1) module la transition vers l’élongation productive, et comment il est régulé ou détourné par les complexes 7SK/HEXIM1 et la protéine virale HIV-1 TAT.

### **Contributions**
- Construction de **bibliothèques mutantes** (HEXIM1, Cyclin T1, TAT) par mutagenèse dirigée.  
- Criblage **double-hybride inverse** pour cartographier les surfaces d’interaction.  
- Génération de **lignes cellulaires stables** et transgenèse dans *C. elegans*.  
- Biologie moléculaire et cellulaire : co-IP, western/southern, extraction ADN/ARN, imagerie.  
- **Microscopie confocale, FRAP, FLIP, FRET** pour analyser la dynamique subcellulaire.  
- Définition d’un modèle de régulation fine des interactions P-TEFb / partenaires.

### **Scripts associés**
- `scripts/thesis/confocal_processing.py`  
- `scripts/thesis/frap_analysis.py`

---

# **2. Postdoc CONICET — Buenos Aires (2013–2015)**

## **2.1. Usage évolutif des acides aminés (MBE 2014)**

<p align="center">
  <img src="figures/aa_usage.png" width="520px">
</p>

### **Objectif scientifique**  
Quantifier comment le coût biosynthétique, la stabilité chimique et la diversité séquentielle contraignent l’usage des acides aminés dans les protéomes.

### **Contributions**
- Extraction reproductible de données PaxDB pour 17 espèces modèles.  
- Pondération des fréquences d’acides aminés par les niveaux d’abondance.  
- Analyse multi-objectifs (coût ↔ stabilité ↔ diversité).  
- Comparaison inter-espèces et modélisation évolutive du compromis protéomique.  
- Contribution aux figures et à la rédaction de l’article (MBE 2014).

### **Scripts associés**
- `scripts/amino_acids/load_paxdb_data.py`  
- `scripts/amino_acids/aa_composition_analysis.R`

---

## **2.2. Structure & dynamique des répétitions Ankyrine (PLOS Comp Biol 2015)**

<p align="center">
  <img src="figures/ankyrin_energy_landscape.png" width="520px">
</p>

### **Objectif scientifique**  
Comprendre la modularité structurale des répétitions Ankyrine et les déterminants énergétiques du repliement coopératif.

### **Contributions**
- Alignements structuraux sur de larges ensembles de répétitions.  
- Comparaison de trois familles d’HMM (Pfam générique, modèle structuro-dérivé, HMM segmentés).  
- Tessellation protéique pour cartographier les zones de stabilité/déstabilisation.  
- Participation à la rédaction et validation des analyses (PLOS Comp Biol 2015).

### **Scripts associés**
- `scripts/ankyrin/structural_alignment.py`  
- `scripts/ankyrin/hmm_comparison.R`

---

## **2.3. Modularité fonctionnelle Ankyrine–SLiMs (projet personnel)**

<p align="center">
  <img src="figures/ankyrin_modularity.png" width="520px">
</p>

### **Objectif scientifique**  
Identifier des modules fonctionnels combinant domaines Ankyrine et motifs linéaires désordonnés (SLiMs).

### **Contributions**
- Pipeline complet automatisé : récupération → parsing → nettoyage → détection SLiMs.  
- Alignements HMM / conservation évolutive (HMMER, CD-HIT, BLAST).  
- Enrichissements fonctionnels (GO, domaines).  
- Figures Python/R et visualisations structurales.

### **Scripts associés**
- `scripts/ankyrin_modularity/fetch_sequences.py`  
- `scripts/ankyrin_modularity/slim_detection.py`

---

# **3. Ingénierie logicielle — Airbus Defence & Space (2017–2019)**  
## **Développement et maintenance de systèmes critiques**

<p align="center">
  <img src="figures/airbus_systems.png" width="520px">
</p>

### **Objectif scientifique / technique**  
Garantir la stabilité et la traçabilité de systèmes critiques de gestion de configuration utilisés dans l’aérospatial.

### **Contributions**
- Développement et maintenance de modules Java JEE, Oracle, Windchill.  
- Scripts d’automatisation (shell) pour chaînes d’intégration interne.  
- Migration et refactorisation de systèmes complexes.  
- Documentation complète, normes industrielles, cycle en V / Agile.  

### **Scripts associés**
- `scripts/airbus/config_migration.sh`  
- `scripts/airbus/windchill_automation.py`

---

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

# **5. Ressources transversales**

- `environments/conda_bulk.yml` – environnement RNAseq  
- `environments/conda_singlecell.yml` – environnement single-cell  
- `environments/conda_structural.yml` – environnement bioinfo structurale  

---

# **🔚 Fin du document**

Merci de votre lecture.  
Je reste disponible pour toute question scientifique ou collaboration.



























# 🧬 Portfolio scientifique – Nina Verstraete

Bienvenue sur mon espace de recherche et de développement.  
Ce dépôt rassemble mes travaux en **biologie moléculaire, bioinformatique, modélisation et ingénierie logicielle**, menés entre 2008 et 2023 dans différents contextes académiques et industriels.

---

## 🧭 Parcours scientifique et technique

| Période | Institution | Domaine principal | Dossier |
|----------|--------------|------------------|----------|
| **2008–2012** | CNRS / ENS Paris | Biologie moléculaire, régulation transcriptionnelle | [`01_PhD_CNRS`](./01_PhD_CNRS) |
| **2013–2015** | CONICET – Univ. de Buenos Aires | Bioinformatique structurale et évolutive | [`02_Postdoc_CONICET`](./02_Postdoc_CONICET) |
| **2017–2019** | Capgemini / Airbus Defence & Space | Ingénierie logicielle, systèmes critiques | [`03_Industry_AIRBUS`](./03_Industry_AIRBUS) |
| **2020–2023** | INSERM – CRCT Toulouse | Bioinformatique en cancérologie et immunologie | [`04_Postdoc_INSERM`](./04_Postdoc_INSERM) |

---

## 🧩 Domaines d’expertise

- **Biologie moléculaire & cellulaire :** transcription, signalisation, microscopie, mutagenèse dirigée  
- **Bioinformatique :** RNAseq (bulk & single-cell), analyses multi-omics, modélisation d’écosystèmes tumoraux  
- **Programmation & automatisation :** Python, R, Bash, Snakemake, OpenMOLE, NetLogo  
- **Modélisation et simulation :** modèles multi-agents, réseaux dynamiques, analyses de stabilité  
- **Reproductibilité scientifique :** gestion de pipelines, documentation FAIR, visualisation des données  
- **Ingénierie logicielle :** Java JEE, SQL, XML, CI/CD, Scrum, documentation technique  

---

## 🧱 Structure du dépôt

├── 01_PhD_CNRS/ ← Thèse CNRS/ENS Paris
│ └── transcription_regulation_hexim_tat.md
│
├── 02_Postdoc_CONICET/ ← Bioinformatique structurale (Buenos Aires)
│ ├── aa_usage_evolution.md
│ ├── ankyrin_structure_dynamics.md
│ └── ankyrin_modularity.md
│
├── 03_Industry_AIRBUS/ ← Ingénierie logicielle & systèmes critiques
│ └── README.md
│
├── 04_Postdoc_INSERM/ ← Bioinformatique et modélisation tumorale (Toulouse)
│ ├── tumor_ecosystem_modeling.md
│ ├── macrophage_polarization.md
│ ├── immunotherapy_prediction.md
│ ├── tumor_microenvironment_LungPredict.md
│ ├── covid_network_medicine.md
│ ├── sarscov2_systemic_effects.md
│ └── clonal_dynamics_scRNAseq.md
│
├── figures_visuals/ ← Schémas, diagrammes, visualisations
└── methods_and_templates/ ← Outils de reproductibilité (workflows, checklists, environnements)


---

## 🧬 Focus thématique

### 🔹 Recherche académique
Exploration de la régulation génique, des contraintes structurales et des réseaux multi-échelles :  
- Modélisation multi-agents du microenvironnement tumoral  
- Prédiction de la réponse à l’immunothérapie (RNAseq, GEMDECAN)  
- Étude des contraintes évolutives sur les acides aminés  
- Détection de motifs fonctionnels dans les protéines répétées

### 🔹 Ingénierie et pipelines
- Conception de workflows reproductibles (Snakemake, OpenMOLE, Nextflow)  
- Déploiement d’environnements conda / Docker  
- Automatisation des traitements RNAseq et scRNAseq  
- Visualisation scientifique avec Python, R et matplotlib

---

## 🔗 Ressources

- **Visualisations :** [`figures_visuals`](./figures_visuals)  
- **Méthodes et outils FAIR :** [`methods_and_templates`](./methods_and_templates)  
- **Publications :** disponibles dans les fichiers projets correspondants  
- **Contact :** nina.verstraete [at] gmail.com  

---

### 🧭 Frise chronologique simplifiée

2008 ────▶ 2012 CNRS / ENS Paris → Biologie moléculaire
│
2013 ────▶ 2015 CONICET Buenos Aires → Bioinformatique structurale
│
2017 ────▶ 2019 Airbus / Capgemini → Ingénierie logicielle
│
2020 ────▶ 2023 INSERM / CRCT → Bioinformatique translationnelle


---

### ✳️ Objectif de ce dépôt

Ce portfolio vise à rendre visibles mes travaux interdisciplinaires à la croisée de la **biologie**, de la **modélisation numérique** et du **développement logiciel**, dans une logique d’ouverture, de pédagogie et de reproductibilité.
