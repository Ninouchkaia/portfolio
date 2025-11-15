# LungPredict — TF Analysis Pipeline

This directory implements a structured and reproducible workflow for analysing the LungPredict cohort at the transcriptional and transcription factor (TF) activity level.

The pipeline integrates:
- unified clinical + molecular annotations,
- gene expression (TPM),
- pathway- and expression-based clustering,
- TF activity inference using DoRothEA + viper (TCGA & GTEx regulons),
- annotated heatmaps,
- regulon comparison,
- PCA.

It is organised into a single master script and several small utility modules.

---
```text
## 1. Directory Layout
TF_analysis/
├── data/
│ ├── LP_FFPE_STAR_RSEM_TPM.txt
│ ├── clinic_data_v2_clean.csv
│ ├── ReactomeClustersAllPatients.csv
│ ├── DeconvCancerClusters.txt
│ ├── expression_based_patient_clusters_resCutCorrected.txt
│ └── full_annotations_with_clusters_corrected1.txt
│
├── scripts/
│ ├── lungpredict_tf_pipeline.R
│ ├── utils_logging.R
│ ├── utils_annotations.R
│ ├── utils_heatmaps.R
│ └── utils_tf_activity.R
│
├── results/
│ ├── heatmaps_expression/
│ ├── heatmaps_TF/
│ ├── comparisons_tcga_gtex/
│ ├── pca/
│ └── logs/
```

All inputs go in `data/`; all derived outputs are saved under `results/`.

---

## 2. Required Input Data

| File | Description |
|------|-------------|
| **LP_FFPE_STAR_RSEM_TPM.txt** | TPM matrix, Gene × Patient |
| **clinic_data_v2_clean.csv** | Clinical + mutational annotations |
| **ReactomeClustersAllPatients.csv** | Pathway-level clustering |
| **DeconvCancerClusters.txt** | Deconvolution-based clustering |
| **expression_based_patient_clusters_resCutCorrected.txt** | Final expression clusters |
| **full_annotations_with_clusters_corrected1.txt** | Final merged annotation table |

The unified annotation table aggregates clinical and clustering metadata for all downstream visualisations.

---

## 3. Workflow Overview

### Step 1 — Load expression + annotation  
The expression matrix and unified annotations are loaded and synchronised.

### Step 2 — Optional clinical filtering  
Filtering criteria (diagnosis, location, KRAS/EGFR/STK11 status, metastatic state…) can be applied directly through parameters.

### Step 3 — Expression analysis  
- expression heatmap (gene × patient),  
- correlation matrix (patient × patient).  

### Step 4 — TF activity inference  
TF activities are computed using:
- **DoRothEA TCGA** regulons  
- **DoRothEA GTEx** regulons  
via `viper`.

### Step 5 — TF-level visualisation  
Heatmaps:
- TF × patient (all TFs),
- TF × patient (filtered by variance),
- TF-based patient correlation.

### Step 6 — TCGA vs GTEx comparison  
The two regulon sources are compared at:
- patient level,
- TF level,
- difference matrix level.

### Step 7 — PCA  
Optional PCA is applied to the patient expression profiles.

---

## 4. Pipeline Diagram

```mermaid
graph TD;;
A[Expression TPM] --> C[Synchronise];
B[Unified annotations] --> C;
C --> D[Filtering];
D --> E[Expression heatmaps];
D --> F[TF activity\n(TCGA + GTEx)];
F --> G[TF heatmaps];
F --> H[TCGA vs GTEx\ncomparison];
D --> I[PCA];
```

All inputs go in `data/`; all derived outputs are saved under `results/`.

---

## 2. Required Input Data

| File | Description |
|------|-------------|
| **LP_FFPE_STAR_RSEM_TPM.txt** | TPM matrix, Gene × Patient |
| **clinic_data_v2_clean.csv** | Clinical + mutational annotations |
| **ReactomeClustersAllPatients.csv** | Pathway-level clustering |
| **DeconvCancerClusters.txt** | Deconvolution-based clustering |
| **expression_based_patient_clusters_resCutCorrected.txt** | Final expression clusters |
| **full_annotations_with_clusters_corrected1.txt** | Final merged annotation table |

The unified annotation table aggregates clinical and clustering metadata for all downstream visualisations.

---

## 3. Workflow Overview

### Step 1 — Load expression + annotation  
The expression matrix and unified annotations are loaded and synchronised.

### Step 2 — Optional clinical filtering  
Filtering criteria (diagnosis, location, KRAS/EGFR/STK11 status, metastatic state…) can be applied directly through parameters.

### Step 3 — Expression analysis  
- expression heatmap (gene × patient),  
- correlation matrix (patient × patient).  

### Step 4 — TF activity inference  
TF activities are computed using:
- **DoRothEA TCGA** regulons  
- **DoRothEA GTEx** regulons  
via `viper`.

### Step 5 — TF-level visualisation  
Heatmaps:
- TF × patient (all TFs),
- TF × patient (filtered by variance),
- TF-based patient correlation.

### Step 6 — TCGA vs GTEx comparison  
The two regulon sources are compared at:
- patient level,
- TF level,
- difference matrix level.

### Step 7 — PCA  
Optional PCA is applied to the patient expression profiles.

---

## 4. Pipeline Diagram

```mermaid
graph TD;;
A[Expression TPM] --> C[Synchronise];
B[Unified annotations] --> C;
C --> D[Filtering];
D --> E[Expression heatmaps];
D --> F[TF activity\n(TCGA + GTEx)];
F --> G[TF heatmaps];
F --> H[TCGA vs GTEx\ncomparison];
D --> I[PCA];
```



## 5. Running the Pipeline (RStudio)

Set the working directory to TF_analysis/:

source("scripts/lungpredict_tf_pipeline.R")


Everything will run end-to-end:

expression heatmaps,

TF activities (TCGA + GTEx),

TF visualisations,

regulon comparisons,

PCA,

logs.

Outputs are written under results/.

6. Script Overview
lungpredict_tf_pipeline.R

Main orchestrator.
Defines paths, loads utilities, executes the complete workflow.

utils_logging.R

Light logging system.

init_logger()

info_print()

debug_print()

utils_annotations.R

Handles:

annotation loading,

sample synchronisation,

clinical filtering,

mutation flags,

age categorisation.

utils_heatmaps.R

Creates publication-ready heatmaps using ComplexHeatmap:

gene × patient,

patient × patient correlations,

TF × patient (with annotation bars).

utils_tf_activity.R

Provides:

TF activity estimation with viper,

regulon selection (A/B confidence),

variance filtering,

TCGA vs GTEx comparison.




## 7. Outputs (results/)
Expression heatmaps

results/heatmaps_expression/

*_gene_by_patient.png

*_patient_correlation.png

TF heatmaps

results/heatmaps_TF/

TF × patient (all TFs)

TF × patient (filtered)

TF correlation

TCGA–GTEx comparison

results/comparisons_tcga_gtex/

patient correlation

TF correlation

difference matrices

PCA

results/pca/

scree plot

PCA plots

variance summary

Logs

results/logs/pipeline.log


## 8. Dependencies

Install:

install.packages(c("tidyverse","ComplexHeatmap","circlize","RColorBrewer","viridis"))
BiocManager::install(c("viper","dorothea","Hmisc"))

## 9. Summary

This pipeline provides:

a unified annotation and expression framework,

gene- and TF-level characterisation of patients,

regulon-based TF inference (TCGA + GTEx),

annotated heatmaps and PCA,

modular and reproducible execution via a single pipeline script.


It enables consistent multi-layer molecular analysis across the LungPredict cohort.



