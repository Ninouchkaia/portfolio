# Predicting Response to Immunotherapy (GEMDECAN Project)

**Affiliation:** INSERM U1037 – Centre de Recherches en Cancérologie de Toulouse (CRCT)  
**Period:** 2020–2023  
**Publication:** [bioRxiv, 2021](https://www.biorxiv.org/content/10.1101/2021.04.09.439207v1)  

---

## Context
Despite major clinical advances, only a fraction of cancer patients respond to immunotherapy.  
The GEMDECAN project aimed to identify gene expression signatures predictive of therapeutic response by integrating bulk RNA-seq data with clinical outcomes and molecular features across cancer types.

The goal was to define reproducible biomarkers guiding patient stratification for immune checkpoint therapies.

---

## Objectives
- Build a standardized analysis pipeline for RNA-seq data from tumor samples.  
- Identify and validate gene signatures correlated with immunotherapy response.  
- Explore functional categories and immune-related pathways underlying predictive features.  

---

## Methods
- **Data:** Bulk RNA-seq datasets from public repositories (TCGA, ICGC, clinical cohorts).  
- **Pipeline design:** Quality control (FastQC, MultiQC), alignment (STAR), quantification (featureCounts), normalization and batch correction.  
- **Analysis:** Differential expression (DESeq2), gene set enrichment, correlation with clinical response scores.  
- **Modeling:** Machine learning approaches (lasso regression, random forest) for predictive feature selection.  
- **Visualization:** Volcano plots, ROC curves, heatmaps of signature genes.  

---

## Contributions
- Co-designed the RNA-seq analysis pipeline and automated preprocessing workflow.  
- Performed expression modeling and differential analysis linking gene expression to therapy response.  
- Contributed to functional interpretation and visual representation of predictive gene sets.  
- Participated in manuscript preparation and figure design.  

---

## Reference

*Xie T.*, et al. *GEM-DeCan: Improving tumor immune microenvironment profiling by the integration of novel gene expression and DNA methylation deconvolution signatures.*  
*bioRxiv*, 2021. [DOI:10.1101/2021.04.09.439207v1](https://www.biorxiv.org/content/10.1101/2021.04.09.439207v1)

---
 
## Deconvolution pipeline

See further details on the deconvolution scripts (R and Python) and the pipeline workflow [here](README_DETAILS.md).





