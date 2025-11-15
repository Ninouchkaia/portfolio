# R Scripts Overview

This directory contains the R scripts used in the differential abundance and
visualization steps of the barcode drug screen analysis.

Two main scripts are used in the modern pipeline:

---

## 1. `deseq2_differential_abundance.R`

Runs DESeq2 on the full dataset (`~ exp + condition`) and produces:

- one differential abundance file per condition vs control  
- a merged **log2 fold-change matrix** (barcodes × conditions)

### **Inputs**

**`results/deseq2_inputs/counts_for_deseq2.tsv`**

- First column: `barcode`  
- Following columns: counts per sample  
- Column names must match `sample_id` in the design file

**`results/deseq2_inputs/design_for_deseq2.tsv`**

Required columns:

- `sample_id`  
- `exp`  
- `condition` (must include `"control"` as reference)

### **Outputs**

- `results/deseq2/deseq2_<condition>_vs_control.tsv`  
  *(one file per treatment condition)*

- `results/deseq2/deseq2_log2fc_matrix.tsv`  
  - column 1: `barcode`  
  - other columns: log2FC values for all conditions

---

## 2. `deseq2_pheatmap_signatures.R`

Reads the merged log2FC matrix and produces a high-quality heatmap similar
to Figure 3 of the manuscript.

### **Inputs**

**`results/deseq2/deseq2_log2fc_matrix.tsv`**  
Output of `deseq2_differential_abundance.R`.

**`results/annotations/colnames_annotated_2023.csv`** *(optional)*  
Annotation file used for column annotation bars.

Minimal format:

```
condition;DrugClass
Gefitinib_050u_exp151121;EGFR_inhibitor
Osimertinib_050u_exp151121;EGFR_inhibitor
...

```

### **Output**

- `results/figures/fig3_barcode_signatures_heatmap.png`

---

# Legacy Scripts

Older preprocessing scripts used during exploratory phases of the analysis are
stored in:

```
scripts/legacy/

```

These include:

- early DESeq2-per-experiment runners  
- old correlation scripts  
- barcode cleaning utilities  
- historical PCA / QC scripts  

They are preserved for reproducibility but are **not part of the modern pipeline**.
```

---
---
---

1. `scripts/deseq2_differential_abundance.R`
   → lance DESeq2 sur tout le jeu de données (`~ exp + condition`) et sort :

   * un fichier **par condition vs control**
   * **une matrice large** de log2FC (barcodes × conditions)

2. `scripts/deseq2_pheatmap_signatures.R`
   → lit la matrice de log2FC et fait la heatmap “Fig.3-like”.

### Inputs attendus pour deseq2_differential_abundance.R

* `results/deseq2_inputs/counts_for_deseq2.tsv`

  * 1ère colonne : `barcode` (ou équivalent)
  * colonnes suivantes : counts par échantillon (noms = `sample_id`)

* `results/deseq2_inputs/design_for_deseq2.tsv`
  Colonnes minimales :

  * `sample_id`
  * `exp`
  * `condition` (avec un niveau `"control"`)

### Outputs clés

* `results/deseq2/deseq2_<condition>_vs_control.tsv` (un par condition)
* `results/deseq2/deseq2_log2fc_matrix.tsv`

  * 1ère colonne : `barcode`
  * colonnes suivantes : log2FC par condition vs contrôle

---

### Inputs attendus pour deseq2_pheatmap_signatures.R

* `results/deseq2/deseq2_log2fc_matrix.tsv`
  (produit par le script précédent)

* `results/annotations/colnames_annotated_2023.csv` *(optionnel)*
  Format minimal :

  ```text
  condition;DrugClass
  Gefitinib_050u_exp151121;EGFR_inhibitor
  Osimertinib_050u_exp151121;EGFR_inhibitor
  ...
  ```

### Output

* `results/figures/fig3_barcode_signatures_heatmap.png`

---