# **📘 Barcoding Drug Screen — Computational Pipeline**

## **Overview**

This repository reproduces the full computational workflow described in the “Barcoding computational analysis” section of the manuscript.
The pipeline processes 520 samples across 11 sequencing runs, filters barcodes, normalizes counts, performs differential abundance with DESeq2, and generates drug-response signatures (log2FC), correlation maps, and drug–drug similarity networks.

### **Main steps**

1. **Merge raw runs** into a single barcode × sample matrix.
2. **Filter unreliable barcodes**

   * keep barcodes detected in **≥5 controls** and **≥5 time-zero samples**.
3. **Handle zeros** (biological absences vs. technical dropouts): replace technical zeros by a small offset (0.01).
4. **Normalize** each sample to a total of **1e6 reads** (CPM-like).
5. **Build design tables** (run, experiment, replicate, condition).
6. **DESeq2 analysis**

   * log2 fold changes (per condition vs control),
   * optional p-value filtering and shrinkage.
7. **Construct log2FC signatures** (conditions × barcodes).
8. **Downstream analyses**

   * heatmaps (drug clusters),
   * correlation matrix (drug–drug similarity),
   * drug similarity network (edges: Pearson ≥ 0.8).
9. **Notebook visualization** of Figures 3–5.

---

# **1. Preprocessing & Barcode Filtering**

```bash
python scripts/01_preprocess_barcodes.py \
  --input_dir data/raw \
  --pattern "Run*.csv" \
  --output_prefix data/processed/barcodes \
  --min_reads 1 \
  --min_controls 5 \
  --min_timezeros 5
```

**What this does**

* merges all runs into `combined_raw_counts.csv`
* removes barcodes absent from *≥5 control* and *≥5 time-zero* samples
* outputs filtered matrices

**Outputs**

* `data/processed/barcodes_combined_raw_counts.csv`
* `data/processed/barcodes_filtered_counts.csv`
* `data/processed/barcodes_filtered_cpm.csv`

---

# **2. QC of Controls**

```bash
python scripts/02_qc_controls_variability.py \
  --counts data/processed/barcodes_filtered_counts.csv \
  --output_prefix results/qc_controls/controls_variability
```

This computes barcode-level control stability:

[
\frac{\max - \min}{\text{mean}}
]

**Outputs**

* `controls_variability_per_barcode.tsv`
* `controls_variability_violin.png`

---

# **3. DESeq2 Inputs (Counts + Design)**

```bash
python scripts/03_build_deseq2_inputs.py \
  --counts data/processed/barcodes_filtered_counts.csv \
  --output_dir results/deseq2_inputs
```

Produces:

* `counts_for_deseq2.tsv`
* `design_for_deseq2.tsv`
  (columns: run, experiment, replicate, condition)

---

# **4. Differential Abundance (DESeq2) — in R**

Run:

```r
source("scripts/deseq2_script.R")
```

This script:

* loads counts + design
* runs DESeq2 with `design = ~ exp + condition` or per-experiment `~ condition`
* shrinks log2FC using **apeglm**
* exports one result CSV per condition

**All DESeq2 outputs should be saved in:**
`results/deseq2/`

---

# **5. Merge DESeq2 Outputs & Build Signatures**

```bash
python scripts/04_correlations_and_networks.py \
  --deseq2_dir results/deseq2 \
  --pattern "*.csv" \
  --output_prefix results/networks/drug_signatures \
  --corr_threshold 0.8
```

Tasks performed:

* merges all log2FC into a single matrix (conditions × barcodes)
* optional p-value filtering
* computes Pearson correlation (condition × condition)
* constructs network (edges: r ≥ threshold)

**Outputs**

* `drug_signatures_logfc_matrix.csv`
* `drug_signatures_correlation_matrix.csv`
* `drug_signatures_correlation_heatmap.png`
* `drug_signatures_network.png`

These correspond to **Figures 3, 4, and 5** in the manuscript.

---

# **6. Notebook Visualization**

* `notebooks/01_visualize_barcode_signatures.ipynb`

This notebook reproduces:

* Figure 3 — logFC heatmap with drug family annotation
* Figure 4 — correlation matrix
* Figure 5 — network graph (Pearson ≥ 0.8)

---

# **7. Pipeline Flowchart**

```mermaid
graph TD;

%% RAW DATA
A0([FASTQ sequencing runs<br>11 runs, 520 samples]) --> A1;
A1[Merge runs → combined_runs.csv] --> B1;

%% BARCODE FILTERING
subgraph B[Barcode filtering]
B1 --> B2[drop_junk_barcodes_2023.py<br>Keep barcodes in ≥5 controls AND ≥5 T0];
B2 --> B3[combined_runs_filtered.csv];
end

%% ZERO HANDLING
subgraph C[Zero-handling]
B3 --> C1[avoid_zero_reads.py<br>replace technical zeros by 0.01];
C1 --> C2[multiply_counts_by_100.py];
C2 --> C3[combined_runs_filtered_avoid_zero_reads_all_replaced_x100.csv];
end

%% NORMALIZATION
C3 --> D1;
subgraph D[Normalization]
D1[normalize_per_millions_2023.py] --> D2[combined_runs_filtered_avoid_zero_reads_scaled.csv];
end

%% DESIGN TABLE
C3 --> E1;
subgraph E[Design table]
E1[build_design_2023.py] --> E2[design_6_2023.csv];
end

%% DESEQ2
E2 --> F1;
subgraph F[DESeq2 (R)]
F1[deseq2_script.R<br>per-experiment & global] --> F2[DESeq2 log2FC tables];
end

%% MERGING
F2 --> G1;
subgraph G[Aggregation]
G1[merge DESeq2 outputs] --> G2[logFC filtering + fillna<br>merged_logfc_pval_filtered_deseq2_2023_fillna.csv];
end

%% VISUALIZATION
G2 --> H1;
subgraph H[Visualization]
H1[Heatmaps → Fig3] --> H2[Corr matrix → Fig4];
H2 --> H3[Drug network → Fig5];
end
```

---
