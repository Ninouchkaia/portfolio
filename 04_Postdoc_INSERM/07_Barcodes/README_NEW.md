# Barcode-based clonal dynamics analysis

Raw sequencing data from 11 runs (520 samples: 40 controls, 12 time-zero, 468 drug-treated conditions) were parsed and merged into a single barcode–sample count matrix using custom Python scripts. Only barcodes with at least 1 read detected in ≥5 control samples and ≥5 time-zero samples were kept, resulting in 12 305 high-confidence barcodes. Missing values (barcodes absent from a given sample) were treated as 0 reads.

To avoid numerical issues in downstream analyses while preserving count structure, zero counts were replaced by a small offset (0.01 reads) and all counts were multiplied by 100, producing an integer count matrix (combined_runs_filtered_avoid_zero_reads_all_replaced_x100.csv) used as input for DESeq2. A sample annotation table (design_6_2023.csv) was generated in Python from the column names, encoding run, experiment, replicate and condition (drug_dose vs control).

Differential barcode abundance between each treatment and its matched control was quantified with DESeq2 (R) using a negative binomial GLM. For the global analysis, we used a design ~ exp + condition, with “control” as reference. For per-experiment analyses, DESeq2 was run separately on each experiment-specific count matrix with design ~ condition, again using the experiment’s control as reference level. For each contrast, we computed log2 fold-changes (log2FC) and p-values, and applied apeglm shrinkage of log2FC to stabilize estimates for low-count barcodes.

Shrunk log2FC values from all experiments were merged into a single barcode × condition matrix (merged_unfiltered_deseq2_with_expdate_2023.csv). This matrix was then:

- filtered on statistical significance (typically p-value < 0.05 with optional |log2FC| thresholds depending on the analysis),

- used to build a barcode signature heatmap (barcodes in rows, conditions in columns) with hierarchical clustering of conditions and drug annotations (mechanism of action, e.g. EGFR inhibitors, chemotherapies, epigenetic modulators),

- converted into a drug–drug correlation matrix (Pearson correlation on log2FC profiles) for (i) a clustered correlation heatmap and (ii) a force-directed drug similarity network (edges drawn for |r| > 0.8, node colours reflecting drug class).











# Barcode drug screening – Pipeline Python + DESeq2

## Résumé

Ce pipeline reproduit l’analyse décrite dans la section “Barcoding computational analysis” :

- Fusion des 11 runs de séquençage (~520 échantillons).
- Filtrage des barcodes :
  - barcodes détectés dans **< 5 contrôles** et/ou **< 5 Temps0** sont retirés.
- Normalisation within-sample à **1e6** reads par échantillon (CPM-like).
- Analyse différentielle par condition vs contrôle avec **DESeq2** (R).
- Construction de signatures de log2FC par condition.
- Corrélations drogue–drogue et réseau de drogues (Pearson).

---

## Étape 1 – Préprocessing / filtrage

```bash
python scripts/01_preprocess_barcodes.py \
  --input_dir data/raw \
  --pattern "Run*.csv" \
  --output_prefix data/processed/barcodes \
  --min_reads 1 \
  --min_controls 5 \
  --min_timezeros 5
````

Sorties principales :

* `data/processed/barcodes_combined_raw_counts.csv`
* `data/processed/barcodes_filtered_counts.csv`
* `data/processed/barcodes_filtered_cpm.csv`

---

## Étape 2 – QC des contrôles

```bash
python scripts/02_qc_controls_variability.py \
  --counts data/processed/barcodes_filtered_counts.csv \
  --output_prefix results/qc_controls/controls_variability
```

Sorties :

* `results/qc_controls/controls_variability_per_barcode.tsv`
* `results/qc_controls/controls_variability_violin.png`

Permet de vérifier la stabilité des contrôles via la métrique `(max-min)/moyenne`.

---

## Étape 3 – Préparation des inputs DESeq2

```bash
python scripts/03_build_deseq2_inputs.py \
  --counts data/processed/barcodes_filtered_counts.csv \
  --output_dir results/deseq2_inputs
```

Sorties :

* `results/deseq2_inputs/counts_for_deseq2.tsv`
* `results/deseq2_inputs/design_for_deseq2.tsv`

À utiliser dans R / DESeq2 pour calculer log2FC et p-values.

---

## Étape 4 – Fusion des résultats DESeq2 et réseaux

Une fois les résultats DESeq2 exportés (un fichier par condition vs contrôle, avec une colonne `log2FoldChange`) dans `results/deseq2/` :

```bash
python scripts/04_correlations_and_networks.py \
  --deseq2_dir results/deseq2 \
  --pattern "*.tsv" \
  --output_prefix results/networks/drug_signatures \
  --corr_threshold 0.8
```

Sorties :

* `results/networks/drug_signatures_logfc_matrix.csv`
* `results/networks/drug_signatures_correlation_matrix.csv`
* `results/networks/drug_signatures_correlation_heatmap.png`
* `results/networks/drug_signatures_network.png`

Ces résultats correspondent aux Figures 3–5 décrites dans le manuscrit : clustering des signatures, heatmap de corrélations, et réseau de drogues.

---


### Barcoding analysis (drug screen)

- `notebooks/01_visualize_barcode_signatures.ipynb`  
  End-to-end visualization of barcode signatures (Figures 3–5).
- `scripts/05_generate_figures_for_paper.py`  
  Regenerates the paper figures from processed DESeq2 outputs.



## Flowchart
  
```mermaid
graph TD;;

%% =========================
%% RAW DATA
%% =========================
A0([FASTQ sequencing runs<br>11 runs, 520 samples]) --> A1

A1[Python: merge all runs<br><code>combined_runs.csv</code>] --> A2

%% =========================
%% BARCODE QC + FILTER
%% =========================

subgraph B[Barcode filtering & QC]
    A2 --> B1[Python: remove junk barcodes<br>keep barcodes detected in ≥5 controls AND ≥5 time-zero samples<br><code>drop_junk_barcodes_2023.py</code>]
    B1 --> B2[Output:<br><code>combined_runs_filtered.csv</code>]
end

B2 --> C1

%% =========================
%% ZERO-READ HANDLING
%% =========================

subgraph C[Handling zero counts]
    C1[Replace missing reads with 0] --> C2[Python: replace zeros with small offset (0.01)<code>avoid_zero_reads.py</code>]
    C2 --> C3[Alternative: replace ALL zeros (not only problematic ones)<code>avoid_zero_reads_by_replacing_all_zeros.py</code>]
    C3 --> C4[Multiply all counts ×100 to avoid decimals<br><code>multiply_counts_by_100.py</code>]
    C4 --> C5[Output:<code>combined_runs_filtered_avoid_zero_reads_all_replaced_x100.csv</code>]
end

C5 --> D1

%% =========================
%% OPTIONAL NORMALIZATION
%% =========================
subgraph D[Optional normalization]
    D1[Python: per-sample reads → million scaling<br><code>normalize_per_millions_2023.py</code>]
    D1 -->|used for PCA & QC only| D2[Output:<br><code>combined_runs_filtered_avoid_zero_reads_scaled.csv</code>]
end

C5 --> E1

%% =========================
%% ANNOTATION
%% =========================
subgraph E[Annotation]
    E1[Python: build design table<br><code>build_design_2023.py</code>] --> E2
    E2[Output:<br><code>design_6_2023.csv</code><br>(run, exp, replicate, condition)]
end

E2 --> F1

%% =========================
%% DESEQ2
%% =========================

subgraph F[Differential abundance (R / DESeq2)]
    F1[Global DESeq2:<br>design = ~ exp + condition<br>ref = control<br><code>deseq2_script.R</code>] --> F2
    F1 --> F3

    F2[Compute shrunk log2FC (apeglm)<br><code>lfcShrink</code> per condition] --> F4

    F3[Per-experiment DESeq2 (11 experiments):<br>design = ~ condition<br>ref = control or CtrlMs_000u<br><code>deseq2_script.R</code>] --> F4

    F4[Output: many DESeq2 result CSV files<br>(one per condition per experiment)]
end

F4 --> G1

%% =========================
%% MERGE DESEQ2 RESULTS
%% =========================

subgraph G[Aggregation & filtering]
    G1[Python/R: merge all DESeq2 results<br><code>merged_unfiltered_deseq2_with_expdate_2023.csv</code>] --> G2
    G2[Filter on p-value / remove NA<br><code>merged_logfc_pval_filtered_deseq2_2023.csv</code>] --> G3
    G3[Fill missing logFC (optional)<br><code>merged_logfc_pval_filtered_deseq2_2023_fillna.csv</code>]
end

G3 --> H1

%% =========================
%% VISUALIZATION ANALYSES
%% =========================

subgraph H[Visualization & high-level analysis]
    H1[Heatmap of barcode signatures<br>(barcodes × conditions)<br>Fig.3<br><code>pheatmap_2023.R</code>] --> H2

    H2[Correlation matrix of conditions<br>(Pearson on log2FC)<br>Fig.4<br><code>correl_deseq2_2023.py</code>] --> H3

    H3[Drug similarity network<br>(edges r ≥ 0.8)<br>Fig.5<br><code>networkx</code> in notebook]

    H1 --> H4[PCA QC (runs/experiments/conditions)<br><code>FactoMineR</code> / <code>prcomp</code>]
end

```



