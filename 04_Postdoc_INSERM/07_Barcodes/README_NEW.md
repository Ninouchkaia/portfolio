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

```

---
