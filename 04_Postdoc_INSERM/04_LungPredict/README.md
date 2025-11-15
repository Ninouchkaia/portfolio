```text
LungPredict/
└── TF_analysis/
    ├── data/
    │   ├── LP_FFPE_STAR_RSEM_TPM.txt
    │   ├── clinic_data_v2_clean.csv
    │   ├── expression_based_patient_clusters2.txt  (ou a.txt)
    │   ├── ReactomeClustersAllPatients.csv
    │   ├── DeconvCancerClusters.txt
    │   └── full_annotations_with_clusters_corrected1.txt  (output)
    │
    ├── results/
    │   ├── heatmaps/
    │   ├── tf_activities/
    │   └── pca/
    │
    └── scripts/
        └── lungpredict_tf_pipeline.R    ← script principal (avec fonctions)
```

```text
lungpredict/
│
├── data/
│   ├── expression/                     # TPM matrices, filtered versions
│   ├── annotations/                    # clinical, clusters, merged annotation panel
│   └── regulons/                       # dorothea regulons if you want to store them
│
├── scripts/
│   ├── 01_load_and_clean_data.R        # imports, filtering, annotation construction
│   ├── 02_compute_expression_clusters.R
│   ├── 03_compute_tf_activities.R
│   ├── 04_compare_regulon_sources.R    # TCGA vs GTEX
│   ├── 05_plot_heatmaps_expression.R
│   ├── 06_plot_heatmaps_tf.R
│   ├── 07_pca_analysis.R
│   ├── utils_color_annotations.R       # Color schemes + annotation builders
│   └── utils_io.R                      # I/O helpers (get_expression(), save_table(), etc.)
│
├── results/
│   ├── clusters/
│   ├── heatmaps/
│   ├── tf_activities/
│   └── pca/
│
└── README.md
```
