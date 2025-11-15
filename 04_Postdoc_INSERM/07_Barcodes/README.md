
```
2022_Barcodes/
├── data/
│   ├── raw/              # Run*.csv sortis du séquenceur
│   └── processed/        # matrices fusionnées / filtrées
├── results/
│   ├── qc_controls/      # QC des contrôles
│   ├── deseq2_inputs/    # inputs pour R/DESeq2
│   └── networks/         # corrélations / réseaux
├── scripts/
│   ├── 01_preprocess_barcodes.py
│   ├── 02_qc_controls_variability.py
│   ├── 03_build_deseq2_inputs.py
│   └── 04_correlations_and_networks.py
└── README_barcode_pipeline.md
```


```R
library(DESeq2)
counts <- read.table("results/deseq2_inputs/counts_for_deseq2.tsv",
                     header=TRUE, row.names=1, sep="\t", check.names=FALSE)
design <- read.table("results/deseq2_inputs/design_for_deseq2.tsv",
                     header=TRUE, sep="\t", stringsAsFactors=TRUE)

# Exemple : analyse par expérience
for (xp in unique(design$exp)) {
  cat("### Experience", xp, "\n")
  idx <- design$exp == xp & !design$is_timezero  # à adapter si tu veux utiliser time0
  dds <- DESeqDataSetFromMatrix(
    countData = counts[, idx],
    colData   = design[idx, ],
    design    = ~ condition
  )
  dds <- DESeq(dds)

  # Comparaisons vs controle
  cond_levels <- levels(design$condition)
  cond_levels <- cond_levels[cond_levels != "control" & cond_levels != "time0"]

  for (cond in cond_levels) {
    res <- results(dds, contrast=c("condition", cond, "control"))
    out_file <- paste0("results/deseq2/", xp, "_", cond, "_vs_control.tsv")
    dir.create("results/deseq2", showWarnings=FALSE, recursive=TRUE)
    write.table(as.data.frame(res), file=out_file, sep="\t", quote=FALSE)
  }
}
```