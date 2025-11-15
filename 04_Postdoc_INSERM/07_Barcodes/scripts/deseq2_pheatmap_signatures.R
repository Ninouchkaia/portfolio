#!/usr/bin/env Rscript

# Heatmap of barcode signatures (Fig.3-like)
# from DESeq2 log2FC matrix

suppressPackageStartupMessages({
  library(readr)
  library(pheatmap)
  library(fs)
})

# ---- Paths ----
logfc_file <- "results/deseq2/deseq2_log2fc_matrix.tsv"
annot_file <- "results/annotations/colnames_annotated_2023.csv"  # optionnel
out_dir    <- "results/figures"
out_file   <- file.path(out_dir, "fig3_barcode_signatures_heatmap.png")

dir_create(out_dir, recurse = TRUE)

message("Reading log2FC matrix: ", logfc_file)
logfc <- read_tsv(logfc_file)

if (ncol(logfc) < 2) {
  stop("log2FC matrix must have at least one condition column.")
}

# First column is barcode
barcode_ids <- logfc[[1]]
mat <- as.matrix(logfc[, -1])
rownames(mat) <- barcode_ids

# ---- Z-score per condition (column-wise scaling) ----
# scale rows -> t(scale(t())) to scale across conditions per barcode
mat_scaled <- t(scale(t(mat)))

# Replace NA (e.g. barcodes with all NA) by 0 for plotting
mat_scaled[is.na(mat_scaled)] <- 0

# ---- Annotation (optionnelle) ----
annotation_col <- NULL

if (file.exists(annot_file)) {
  message("Reading annotation: ", annot_file)
  annot <- read_delim(annot_file, delim = ";", col_types = cols())
  
  # On suppose que la 1ère colonne est le nom de la condition
  cond_colname <- names(annot)[1]
  rownames(annot) <- annot[[1]]
  annot[[1]] <- NULL
  
  # On aligne sur les colonnes de la matrice
  common <- intersect(colnames(mat_scaled), rownames(annot))
  if (length(common) == 0) {
    message("Warning: no common column names between log2FC matrix and annotation. Heatmap will be drawn without annotations.")
  } else {
    annotation_col <- annot[common, , drop = FALSE]
    # Réordonner mat_scaled pour matcher l’annotation
    mat_scaled <- mat_scaled[, common, drop = FALSE]
  }
} else {
  message("No annotation file found, continuing without condition annotations.")
}

# ---- Plot heatmap ----

message("Drawing heatmap: ", out_file)
png(out_file, width = 2400, height = 1800, res = 300)

pheatmap(
  mat_scaled,
  show_rownames   = FALSE,
  fontsize_col    = 6,
  clustering_method = "complete",
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  annotation_col  = annotation_col,
  scale           = "none",
  border_color    = NA
)

dev.off()
message("Done.")
