#!/usr/bin/env Rscript

# Differential abundance analysis on barcodes
# using DESeq2 with design ~ exp + condition

suppressPackageStartupMessages({
  library(DESeq2)
  library(apeglm)
  library(readr)
  library(dplyr)
  library(fs)
})

# ---- Paths ----
counts_file <- "results/deseq2_inputs/counts_for_deseq2.tsv"
design_file <- "results/deseq2_inputs/design_for_deseq2.tsv"
output_dir  <- "results/deseq2"

dir_create(output_dir, recurse = TRUE)

message("Reading counts: ", counts_file)
cts <- read_tsv(counts_file)

message("Reading design: ", design_file)
coldata <- read_tsv(design_file)

# ---- Prepare DESeq2 objects ----

# 1st column = barcode id
if (!"barcode" %in% tolower(names(cts)[1])) {
  message("Assuming first column is barcode ID.")
}

barcode_ids <- cts[[1]]
cts_mat <- as.matrix(cts[, -1])
rownames(cts_mat) <- barcode_ids

# Design must have a column with sample IDs matching colnames of counts
# Suppose the column is called 'sample'
if (!"sample" %in% names(coldata)) {
  stop("Design file must contain a column named 'sample' matching column names of counts.")
}

if (!all(colnames(cts_mat) == coldata$sample_id)) {
  stop("Column names of counts do not match 'sample_id' in design file (same order expected).")
}

coldata$condition <- factor(coldata$condition)
coldata$exp       <- factor(coldata$exp)

message("Building DESeqDataSet...")
dds <- DESeqDataSetFromMatrix(
  countData = round(cts_mat),
  colData   = coldata,
  design    = ~ exp + condition
)

# Relevel control
if (!"control" %in% levels(dds$condition)) {
  stop("No level named 'control' found in 'condition' factor. Please check design file.")
}
dds$condition <- relevel(dds$condition, ref = "control")

message("Running DESeq()...")
dds <- DESeq(dds)

message("Available coefficients:")
print(resultsNames(dds))

# ---- Extract results for all conditions vs control ----

res_list <- list()

for (coef_name in resultsNames(dds)) {
  # We keep only contrasts of the form condition_<COND>_vs_control
  if (grepl("^condition_", coef_name) && grepl("_vs_control$", coef_name)) {
    cond <- sub("^condition_", "", coef_name)
    cond <- sub("_vs_control$", "", cond)
    message("Processing contrast: ", coef_name, " (condition = ", cond, ")")
    
    res <- lfcShrink(dds, coef = coef_name, type = "apeglm")
    res_df <- as.data.frame(res)
    res_df$barcode <- rownames(res_df)
    
    # Save one file per contrast
    out_file <- file.path(output_dir, paste0("deseq2_", cond, "_vs_control.tsv"))
    write_tsv(res_df, out_file)
    message("Written: ", out_file)
    
    # Keep only log2FC for the merged matrix
    res_list[[cond]] <- res_df %>%
      select(barcode, log2FoldChange) %>%
      rename(!!cond := log2FoldChange)
  }
}

# ---- Build wide log2FC matrix (barcodes × conditions) ----

if (length(res_list) == 0) {
  stop("No contrasts of the form condition_*_vs_control were found in resultsNames(dds).")
}

message("Merging log2FC across conditions...")
merged <- Reduce(function(x, y) full_join(x, y, by = "barcode"), res_list)

# Clean and sort
logfc_mat <- merged %>%
  arrange(barcode)

logfc_file <- file.path(output_dir, "deseq2_log2fc_matrix.tsv")
write_tsv(logfc_mat, logfc_file)
message("Written merged log2FC matrix: ", logfc_file)

message("Done.")
