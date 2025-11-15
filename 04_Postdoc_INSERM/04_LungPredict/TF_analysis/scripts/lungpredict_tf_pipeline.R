# scripts/lungpredict_tf_pipeline.R

suppressPackageStartupMessages({
  library(tidyverse)
  library(ComplexHeatmap)
})

source("scripts/utils_logging.R")
source("scripts/utils_annotations.R")
source("scripts/utils_tf_activity.R")
source("scripts/utils_heatmaps.R")

#-------------------------
# Config
#-------------------------

expr_path   <- file.path("data", "LP_FFPE_STAR_RSEM_TPM.txt")
annot_path  <- file.path("data", "full_annotations_with_clusters_corrected1.txt")

out_dir     <- "results"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(out_dir, "heatmaps_expression"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(out_dir, "heatmaps_TF"),         showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(out_dir, "comparisons_tcga_gtex"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(out_dir, "pca"),                 showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(out_dir, "logs"),                showWarnings = FALSE, recursive = TRUE)

init_logger(file.path(out_dir, "logs", "pipeline.log"))
log_info("Pipeline LungPredict TF started.")

#-------------------------
# 1. Chargement des données
#-------------------------

log_info("Chargement expression : ", expr_path)
expression_raw <- load_expression(expr_path)

log_info("Chargement annotations : ", annot_path)
annotations_raw <- load_full_annotations(annot_path)

aligned <- align_expression_annotations(expression_raw, annotations_raw)
expression  <- aligned$expression
annotations <- aligned$annotations

ann_colors <- build_sample_annotation(annotations)

#-------------------------
# 2. Heatmaps basées sur l'expression
#-------------------------

log_info("Heatmap gene x patient.")
plot_gene_expression_heatmap(
  expression,
  ann_colors,
  file = file.path(out_dir, "heatmaps_expression", "expression_gene_by_patient.png")
)

log_info("Heatmap corrélation patients (expression).")
plot_expression_patient_corr(
  expression,
  ann_colors,
  file = file.path(out_dir, "heatmaps_expression", "expression_patient_correlation.png")
)

#-------------------------
# 3. Activités de TF (TCGA)
#-------------------------

log_info("Calcul des activités de TF (TCGA).")
tf_tcga <- compute_tf_activity(expression, source = "tcga")

log_info("Filtrage des TF (sd > 2).")
tf_tcga_filt <- filter_tf_by_sd(tf_tcga, sd_threshold = 2)

log_info("Heatmap TF x patient (non filtré).")
plot_tf_heatmap(
  tf_tcga,
  ann_colors,
  title = "TF activities (TCGA, all TFs)",
  file  = file.path(out_dir, "heatmaps_TF", "TF_tcga_all.png")
)

log_info("Heatmap TF x patient (filtré).")
plot_tf_heatmap(
  tf_tcga_filt,
  ann_colors,
  title = "TF activities (TCGA, sd>2)",
  file  = file.path(out_dir, "heatmaps_TF", "TF_tcga_filtered.png")
)

log_info("Heatmap corrélation patients (TF).")
plot_tf_patient_corr(
  tf_tcga_filt,
  ann_colors,
  file = file.path(out_dir, "heatmaps_TF", "TF_tcga_patient_correlation.png")
)

log_info("Heatmap corrélation TF x TF.")
plot_tf_tf_corr(
  tf_tcga_filt,
  file = file.path(out_dir, "heatmaps_TF", "TF_tcga_TF_correlation.png")
)

log_info("Pipeline LungPredict TF terminé.")


#--------------------------------------------
# 4. Comparaison TCGA vs GTEx
#--------------------------------------------

log_info("Calcul TF activities (GTEx).")
tf_gtex <- compute_tf_activity(expression, source = "gtex")

# Harmonise les TF (garde uniquement TF communs TCGA/GTEx)
common_tfs <- intersect(rownames(tf_tcga), rownames(tf_gtex))
tf_tcga_c <- tf_tcga[common_tfs, , drop = FALSE]
tf_gtex_c <- tf_gtex[common_tfs, , drop = FALSE]

# (1) Corrélation TF × TF entre TCGA et GTEx
log_info("Corrélation TF × TF entre TCGA et GTEx.")
cor_tf_tf <- cor(tf_tcga_c, tf_gtex_c)

png(file.path(out_dir, "comparisons_tcga_gtex", "tcga_gtex_TFcor.png"),
    width = 2000, height = 1800, res = 200)
Heatmap(cor_tf_tf,
        name = "TCGA–GTEx TF corr",
        show_row_names = FALSE,
        show_column_names = FALSE)
dev.off()

# (2) Matrice des différences TF (TCGA − GTEx)
log_info("Matrice des différences TF (TCGA − GTEx).")
diff_mat <- tf_tcga_c - tf_gtex_c

png(file.path(out_dir, "comparisons_tcga_gtex", "tcga_gtex_diff_matrix.png"),
    width = 2000, height = 1800, res = 200)
Heatmap(diff_mat,
        name = "TCGA–GTEx diff",
        show_row_names = FALSE,
        show_column_names = TRUE,
        column_names_gp = gpar(fontsize = 6))
dev.off()


#--------------------------------------------
# 5. PCA (sur l'expression)
#--------------------------------------------

log_info("PCA sur l'expression.")

res_pca <- prcomp(t(expression), scale. = TRUE)

# 5.1 Scree plot
png(file.path(out_dir, "pca", "pca_scree.png"),
    width = 2000, height = 1600, res = 200)
plot(res_pca, type = "l", main = "Scree plot")
dev.off()

# 5.2 Variances expliquées
var_expl <- data.frame(
  PC = paste0("PC", 1:length(res_pca$sdev)),
  Variance = res_pca$sdev^2 / sum(res_pca$sdev^2)
)
write.csv(var_expl,
          file = file.path(out_dir, "pca", "pca_variance.csv"),
          row.names = FALSE)

# 5.3 PCA individus
png(file.path(out_dir, "pca", "pca_individuals.png"),
    width = 2000, height = 1600, res = 200)
plot(res_pca$x[,1], res_pca$x[,2],
     xlab = "PC1", ylab = "PC2",
     pch = 19, col = "black",
     main = "PCA – Individuals")
text(res_pca$x[,1], res_pca$x[,2],
     labels = rownames(res_pca$x), pos = 3, cex = 0.6)
dev.off()

# 5.4 PCA avec annots (couleur par cluster génique si dispo)
if ("GE_Cluster" %in% colnames(annotations)) {
  
  cols <- c("1" = "chocolate1", "2" = "chocolate4", "3" = "chocolate3")
  cluster_col <- cols[as.character(annotations$GE_Cluster)]
  
  png(file.path(out_dir, "pca", "pca_individuals_with_annots.png"),
      width = 2000, height = 1600, res = 200)
  
  plot(res_pca$x[,1], res_pca$x[,2],
       xlab = "PC1", ylab = "PC2",
       pch = 19, col = cluster_col,
       main = "PCA – colored by GE_Cluster")
  
  legend("topright", legend = names(cols), col = cols, pch = 19)
  text(res_pca$x[,1], res_pca$x[,2],
       labels = rownames(res_pca$x), pos = 3, cex = 0.6)
  
  dev.off()
}

