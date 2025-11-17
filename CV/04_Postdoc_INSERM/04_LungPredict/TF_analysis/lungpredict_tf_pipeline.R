############################################################
# LungPredict TF analysis pipeline
# - builds full annotations with clusters
# - loads expression
# - computes TF activities (TCGA / GTEx)
# - generates expression & TF heatmaps
# - optional: compares TCGA vs GTEx, PCA
############################################################

suppressPackageStartupMessages({
  library(tidyverse)
  library(ComplexHeatmap)
  library(circlize)
  library(viper)
  library(dorothea)
  library(Hmisc)
})

############################################################
# 0. CONFIGURATION
############################################################

# Racine de ton projet (à adapter si besoin)
base_dir <- "A:/Downloads/Projects/workFromHome/Projects/LungPredict/TF_analysis"

# Fichiers d'entrée principaux
expression_file      <- file.path(base_dir, "LP_FFPE_STAR_RSEM_TPM.txt")
clinical_file        <- file.path(base_dir, "clinic_data_v2_clean.csv")
exp_clusters_file    <- file.path(base_dir, "expression_based_patient_clusters2.txt") # ou "a.txt"
reactome_clusters_file <- file.path(base_dir, "ReactomeClustersAllPatients.csv")
deconv_clusters_file <- file.path(base_dir, "DeconvCancerClusters.txt")

# Annotation fusionnée (output)
full_annotations_file <- file.path(base_dir, "full_annotations_with_clusters_corrected1.txt")

# Dossier de résultats
results_dir         <- file.path(base_dir, "results")
heatmap_dir         <- file.path(results_dir, "heatmaps")
tf_dir              <- file.path(results_dir, "tf_activities")
pca_dir             <- file.path(results_dir, "pca")

dir.create(results_dir, showWarnings = FALSE)
dir.create(heatmap_dir, showWarnings = FALSE)
dir.create(tf_dir, showWarnings = FALSE)
dir.create(pca_dir, showWarnings = FALSE)

############################################################
# 1. UTILITAIRES
############################################################

load_expression <- function(path) {
  read.delim(path, sep = "\t", check.names = FALSE) %>%
    column_to_rownames(var = "Gene") %>%
    as.matrix()
}

add_mutation_flags <- function(clinical_df) {
  egfr_mutant <- clinical_df$sample[clinical_df$mutation == "EGFR"]
  kras_mutant <- clinical_df$sample[grep("KRAS", clinical_df$mutation)]
  stk11_mutant <- clinical_df$sample[clinical_df$mutation == "STK11"]
  `%ni%` <- Negate(`%in%`)
  
  clinical_df %>%
    mutate(
      KRAS = case_when(
        sample %in% kras_mutant ~ "yes",
        sample %ni% kras_mutant ~ "no"
      ),
      EGFR = case_when(
        sample %in% egfr_mutant ~ "yes",
        sample %ni% egfr_mutant ~ "no"
      ),
      STK11 = case_when(
        sample %in% stk11_mutant ~ "yes",
        sample %ni% stk11_mutant ~ "no"
      )
    )
}

add_age_category <- function(clinical_df, age_col = "age") {
  clinical_df %>%
    mutate(age_category = case_when(
      .data[[age_col]] < 85 & .data[[age_col]] >= 80 ~ "80-85",
      .data[[age_col]] < 80 & .data[[age_col]] >= 75 ~ "75-80",
      .data[[age_col]] < 75 & .data[[age_col]] >= 70 ~ "70-75",
      .data[[age_col]] < 70 & .data[[age_col]] >= 65 ~ "65-70",
      .data[[age_col]] < 65 & .data[[age_col]] >= 60 ~ "60-65",
      .data[[age_col]] < 60 & .data[[age_col]] >= 55 ~ "55-60",
      .data[[age_col]] < 55 & .data[[age_col]] >= 50 ~ "50-55",
      .data[[age_col]] < 50 & .data[[age_col]] >= 45 ~ "45-50",
      .data[[age_col]] < 45 & .data[[age_col]] >= 40 ~ "40-45",
      TRUE ~ NA_character_
    ))
}

build_greyscale <- function(n = 10) {
  grey.colors(n, rev = TRUE)
}

############################################################
# 2. CONSTRUCTION DES ANNOTATIONS COMPLÈTES
############################################################

build_full_annotations <- function(clinical_file,
                                   exp_clusters_file,
                                   reactome_clusters_file,
                                   deconv_clusters_file,
                                   out_file) {
  
  clinical_features <- read.table(clinical_file, header = TRUE, sep = ",")
  
  clinical_features <- clinical_features %>%
    add_mutation_flags() %>%
    add_age_category(age_col = "age")
  
  # Début DF annotations de base
  annotations <- data.frame(
    sample      = clinical_features$sample,
    Stage       = clinical_features$Stage,
    Surgery     = clinical_features$Surgery,
    diagnostic  = clinical_features$diagnostic,
    Metastatic  = clinical_features$diagnoostic_VI_metastasis,
    Sexe        = clinical_features$sexe,
    Smoker      = clinical_features$smoker,
    Localisation= clinical_features$localisation_ponction,
    Age         = clinical_features$age_category,
    KRAS        = clinical_features$KRAS,
    EGFR        = clinical_features$EGFR,
    STK11       = clinical_features$STK11,
    stringsAsFactors = FALSE
  )
  
  rownames(annotations) <- annotations$sample
  
  # Filtrage patients (comme dans tes scripts)
  filtered_annotations <- annotations %>%
    filter(
      sample != "LP.01.56B",
      sample != "LP.01.69",
      Localisation != "extra-thoracic"
    )
  
  # Clusters expression
  expClusters <- read.table(exp_clusters_file, header = TRUE, sep = "\t")
  colnames(expClusters)[1] <- "sample"
  
  # Clusters reactome
  reactomeClusters <- read.table(reactome_clusters_file, header = TRUE, sep = ",")
  colnames(reactomeClusters)[1] <- "sample"
  
  # Clusters deconv
  deconvClusters <- read.table(deconv_clusters_file, header = TRUE, sep = "\t")
  colnames(deconvClusters)[1] <- "sample"
  
  # Jointures successives (équivalent à tes right_join)
  out <- right_join(expClusters, filtered_annotations, by = "sample")
  colnames(out)[2] <- "GE_Cluster"
  
  out <- right_join(reactomeClusters, out, by = "sample")
  colnames(out)[2] <- "reactomeCluster"
  
  out2 <- right_join(deconvClusters[, 1:2], out, by = "sample")
  out3 <- right_join(deconvClusters[, c(1, 3)], out2, by = "sample")
  
  # Ton script enlève la colonne 12 (ancienne colonne redondante)
  out4 <- out3[, -12]
  
  write.table(out4, file = out_file, sep = "\t", quote = FALSE, row.names = FALSE)
  
  message("Full annotations written to: ", out_file)
  invisible(out4)
}

############################################################
# 3. CHARGER EXPRESSION + ANNOTATIONS COMPLÈTES
############################################################

load_annotations_with_clusters <- function(full_annotations_file, expression_mat) {
  annotations <- read.table(full_annotations_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  
  rownames(annotations) <- annotations$sample
  
  # Aligner avec l’expression
  selector <- intersect(annotations$sample, colnames(expression_mat))
  expression_mat <- expression_mat[, selector]
  annotations <- annotations[selector, ]
  
  list(
    expression  = expression_mat,
    annotations = annotations
  )
}

build_ann_colors_full <- function(annotations) {
  greyscale <- build_greyscale(10)
  
  HeatmapAnnotation(
    Stage          = annotations$Stage,
    Surgery        = annotations$Surgery,
    Diagnostic     = annotations$diagnostic,
    Metastatic     = annotations$Metastatic,
    Sex            = annotations$Sexe,
    Smoker         = annotations$Smoker,
    Age            = annotations$Age,
    KRAS           = annotations$KRAS,
    EGFR           = annotations$EGFR,
    STK11          = annotations$STK11,
    Cancer_content = annotations$Cancer.clust,
    Deconv_clust   = annotations$Dclust,
    Reactome_clust = annotations$reactomeCluster,
    Gene_clust     = annotations$GE_Cluster,
    col = list(
      Stage = c(
        "IA-1" = greyscale[1], "IA-2" = greyscale[2], "IA-3" = greyscale[3],
        "IB"   = greyscale[4], "IIA"  = greyscale[5], "IIB" = greyscale[6],
        "IIIA" = greyscale[7], "IIIB" = greyscale[8], "IV-A" = greyscale[9],
        "IV-B" = greyscale[10]
      ),
      KRAS = c("yes" = "blue", "no" = "white"),
      EGFR = c("yes" = "red", "no" = "white"),
      STK11 = c("yes" = "green", "no" = "white"),
      Surgery = c("No" = "darkmagenta", "Yes" = "darkgoldenrod1"),
      Diagnostic = c(
        "Adenocarcinoma"       = "blue",
        "epidermoid carcinomas"= "pink",
        "Neuroendocrine tumors"= "black",
        "Other"                = "white"
      ),
      Metastatic = c("Yes" = "red", "No" = "white"),
      Sex = c("M" = "blue", "F" = "pink"),
      Smoker = c("No" = "black", "Yes" = "red", "Old" = "salmon", "Unknown" = "white"),
      Age = c(
        "40-45" = greyscale[1], "45-50" = greyscale[2], "50-55" = greyscale[3],
        "55-60" = greyscale[4], "60-65" = greyscale[5], "65-70" = greyscale[6],
        "70-75" = greyscale[7], "75-80" = greyscale[8], "80-85" = greyscale[9]
      ),
      Cancer_content = c("high" = "lightskyblue", "low" = "lightyellow"),
      Deconv_clust = c("c1" = "antiquewhite", "c2" = "antiquewhite4"),
      Reactome_clust = c("1" = "blue", "2" = "pink"),
      Gene_clust = c("1" = "chocolate1", "2" = "chocolate4", "3" = "chocolate3")
    )
  )
}

############################################################
# 4. HEATMAPS EXPRESSION
############################################################

plot_expression_heatmaps <- function(expression_mat, annotations, ann_colors, prefix = "expression") {
  # Heatmap gène x patient
  ht1 <- Heatmap(
    expression_mat,
    cluster_rows        = FALSE,
    cluster_row_slices  = FALSE,
    show_row_dend       = FALSE,
    top_annotation      = ann_colors,
    name                = "Expression",
    show_row_names      = FALSE,
    show_heatmap_legend = TRUE,
    row_names_gp        = gpar(fontsize = 7),
    column_names_gp     = gpar(fontsize = 7)
  )
  png(file.path(heatmap_dir, paste0(prefix, "_gene_by_patient.png")), width = 1600, height = 1200, res = 150)
  draw(ht1)
  dev.off()
  
  # Corr patient x patient (test de significativité via Hmisc::rcorr)
  cor_res <- rcorr(as.matrix(expression_mat))
  cor_mat <- cor_res$r
  cor_mat[cor_res$P > 0.05] <- 0
  
  colors <- colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))
  
  ht2 <- Heatmap(
    cor_mat, col = colors,
    top_annotation      = ann_colors,
    name                = "GE-PatientCorr",
    show_row_names      = TRUE,
    show_heatmap_legend = TRUE,
    row_names_gp        = gpar(fontsize = 7),
    column_names_gp     = gpar(fontsize = 7)
  )
  png(file.path(heatmap_dir, paste0(prefix, "_patient_correlation.png")), width = 1600, height = 1200, res = 150)
  draw(ht2)
  dev.off()
}

############################################################
# 5. ACTIVITÉS DE TF (TCGA & GTEX) + HEATMAPS
############################################################

compute_tf_activities_tcga_gtex <- function(expression_mat) {
  # TCGA
  data("dorothea_hs_pancancer", package = "dorothea")
  regulons_tcga <- dorothea_hs_pancancer %>%
    filter(confidence %in% c("A", "B"))
  
  tf_tcga <- run_viper(
    expression_mat,
    regulons_tcga,
    options = list(
      method      = "scale",
      minsize     = 4,
      eset.filter = FALSE,
      cores       = 1,
      verbose     = FALSE
    )
  )
  
  # GTEX
  data("dorothea_hs", package = "dorothea")
  regulons_gtex <- dorothea_hs %>%
    filter(confidence %in% c("A", "B"))
  
  tf_gtex <- run_viper(
    expression_mat,
    regulons_gtex,
    options = list(
      method      = "scale",
      minsize     = 4,
      eset.filter = FALSE,
      cores       = 1,
      verbose     = FALSE
    )
  )
  
  list(tcga = tf_tcga, gtex = tf_gtex)
}

plot_tf_heatmaps <- function(tf_mat, annotations, ann_colors, prefix = "TF_tcga", sd_threshold = 2) {
  # Corr patient x patient (TF)
  cor_pat <- cor(tf_mat)
  ht_pat <- Heatmap(
    cor_pat,
    top_annotation      = ann_colors,
    name                = paste0(prefix, "_patientCorr"),
    show_row_names      = TRUE,
    show_heatmap_legend = TRUE,
    row_names_gp        = gpar(fontsize = 7),
    column_names_gp     = gpar(fontsize = 7)
  )
  png(file.path(heatmap_dir, paste0(prefix, "_patientCorr.png")), width = 1600, height = 1200, res = 150)
  draw(ht_pat)
  dev.off()
  
  # Heatmap TF x patient (non filtré)
  ht_tf_full <- Heatmap(
    tf_mat,
    top_annotation      = ann_colors,
    name                = paste0(prefix, "_TFxPatient"),
    show_row_names      = TRUE,
    show_heatmap_legend = TRUE,
    row_names_gp        = gpar(fontsize = 7),
    column_names_gp     = gpar(fontsize = 7)
  )
  png(file.path(heatmap_dir, paste0(prefix, "_TFxPatient_full.png")), width = 1600, height = 1200, res = 150)
  draw(ht_tf_full)
  dev.off()
  
  # Filtrer TF sur variance
  TFstd <- apply(tf_mat, 1, sd)
  tf_filt <- tf_mat[TFstd > sd_threshold, , drop = FALSE]
  
  ht_tf_filt <- Heatmap(
    tf_filt,
    top_annotation      = ann_colors,
    name                = paste0(prefix, "_TFxPatient_filtered"),
    show_row_names      = TRUE,
    show_heatmap_legend = TRUE,
    row_names_gp        = gpar(fontsize = 7),
    column_names_gp     = gpar(fontsize = 7)
  )
  png(file.path(heatmap_dir, paste0(prefix, "_TFxPatient_filtered.png")), width = 1600, height = 1200, res = 150)
  draw(ht_tf_filt)
  dev.off()
  
  invisible(tf_filt)
}

compare_tcga_gtex <- function(tf_tcga, tf_gtex, prefix = "TCGA_GTEX") {
  tf_tcga_df <- as.data.frame(tf_tcga)
  tf_gtex_df <- as.data.frame(tf_gtex)
  
  # Filtrer pour n’avoir que les TF communs
  common_tfs <- intersect(rownames(tf_tcga_df), rownames(tf_gtex_df))
  tf_tcga_df <- tf_tcga_df[common_tfs, ]
  tf_gtex_df <- tf_gtex_df[common_tfs, ]
  
  # Corr patient x patient
  tcga_gtex_patientcor <- cor(tf_tcga_df, tf_gtex_df)
  png(file.path(heatmap_dir, paste0(prefix, "_patient_correlation.png")), width = 1600, height = 1200, res = 150)
  heatmap(tcga_gtex_patientcor)
  dev.off()
  
  # Corr TF x TF
  tcga_gtex_tfcor <- cor(t(tf_tcga_df), t(tf_gtex_df))
  png(file.path(heatmap_dir, paste0(prefix, "_TF_correlation.png")), width = 1600, height = 1200, res = 150)
  heatmap(tcga_gtex_tfcor)
  dev.off()
  
  # Différences par TF/patient (sauf SOX2 comme dans ton script)
  diff <- tf_tcga_df - tf_gtex_df
  diff2 <- diff[rownames(diff) != "SOX2", ] %>% as.matrix()
  ht_diff <- Heatmap(
    diff2,
    name            = "TCGA-GTEX",
    row_names_gp    = gpar(fontsize = 7),
    column_names_gp = gpar(fontsize = 7)
  )
  png(file.path(heatmap_dir, paste0(prefix, "_diff_TCGA_minus_GTEX.png")), width = 1600, height = 1200, res = 150)
  draw(ht_diff)
  dev.off()
}

############################################################
# 6. PCA SUR L’EXPRESSION
############################################################

run_pca_expression <- function(expression_mat, prefix = "expression") {
  res.pca <- prcomp(t(expression_mat), scale. = TRUE) # PCA sur les patients
  
  # Variance expliquée
  sink(file.path(pca_dir, paste0(prefix, "_pca_variance.txt")))
  print(summary(res.pca))
  sink()
  
  # Tu peux ajouter ici fviz_eig / fviz_pca_ind si tu veux, mais ça demande factoextra
  invisible(res.pca)
}

############################################################
# 7. MAIN : PIPELINE HAUT NIVEAU
############################################################

# 7.1 (optionnel) reconstruire full_annotations_with_clusters_corrected1.txt
# build_full_annotations(
#   clinical_file        = clinical_file,
#   exp_clusters_file    = exp_clusters_file,
#   reactome_clusters_file = reactome_clusters_file,
#   deconv_clusters_file = deconv_clusters_file,
#   out_file             = full_annotations_file
# )

# 7.2 Charger expression + annotations complètes
expression <- load_expression(expression_file)

loaded <- load_annotations_with_clusters(full_annotations_file, expression)
expression <- loaded$expression
annotations <- loaded$annotations
ann_colors  <- build_ann_colors_full(annotations)

# 7.3 Heatmaps basés sur l’expression
plot_expression_heatmaps(expression, annotations, ann_colors, prefix = "LP_expression")

# 7.4 Activités de TF (TCGA + GTEx)
tf_list <- compute_tf_activities_tcga_gtex(expression)
tf_tcga <- tf_list$tcga
tf_gtex <- tf_list$gtex

# Sauvegarder les matrices si tu veux
write.csv(as.data.frame(tf_tcga), file.path(tf_dir, "TF_activities_TCGA.csv"))
write.csv(as.data.frame(tf_gtex), file.path(tf_dir, "TF_activities_GTEX.csv"))

# 7.5 Heatmaps TF (TCGA)
tf_tcga_filtered <- plot_tf_heatmaps(tf_tcga, annotations, ann_colors, prefix = "TF_TCGA", sd_threshold = 2)

# 7.6 Comparaison TCGA vs GTEx
compare_tcga_gtex(tf_tcga, tf_gtex, prefix = "TCGA_GTEX")

# 7.7 PCA sur l’expression
run_pca_expression(expression, prefix = "LP")
