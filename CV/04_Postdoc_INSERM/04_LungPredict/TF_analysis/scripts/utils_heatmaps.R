# scripts/utils_heatmaps.R

suppressPackageStartupMessages({
  library(ComplexHeatmap)
  library(circlize)
  library(Hmisc)   # pour rcorr
  library(grid)
})

# Heatmap expression gene x patient
plot_gene_expression_heatmap <- function(expression,
                                         ann,
                                         file = NULL) {

  if (!is.null(file)) png(file, width = 2000, height = 1600, res = 200)

  ht <- Heatmap(
    expression,
    cluster_rows        = FALSE,
    cluster_row_slices  = FALSE,
    show_row_dend       = FALSE,
    top_annotation      = ann,
    name                = "Expression",
    show_row_names      = FALSE,
    show_column_names   = TRUE,
    row_names_gp        = gpar(fontsize = 7),
    column_names_gp     = gpar(fontsize = 7)
  )

  draw(ht)
  if (!is.null(file)) dev.off()
}

# Corrélation patients basée sur l'expression (rcorr + p-value)
plot_expression_patient_corr <- function(expression,
                                         ann,
                                         file = NULL) {

  if (!is.null(file)) png(file, width = 2000, height = 1800, res = 200)

  cor_rcorr <- rcorr(as.matrix(expression))
  r_mat     <- cor_rcorr$r
  r_mat[cor_rcorr$P > 0.05] <- 0

  col_fun <- colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))

  ht <- Heatmap(
    r_mat,
    col               = col_fun,
    top_annotation    = ann,
    name              = "GE-PatientCorr",
    show_row_names    = TRUE,
    show_column_names = TRUE,
    row_names_gp      = gpar(fontsize = 7),
    column_names_gp   = gpar(fontsize = 7)
  )

  draw(ht)
  if (!is.null(file)) dev.off()
}

# Heatmap TF x patient
plot_tf_heatmap <- function(tf_matrix,
                            ann,
                            title = "TF activities",
                            file  = NULL) {

  if (!is.null(file)) png(file, width = 2000, height = 1800, res = 200)

  ht <- Heatmap(
    tf_matrix,
    top_annotation    = ann,
    name              = title,
    show_row_names    = TRUE,
    show_column_names = TRUE,
    row_names_gp      = gpar(fontsize = 7),
    column_names_gp   = gpar(fontsize = 7)
  )

  draw(ht)
  if (!is.null(file)) dev.off()
}

# Corrélation patient basée sur TF activities
plot_tf_patient_corr <- function(tf_matrix,
                                 ann,
                                 file = NULL) {

  if (!is.null(file)) png(file, width = 2000, height = 1800, res = 200)

  cor_mat <- cor(tf_matrix)

  ht <- Heatmap(
    cor_mat,
    top_annotation    = ann,
    name              = "TF-PatientCorr",
    show_row_names    = TRUE,
    show_column_names = TRUE,
    row_names_gp      = gpar(fontsize = 7),
    column_names_gp   = gpar(fontsize = 7)
  )

  draw(ht)
  if (!is.null(file)) dev.off()
}

# Corrélation TF x TF
plot_tf_tf_corr <- function(tf_matrix,
                            file = NULL) {

  if (!is.null(file)) png(file, width = 2000, height = 1800, res = 200)

  cor_mat <- cor(t(tf_matrix))

  ht <- Heatmap(
    cor_mat,
    name              = "TF-TF corr",
    show_row_names    = TRUE,
    show_column_names = TRUE,
    row_names_gp      = gpar(fontsize = 7),
    column_names_gp   = gpar(fontsize = 7)
  )

  draw(ht)
  if (!is.null(file)) dev.off()
}
