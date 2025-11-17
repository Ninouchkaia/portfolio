# scripts/utils_tf_activity.R

suppressPackageStartupMessages({
  library(tidyverse)
  library(viper)
  library(dorothea)
})

# Charge les régulons Dorothea (TCGA ou GTEx)
get_dorothea_regulons <- function(source = c("tcga", "gtex"),
                                  confidence = c("A", "B")) {

  source <- match.arg(source)

  if (source == "tcga") {
    data("dorothea_hs_pancancer", package = "dorothea")
    regulons_raw <- dorothea_hs_pancancer
  } else {
    data("dorothea_hs", package = "dorothea")
    regulons_raw <- dorothea_hs
  }

  regulons_raw %>%
    dplyr::filter(confidence %in% confidence)
}

# Calcule les activités de TF avec viper
compute_tf_activity <- function(expression,
                                source = c("tcga", "gtex"),
                                confidence = c("A", "B"),
                                minsize = 4,
                                cores   = 1) {

  regulons_tbl <- get_dorothea_regulons(source, confidence)

  tf_activities <- run_viper(
    expression,
    regulons_tbl,
    options = list(
      method      = "scale",
      minsize     = minsize,
      eset.filter = FALSE,
      cores       = cores,
      verbose     = FALSE
    )
  )

  as.matrix(tf_activities)
}

# Filtre les TFs peu variables
filter_tf_by_sd <- function(tf_matrix, sd_threshold = 2) {
  sds <- apply(tf_matrix, 1, sd)
  tf_matrix[sds > sd_threshold, , drop = FALSE]
}
