# scripts/utils_annotations.R

suppressPackageStartupMessages({
  library(tidyverse)
  library(ComplexHeatmap)
})

#-------------------------
# Chargement des données
#-------------------------

load_expression <- function(expr_path) {
  expr <- read.delim(expr_path, check.names = FALSE)
  if (!"Gene" %in% colnames(expr)) {
    stop("La colonne 'Gene' est absente de ", expr_path)
  }
  expr %>%
    column_to_rownames("Gene") %>%
    as.matrix()
}

load_full_annotations <- function(annot_path) {
  read.delim(annot_path, check.names = FALSE)
}

# Aligne expression et annotations sur les mêmes patients
align_expression_annotations <- function(expression,
                                         annotations,
                                         sample_col = "sample") {

  if (!sample_col %in% colnames(annotations)) {
    stop("Colonne d'échantillon '", sample_col,
         "' absente des annotations.")
  }

  selector <- intersect(annotations[[sample_col]], colnames(expression))
  if (length(selector) == 0) {
    stop("Aucun échantillon commun entre expression et annotations.")
  }

  expr_f <- expression[, selector, drop = FALSE]
  ann_f  <- annotations[match(selector, annotations[[sample_col]]),
                        , drop = FALSE]
  rownames(ann_f) <- ann_f[[sample_col]]

  list(
    expression  = expr_f,
    annotations = ann_f
  )
}

#-------------------------
# Construction HeatmapAnnotation
#-------------------------

build_sample_annotation <- function(annotations) {

  greyscale <- grey.colors(10, rev = TRUE)

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
        "IB"   = greyscale[4], "IIA"  = greyscale[5], "IIB"  = greyscale[6],
        "IIIA" = greyscale[7], "IIIB" = greyscale[8],
        "IV-A" = greyscale[9], "IV-B" = greyscale[10]
      ),
      Surgery = c("No" = "darkmagenta", "Yes" = "darkgoldenrod1"),
      Diagnostic = c(
        "Adenocarcinoma"       = "blue",
        "epidermoid carcinomas"= "pink",
        "Neuroendocrine tumors"= "black",
        "Other"                = "white"
      ),
      Metastatic = c("Yes" = "red", "No" = "white"),
      Sex   = c("M" = "blue", "F" = "pink"),
      Smoker = c("No" = "black", "Yes" = "red",
                 "Old" = "salmon", "Unknown" = "white"),
      Age = c(
        "40-45" = greyscale[1], "45-50" = greyscale[2],
        "50-55" = greyscale[3], "55-60" = greyscale[4],
        "60-65" = greyscale[5], "65-70" = greyscale[6],
        "70-75" = greyscale[7], "75-80" = greyscale[8],
        "80-85" = greyscale[9]
      ),
      KRAS = c("yes" = "blue",  "no" = "white"),
      EGFR = c("yes" = "red",   "no" = "white"),
      STK11= c("yes" = "green", "no" = "white"),
      Cancer_content = c("high" = "lightskyblue", "low" = "lightyellow"),
      Deconv_clust = c("c1" = "antiquewhite", "c2" = "antiquewhite4"),
      Reactome_clust = c("1" = "blue", "2" = "pink"),
      Gene_clust = c("1" = "chocolate1", "2" = "chocolate4", "3" = "chocolate3")
    )
  )
}
