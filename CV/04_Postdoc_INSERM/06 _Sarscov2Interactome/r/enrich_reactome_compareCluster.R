#!/usr/bin/env Rscript

# Usage: Rscript enrich_reactome_compareCluster.R genes_list_file out_tsv

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: enrich_reactome_compareCluster.R <genes_list.txt> <out.tsv>")
}

genes_list_file <- args[1]
out_tsv <- args[2]

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(ReactomePA)
  library(reactome.db)
  library(readr)
})

# 1. Lire genes_list.txt (même format que dans tes scripts)
gcRaw <- read.table(
  genes_list_file,
  header = FALSE,
  sep = "\t",
  col.names = paste0("V", seq_len(max(count.fields(genes_list_file, sep = "\t")))),
  dec = ".",
  fill = TRUE,
  stringsAsFactors = FALSE
)

# 2. Mapping long→short names (fichier à garder dans le même dossier que ce script)
virus_names_mapping <- read.table("virus-host_PPI_info_nina.txt",
                                  header = FALSE,
                                  sep = "\t",
                                  stringsAsFactors = FALSE)
short_names <- virus_names_mapping$V2
names(short_names) <- virus_names_mapping$V1
virus_name_mapping_function <- function(virus_name) short_names[virus_name]

gc <- gcRaw
gc$V1 <- sapply(gcRaw$V1, virus_name_mapping_function)

# 3. Transposition pour obtenir une liste nommée de vecteurs de gènes
df <- data.frame(t(gc[-1]))
colnames(df) <- gc[, 1]

# 4. Mapping Gene symbol → Entrez
mart_export <- read.table("mart_export.txt",
                          sep = "\t",
                          header = TRUE,
                          stringsAsFactors = FALSE)
mart_export_curated <- mart_export[!is.na(mart_export$NCBI.gene..formerly.Entrezgene..ID), ]

entrez_names <- mart_export_curated$NCBI.gene..formerly.Entrezgene..ID
names(entrez_names) <- mart_export_curated$Gene.name

to_entrez <- function(gene_name) entrez_names[gene_name]
df_mapped <- as.data.frame(lapply(df, to_entrez))

# 5. Enrichissement Reactome
ck <- compareCluster(
  geneCluster = df_mapped,
  fun = "enrichPathway",
  pvalueCutoff = 0.005
)

# 6. Export TSV
ck_df <- as.data.frame(ck, stringsAsFactors = FALSE)
write_tsv(ck_df, out_tsv)
