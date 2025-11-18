
#1 fetch the gene lists from file (mine was formatted in an atypical way, so i had to manipulate it a bit)
# check how the data looks like, especially for step #10 dfNinaMapped is how it should be formatted
gcSampleNinaRaw <- read.table("genes_list.txt", header = FALSE, sep = "\t", col.names = paste0("V",seq_len(max(count.fields("genes_list.txt", sep = '\t')))),, dec = ".", fill = TRUE, stringsAsFactors=FALSE)

#2 input the virus long to short names mapping file here
virus_names_mapping <- read.table("virus-host_PPI_info_nina.txt", header = FALSE, sep = '\t')
short_names = virus_names_mapping$V2
names(short_names) = virus_names_mapping$V1
#View(short_names)

#3 create a function that converts long names into short names
virus_name_mapping_function <- function(virus_name){short_names[virus_name]}


#4 map the long names to short names
gcSampleNina <- gcSampleNinaRaw
gcSampleNina$V1 <- sapply(gcSampleNinaRaw$V1, virus_name_mapping_function)
View(gcSampleNina$V1)

#5 transpose
dfNina <- data.frame(t(gcSampleNina[-1]))
colnames(dfNina) <- gcSampleNina[, 1]

#6 make mapping file and remove NA from it
mart_export_mappings <- read.table("mart_export.txt", "\t", header=TRUE, stringsAsFactors=FALSE)
mart_export_mappings_curated <- mart_export_mappings[!is.na(mart_export_mappings$NCBI.gene..formerly.Entrezgene..ID),]

#7 make a dict of the mapping file
entrez_names = mart_export_mappings_curated$NCBI.gene..formerly.Entrezgene..ID
names(entrez_names) = mart_export_mappings_curated$Gene.name

#8 create a function that converts gene symbols into entrez ids
my_function <- function(gene_name){entrez_names[gene_name]}

#9 map our gene symbols in our dataframe to entrez ids
dfNinaMapped <- as.data.frame(lapply(dfNina, my_function))

#10 run enrichment analysis and compare clusters together with clusterProfiler and function compareCluster
library(clusterProfiler)
library("ReactomePA")
library("reactome.db")

## pval 0.005
# ck <- compareCluster(geneCluster = dfNinaMapped, fun = "enrichGO", OrgDb='org.Hs.eg.db', pvalueCutoff = 0.005, ont="BP")
# ck <- compareCluster(geneCluster = dfNinaMapped, fun = "enrichGO", OrgDb='org.Hs.eg.db', pvalueCutoff = 0.005, ont="MF")
ck <- compareCluster(geneCluster = dfNinaMapped, fun = "enrichPathway", pvalueCutoff=0.005)

#11 plot results (optional)
library(ggplot2)
dotplot(ck, font.size=3)+ theme(axis.text.x = element_text(angle = 45, hjust = 1))

#12 transform the result as a dataframe
ck_df <- as.data.frame(ck, stringsAsFactors=FALSE)

#13 save results in a text file
library(readr)
write_tsv(ck_df, "enrichPathway-pvalueCutoff0.005.tsv")