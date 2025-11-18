

library(stats)

a <- phyper(4, 247, 18423, 4, lower.tail = TRUE)
?phyper

# fetch the gene lists from file
gcSampleNinaRaw <- read.table("genes_list.txt", header = FALSE, sep = "\t", col.names = paste0("V",seq_len(max(count.fields("genes_list.txt", sep = '\t')))),, dec = ".", fill = TRUE, stringsAsFactors=FALSE)
# gcSampleNinaRaw <- read.table("first_and_second_order_genes_list.txt", header = FALSE, sep = "\t", col.names = paste0("V",seq_len(max(count.fields("first_and_second_order_genes_list.txt", sep = '\t')))),, dec = ".", fill = TRUE, stringsAsFactors=FALSE)
# gcSampleNinaRaw <- read.table("first_and_second_and_third_order_genes_list.txt", header = FALSE, sep = "\t", col.names = paste0("V",seq_len(max(count.fields("first_and_second_and_third_order_genes_list.txt", sep = '\t')))),, dec = ".", fill = TRUE, stringsAsFactors=FALSE)
# gcSampleNina <- read.table("first_and_secondThirdFourth_order_genes_list.txt", header = FALSE, sep = "\t", col.names = paste0("V",seq_len(max(count.fields("first_and_secondThirdFourth_order_genes_list.txt", sep = '\t')))),, dec = ".", fill = TRUE, stringsAsFactors=FALSE)

# gcSampleNina <- read.table("only_second_order_genes_list.txt", header = FALSE, sep = "\t", col.names = paste0("V",seq_len(max(count.fields("only_second_order_genes_list.txt", sep = '\t')))),, dec = ".", fill = TRUE, stringsAsFactors=FALSE)
# gcSampleNina <- read.table("only_Third_order_genes_list.txt", header = FALSE, sep = "\t", col.names = paste0("V",seq_len(max(count.fields("only_Third_order_genes_list.txt", sep = '\t')))),, dec = ".", fill = TRUE, stringsAsFactors=FALSE)
# gcSampleNina <- read.table("only_Fourth_order_genes_list.txt", header = FALSE, sep = "\t", col.names = paste0("V",seq_len(max(count.fields("only_Fourth_order_genes_list.txt", sep = '\t')))),, dec = ".", fill = TRUE, stringsAsFactors=FALSE)
# gcSampleNina <- read.table("only_Fifth_order_genes_list.txt", header = FALSE, sep = "\t", col.names = paste0("V",seq_len(max(count.fields("only_Fifth_order_genes_list.txt", sep = '\t')))),, dec = ".", fill = TRUE, stringsAsFactors=FALSE)


virus_names_mapping <- read.table("virus-host_PPI_info_nina.txt", header = FALSE, sep = '\t')

short_names = virus_names_mapping$V2
names(short_names) = virus_names_mapping$V1

# print(short_names["Chikungunya virus (strain S27-African prototype)"])
View(short_names)

# create a function that converts long names into short names
virus_name_mapping_function <- function(virus_name){short_names[virus_name]}


# map the long names to short names
gcSampleNina <- gcSampleNinaRaw
gcSampleNina$V1 <- sapply(gcSampleNinaRaw$V1, virus_name_mapping_function) #better to use sapply instead of lapply because we're working with a single column ie a vector
View(gcSampleNina$V1)

#transpose
dfNina <- data.frame(t(gcSampleNina[-1]))
colnames(dfNina) <- gcSampleNina[, 1]
# typeof(dfNina$Swinepox[2])

# sum <- sapply(dfNina, function(x) sum(x!=''))
# sum_values <- data.frame(sum)
# barplot(sum_values$sum, main='number of human protein interactors', xlab = names(sum_values))
# # print(colnames(sum_values))
# 
# sum_values$virus <- rownames(sum_values)
# 
# par(mar=c(2, 15, 3, 10)) # 15 line height for bottom margin
# barplot(sum_values$sum, main='number of human protein interactors', names = sum_values$virus, horiz=T , las=1, font.axis=2, font.lab=1, cex.names=0.5
#         )


# fetch mapping entrez - gene symbols (it will take the first mappings if there are synonyms) # check the mappings !!!
mart_export_mappings <- read.table("mart_export.txt", "\t", header=TRUE, stringsAsFactors=FALSE)
# mart_export_mappings <- read.table("mart_export_entrez_uniprot.txt", "\t", header=TRUE, stringsAsFactors=FALSE)

#remove NA from the mapping file
# mart_export_mappings_curated <- mart_export_mappings[!is.na(mart_export_mappings$NCBI.gene..formerly.Entrezgene..ID),]
mart_export_mappings_curated <- mart_export_mappings[!is.na(mart_export_mappings$NCBI.gene..formerly.Entrezgene..ID),]


# make a dict of the mapping file
entrez_names = mart_export_mappings_curated$NCBI.gene..formerly.Entrezgene..ID
names(entrez_names) = mart_export_mappings_curated$Gene.name # check !
# entrez_names = mart_export_mappings_curated$NCBI.gene..formerly.Entrezgene..ID
# names(entrez_names) = mart_export_mappings_curated$UniProtKB.Gene.Name.symbol # check !

print(entrez_names["PRKACA"])


# create a function that converts gene symbols into entrez ids
my_function <- function(gene_name){entrez_names[gene_name]}

# map our gene symbols in our dataframe to entrez ids
dfNinaMapped <- as.data.frame(lapply(dfNina, my_function))

# run enrichment analysis and compare clusters together with clusterProfiler and function compareCluster
library(clusterProfiler)


# ?compareCluster
# ??enrichResult
# ck <- compareCluster(geneCluster = dfNinaMapped, fun = "enrichKEGG") #, pvalueCutoff=0.05)
# ck <- compareCluster(geneCluster = dfNinaMapped, fun = "enrichKEGG", pvalueCutoff=0.1)

# viewKEGG(ck)

# BiocManager::install("ReactomePA")
# BiocManager::install("reactome.db")
library("ReactomePA")
library("reactome.db")

## pval 0.005
# ck <- compareCluster(geneCluster = dfNinaMapped, fun = "enrichGO", OrgDb='org.Hs.eg.db', pvalueCutoff = 0.005, ont="BP")
# ck <- compareCluster(geneCluster = dfNinaMapped, fun = "enrichGO", OrgDb='org.Hs.eg.db', pvalueCutoff = 0.005, ont="MF")
# ck <- compareCluster(geneCluster = dfNinaMapped, fun = "enrichPathway", pvalueCutoff=0.005)


ck <- compareCluster(geneCluster = dfNinaMapped, fun = "enrichPathway", pvalueCutoff=0.05)
# ck <- compareCluster(geneCluster = dfNinaMapped, fun = "enrichPathway", pvalueCutoff=0.1)


# ck <- compareCluster(geneCluster = dfNinaMapped, fun = "enrichGO", OrgDb='org.Hs.eg.db', pvalueCutoff = 0.05, ont="BP")
# ck <- compareCluster(geneCluster = dfNinaMapped, fun = "enrichGO", OrgDb='org.Hs.eg.db', pvalueCutoff = 0.1, ont="BP")
# ck <- compareCluster(geneCluster = dfNinaMapped, fun = "enrichGO", OrgDb='org.Hs.eg.db', pvalueCutoff = 0.1, ont="MF")
# ck <- compareCluster(geneCluster = dfNinaMapped, fun = "enrichGO", OrgDb='org.Hs.eg.db', pvalueCutoff = 0.05, ont="MF")




# BiocManager::install("DOSE")
# library("DOSE")
# ck <- compareCluster(geneCluster = dfNinaMapped, fun = "enrichDO")

# ?compareCluster # which fun?

library(ggplot2)
dotplot(ck, font.size=3)+ theme(axis.text.x = element_text(angle = 45, hjust = 1))
  #theme(axis.text.x = element_text(angle = 90))#, split = TRUE) #,showCategory = 50)

# ?dotplot

# transform the result as a dataframe
ck_df <- as.data.frame(ck, stringsAsFactors=FALSE)

# install.packages("readr")

library(readr)
# write_tsv(ck_df, "enrichPathway-3rdorder-pvalueCutoff0.005.tsv")
# write_tsv(ck_df, "enrichGOBP-3dorder-pvalueCutoff0.005.tsv")
# write_tsv(ck_df, "enrichGOMF-3dorder-pvalueCutoff0.005.tsv")

# write_tsv(ck_df, "enrichPathway-1storder-pvalueCutoff0.05.tsv")



#enrichGO_MF_pval0005_Order2Only.tsv")
# ?write_tsv

# extract the data from the compareResult function and make a matrix to draw a heatmap
descriptionList <- unique(ck_df$Description)
clusterList <- unique(ck_df$Cluster)
# View(descriptionList)
# init df with O values

myClusterMap <- matrix(0,nrow=length(descriptionList),ncol = length(clusterList))
rownames(myClusterMap)<-descriptionList
colnames(myClusterMap) <- clusterList





for ( i in 1:nrow(ck_df)) {
  print(i)
  description <- ck_df$Description[i]
  print(description)
  cluster <- ck_df$Cluster[i]
  print(cluster)
  myClusterMap[rownames(myClusterMap) == description, colnames(myClusterMap) == cluster] <- 1  ## trop fort Miguel (problem with the dots in the column names)
  
}

# View(myClusterMap)
# View(ck_df)

library(pheatmap)

library(viridis)

cor.matrix <- cor(myClusterMap)
pheatmap(cor.matrix, fontsize_row = 5, fontsize_column = 1)


pheatmap(myClusterMap, fontsize_row = 1, cluster_rows = FALSE, clustering_distance_cols = "binary", clustering_method = "complete", color=viridis(100))
# pheatmap(myClusterMap, fontsize_row = 1, cluster_rows = FALSE, clustering_distance_cols = "manhattan", clustering_method = "ward.D2")

?pheatmap

my_clusters <- pheatmap(myClusterMap)
my_clusters
sorted_clusters_4 <- sort(cutree(my_clusters$tree_col, k=4))
View(sorted_clusters_4)

sorted_clusters_8 <- sort(cutree(my_clusters$tree_col, k=8))

View(sorted_clusters_8)


sorted_clusters_10 <- sort(cutree(my_clusters$tree_col, k=10))
View(sorted_clusters_10)

sorted_clusters_20 <- sort(cutree(my_clusters$tree_col, k=20))
View(sorted_clusters_20)

plot(my_clusters$tree_col)
abline(h=4, col="red", lty=2, lwd=2)
?abline



# I need to get my input from python processing 
# as virusA -> gene A
# virusB -> gene B
# virusB -> gene A
# etc

# make a list out of all the viruses and genes that are viral interactors
# myGeneVirusTable <- read.table("virus_gene_table.txt", header=FALSE, sep="\t")
myGeneVirusTable <- read.table("virus_gene_table_third_order.txt", header=FALSE, sep="\t")
myGeneVirusTable$V1 <- sapply(myGeneVirusTable$V1, virus_name_mapping_function) 
View(myGeneVirusTable)

myVirusList <- unique(myGeneVirusTable$V1)
myGenesList <- unique(myGeneVirusTable$V2)
View(myGenesList)

# make empty matrix of size of genes_list
myGenesMatrix <- matrix(0,nrow=length(myGenesList), ncol = length(myVirusList))
rownames(myGenesMatrix)<- myGenesList
colnames(myGenesMatrix) <- myVirusList

#myVirusList[names(myVirusList) %in% "SARS-CoV-2" == FALSE] 

# fill it with 1 if gene is among viral interactors of specific virus
for ( i in 1:nrow(myGeneVirusTable)) {
  gene <- myGeneVirusTable$V2[i]
  virus <- myGeneVirusTable$V1[i]
  myGenesMatrix[rownames(myGenesMatrix) == gene, colnames(myGenesMatrix) == virus] <- 1 
  
}

install.packages("viridis")

library(viridis)

library(pheatmap)
pheatmap(myGenesMatrix, fontsize_row = 1, fonsize_col = 2, treeheight_col = 100, cluster_rows = FALSE, clustering_distance_cols = "binary", clustering_method = "complete", color=viridis(100))
?pheatmap

clusters4 <- sort(cutree(my_clusters$tree_col, k=4))
View(clusters4)

clusters8 <- sort(cutree(my_clusters$tree_col, k=8))
View(clusters8)







  # pheatmap(gcSampleNina)


BiocManager::install("ComplexHeatmap")
library("ComplexHeatmap")
?Heatmap
