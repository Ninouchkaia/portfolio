
Gene expression and TF activities plots (Nina).R 
annotations <- read.table("A:/Downloads/Projects/workFromHome/Projects/LungPredict/TF_analysis/full_annotations_with_clusters_corrected1.txt", sep="	", header=T)


#loading expression file
expression=read.csv("A:/Downloads/Projects/workFromHome/Projects/LungPredict/TF_analysis/LP_FFPE_STAR_RSEM_TPM.txt", sep = "	")
expression=expression%>% column_to_rownames(var='Gene') %>% as.matrix()

selector <- intersect(annotations$sample, colnames(expression))

# filter expression for the annotated patients
expression  <- expression[,selector]

### keep only samples that are also in expression
rownames(annotations) <- annotations$sample

# filtered_annotations <- annotations[annotations$sample %in% colnames(expression),]
filtered_annotations <- annotations[selector,]
annotations <- filtered_annotations

# create the annotation
greyscale <- grey.colors(10, rev = T)

ann_colors <- HeatmapAnnotation(
  Stage = annotations$Stage, 
  Surgery = annotations$Surgery,
  Diagnostic = annotations$diagnostic,
  Metastatic = annotations$Metastatic,
  Sex = annotations$Sexe,
  Smoker = annotations$Smoker,
  Age = annotations$Age,
  KRAS = annotations$KRAS,
  EGFR = annotations$EGFR,
  STK11 = annotations$STK11,
  Cancer_content = annotations$Cancer.clust,
  Deconv_clust = annotations$Dclust,
  Reactome_clust = annotations$reactomeCluster,
  Gene_clust = annotations$GE_Cluster,
  
  col = list(
    Stage = c(
      "IA-1" = greyscale[1], "IA-2" = greyscale[2], "IA-3" = greyscale[3],
      "IB" = greyscale[4], "IIA" = greyscale[5], "IIB" = greyscale[6],
      "IIIA" = greyscale[7], "IIIB" = greyscale[8], "IV-A" = greyscale[9],
      "IV-B" = greyscale[10]
    ),
    KRAS = c("yes" = "blue", "no" = "white"),
    EGFR = c("yes" = "red", "no" = "white"),
    STK11 = c("yes" = "green", "no" = "white"),
    Surgery = c("No" = "darkmagenta", "Yes" = "darkgoldenrod1"),
    Diagnostic = c("Adenocarcinoma" = "blue", "epidermoid carcinomas" = "pink", "Neuroendocrine tumors" = "black", "Other" = "white"),
    Metastatic = c("Yes" = "red", "No" = "white"),
    Sex = c("M" = "blue", "F" = "pink"),
    Smoker = c("No" = "black", "Yes" = "red", "Old" = "salmon", "Unknown" = "white"),
    Age = c(
      "40-45" = greyscale[1], "45-50" = greyscale[2], "50-55" = greyscale[3],
      "55-60" = greyscale[4], "60-65" = greyscale[5], "65-70" = greyscale[6],
      "70-75" = greyscale[7], "75-80" = greyscale[8], "80-85" = greyscale[9]
    ),
    Cancer_content = c("high" = "lightskyblue", "low" = "lightyellow"),
    Deconv_clust = c("c1" = "antiquewhite", "c2" = "antiquewhite4" ),
    Reactome_clust = c( "1" = "blue", "2" = "pink") ,
    Gene_clust = c("1" = "chocolate1", "2" = "chocolate4", "3" = "chocolate3")
  )
)

# #plot the matrix gene x patient
# ?Heatmap
Heatmap(expression,
        cluster_rows = FALSE,
        cluster_row_slices = FALSE,
        show_row_dend = FALSE,
        top_annotation = ann_colors,
        name = "Patient Correlation based on Expression",
        show_row_names = F,
        show_heatmap_legend = T,
        row_names_gp = gpar(fontsize = 7),
        column_names_gp = gpar(fontsize =7),
)

# #plot the matrix patient x patient correlation based on gene expression

library(Hmisc)
cor.matrix.expression2 <- rcorr(as.matrix(expression))

cor.matrix.expression2$r[cor.matrix.expression2$P > 0.05] <- 0
colors <- circlize::colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))

Heatmap(cor.matrix.expression2$r, col = colors,
        top_annotation = ann_colors, 
        name = "GE-PatientCorr",
        show_row_names = T, 
        show_heatmap_legend = T, 
        row_names_gp = gpar(fontsize = 7),
        column_names_gp = gpar(fontsize =7),
)


############################################
################ TF ACTIVITY ###############
############################################

### TCGA reference for regulons
data("dorothea_hs_pancancer", package = "dorothea")
regulons_tcga = dorothea_hs_pancancer %>%
  filter(confidence %in% c("A", "B"))

tf_activities_tcga <- run_viper(expression, regulons_tcga, 
                                options =  list(method = "scale", minsize = 4,eset.filter = FALSE, cores = 1, verbose = FALSE))



## plot patient correlation based on TF activities
cor.matrix.TF.filtered <- cor(tf_activities_tcga)

Heatmap(cor.matrix.TF.filtered, 
        top_annotation = ann_colors, 
        name = "patient correlation based on TF activities",
        show_row_names = T, show_heatmap_legend = T, 
        row_names_gp = gpar(fontsize = 7),
        column_names_gp = gpar(fontsize =7),
)

## plot TF x patient
Heatmap(tf_activities_tcga, 
        top_annotation = ann_colors, 
        name = "patient x TF activities not filtered",
        show_row_names = T, show_heatmap_legend = T, 
        row_names_gp = gpar(fontsize = 7),
        column_names_gp = gpar(fontsize =7),
)

## keep TF that vary a lot only
TFstd <- apply(tf_activities_tcga, 1, sd)
# print(TFstd)
hist(TFstd, breaks=106)

### Remove all TF that have std < 2
my_matrix_TF_filtered <- tf_activities_tcga[TFstd > 2,]


Heatmap(my_matrix_TF_filtered, 
        top_annotation = ann_colors, 
        name = "patient x TF activities not filtered",
        show_row_names = T, show_heatmap_legend = T, 
        row_names_gp = gpar(fontsize = 7),
        column_names_gp = gpar(fontsize =7),
)

Heatmaps_annotations.R 
annotations <- read.table("A:/Downloads/Projects/workFromHome/Projects/LungPredict/TF_analysis/full_annotations_with_clusters_corrected1.txt", sep="	", header=T)



#loading expression file
expression=read.csv("A:/Downloads/Projects/workFromHome/Projects/LungPredict/TF_analysis/LP_FFPE_STAR_RSEM_TPM.txt", sep = "	")
expression=expression%>% column_to_rownames(var='Gene') %>% as.matrix()

selector <- intersect(annotations$sample, colnames(expression))

# filter expression for the annotated patients
expression  <- expression[,selector]

### keep only samples that are also in expression
rownames(annotations) <- annotations$sample

# filtered_annotations <- annotations[annotations$sample %in% colnames(expression),]
filtered_annotations <- annotations[selector,]
annotations <- filtered_annotations

# create the annotation
greyscale <- grey.colors(10, rev = T)

ann_colors <- HeatmapAnnotation(
  Stage = annotations$Stage, 
  Surgery = annotations$Surgery,
  Diagnostic = annotations$diagnostic,
  Metastatic = annotations$Metastatic,
  Sex = annotations$Sexe,
  Smoker = annotations$Smoker,
  Age = annotations$Age,
  KRAS = annotations$KRAS,
  EGFR = annotations$EGFR,
  STK11 = annotations$STK11,
  Cancer_content = annotations$Cancer.clust,
  Deconv_clust = annotations$Dclust,
  Reactome_clust = annotations$reactomeCluster,
  Gene_clust = annotations$GE_Cluster,
  
  col = list(
    Stage = c(
      "IA-1" = greyscale[1], "IA-2" = greyscale[2], "IA-3" = greyscale[3],
      "IB" = greyscale[4], "IIA" = greyscale[5], "IIB" = greyscale[6],
      "IIIA" = greyscale[7], "IIIB" = greyscale[8], "IV-A" = greyscale[9],
      "IV-B" = greyscale[10]
    ),
    KRAS = c("yes" = "blue", "no" = "white"),
    EGFR = c("yes" = "red", "no" = "white"),
    STK11 = c("yes" = "green", "no" = "white"),
    Surgery = c("No" = "darkmagenta", "Yes" = "darkgoldenrod1"),
    Diagnostic = c("Adenocarcinoma" = "blue", "epidermoid carcinomas" = "pink", "Neuroendocrine tumors" = "black", "Other" = "white"),
    Metastatic = c("Yes" = "red", "No" = "white"),
    Sex = c("M" = "blue", "F" = "pink"),
    Smoker = c("No" = "black", "Yes" = "red", "Old" = "salmon", "Unknown" = "white"),
    Age = c(
      "40-45" = greyscale[1], "45-50" = greyscale[2], "50-55" = greyscale[3],
      "55-60" = greyscale[4], "60-65" = greyscale[5], "65-70" = greyscale[6],
      "70-75" = greyscale[7], "75-80" = greyscale[8], "80-85" = greyscale[9]
    ),
    Cancer_content = c("high" = "lightskyblue", "low" = "lightyellow"),
    Deconv_clust = c("c1" = "antiquewhite", "c2" = "antiquewhite4" ),
    Reactome_clust = c( "1" = "blue", "2" = "pink") ,
    Gene_clust = c("1" = "chocolate1", "2" = "chocolate4", "3" = "chocolate3")
  )
)

# #plot the matrix gene x patient
# ?Heatmap
Heatmap(expression,
        cluster_rows = FALSE,
        cluster_row_slices = FALSE,
        show_row_dend = FALSE,
        top_annotation = ann_colors,
        name = "Patient Correlation based on Expression",
        show_row_names = F,
        show_heatmap_legend = T,
        row_names_gp = gpar(fontsize = 7),
        column_names_gp = gpar(fontsize =7),
)

# #plot the matrix patient x patient correlation based on gene expression
# cor.matrix.expression <- cor(expression)

library(Hmisc)
cor.matrix.expression2 <- rcorr(as.matrix(expression))

cor.matrix.expression2$r[cor.matrix.expression2$P > 0.05] <- 0
colors <- circlize::colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))

Heatmap(cor.matrix.expression2$r, col = colors,
        top_annotation = ann_colors, 
        name = "GE-PatientCorr",
        show_row_names = T, 
        show_heatmap_legend = T, 
        row_names_gp = gpar(fontsize = 7),
        column_names_gp = gpar(fontsize =7),
)



# my_matrix <- cor.matrix.expression
# 
# my_matrix_filtered <- my_matrix[selector , selector]
# 
# 
# 
# 
# # annotations <- annotations[,selector]
# # annotations <- annotations[,-1]
# # annotations <- subset(annotations, select = -c(sample) )
# 
# 
# 
# 
# 
# 
# Heatmap(my_matrix_filtered,
#         top_annotation = ann_colors, 
#         name = "Patient Correlation based on Expression",
#         show_row_names = T, 
#         show_heatmap_legend = T, 
#         row_names_gp = gpar(fontsize = 7),
#         column_names_gp = gpar(fontsize =7),
# )


############################################
################ TF ACTIVITY ###############
############################################

### TCGA reference for regulons
data("dorothea_hs_pancancer", package = "dorothea")
regulons_tcga = dorothea_hs_pancancer %>%
  filter(confidence %in% c("A", "B"))

tf_activities_tcga <- run_viper(expression, regulons_tcga, 
                                options =  list(method = "scale", minsize = 4,eset.filter = FALSE, cores = 1, verbose = FALSE))



## plot patient correlation based on TF activities
cor.matrix.TF.filtered <- cor(tf_activities_tcga)

Heatmap(cor.matrix.TF.filtered, 
        top_annotation = ann_colors, 
        name = "patient correlation based on TF activities",
        show_row_names = T, show_heatmap_legend = T, 
        row_names_gp = gpar(fontsize = 7),
        column_names_gp = gpar(fontsize =7),
)

## plot TF x patient
Heatmap(tf_activities_tcga, 
        top_annotation = ann_colors, 
        name = "patient x TF activities not filtered",
        show_row_names = T, show_heatmap_legend = T, 
        row_names_gp = gpar(fontsize = 7),
        column_names_gp = gpar(fontsize =7),
)

## keep TF that vary a lot only
TFstd <- apply(tf_activities_tcga, 1, sd)
# print(TFstd)
hist(TFstd, breaks=106)

### Remove all TF that have std < 2
my_matrix_TF_filtered <- tf_activities_tcga[TFstd > 2,]


Heatmap(my_matrix_TF_filtered, 
        top_annotation = ann_colors, 
        name = "patient x TF activities not filtered",
        show_row_names = T, show_heatmap_legend = T, 
        row_names_gp = gpar(fontsize = 7),
        column_names_gp = gpar(fontsize =7),
)










?Heatmap








# then your heatmap
Heatmap(whatever_data_you_use,
        top_annotation = ann_colors, # annotations just created
        name = "Name on the main legend",
        show_heatmap_legend = T, 
        row_names_gp = gpar(fontsize = 8), column_names_gp = gpar(fontsize = 8) # set size of row and columns names
)

TF_activities.R 
#required packages
# install.packages("pillar")
# install.packages("ggplot2")
# install.packages("dplyr")
# install.packages("forcats")
# BiocManager::install("viper")
#BiocManager::install("dorothea")
# ??dorothea
#install.packages("devtools")
# library(devtools)
# install_github("saezlab/dorothea") # to get the TCGA (pan cancer) ref for regulons

library("tidyverse")
library("viper")
library("dorothea")
library(pheatmap)
library(RColorBrewer)
library(viridis)
library(ComplexHeatmap)

#################################
######## EXPRESSION DATA ########
#################################

#loading expression file
expression=read.csv("A:/Downloads/Projects/workFromHome/Projects/LungPredict/TF_analysis/LP_FFPE_STAR_RSEM_TPM.txt", sep = "	")
expression=expression%>% column_to_rownames(var='Gene') %>% as.matrix()
print(colnames(expression))
### remove patient 69 from expression data
# expression <- subset(expression, select = -c(LP.01.69))

#########################################################
###### FILTER PATIENT BASED ON CLINICAL FEATURES  #######
#########################################################

### create annotations
clinical_features <- read.table("A:/Downloads/Projects/workFromHome/Projects/LungPredict/TF_analysis/clinic_data_v2_clean.csv", header = TRUE, sep = ",")
# add age category
clinical_features <- clinical_features %>%
  mutate(age_category = case_when(
    age < 85 & age >= 80 ~ '80-85',
    age < 80 & age >= 75 ~ '75-80',
    age < 75 & age >= 70 ~ '70-75',
    age < 70 & age >= 65 ~ '65-70',
    age < 65 & age >= 60 ~ '60-65',
    age < 60 & age >= 55 ~ '55-60',
    age < 55 & age >= 50 ~ '50-55',
    age < 50 & age >= 45 ~ '45-50',
    age < 45 & age >= 40 ~ '40-45'
  ))

annotations <- data.frame(clinical_features$sample, 
                          clinical_features$Stage, 
                          clinical_features$mutation, 
                          clinical_features$Surgery, 
                          clinical_features$diagnoostic_VI_metastasis, 
                          clinical_features$sexe, 
                          clinical_features$smoker, 
                          clinical_features$localisation_ponction,
                          clinical_features$age_category)
rownames(annotations) <- annotations$clinical_features.sample


colnames(annotations) <- c('sample', 'Stage', 'Mutation', 'Surgery', 'Metastatic', 'Sexe', 'Smoker', 'Localisation', 'Age')


### remove patient 56B from annotation data
filtered_annotations <- annotations[annotations$sample != 'LP.01.56B',]
### remove patient 69 from annotation data
filtered_annotations <- filtered_annotations[filtered_annotations$sample != 'LP.01.69',]
### remove extra-thoracic samples from annotation data
filtered_annotations <- filtered_annotations[filtered_annotations$Localisation != "extra-thoracic",]
### keep only samples that are also in expression data
filtered_annotations <- filtered_annotations[filtered_annotations$sample %in% colnames(expression),]
### filter expression data columns for samples that are in the filtered annotations rows
expression <- expression[,colnames(expression) %in% rownames(filtered_annotations)]


############################################
################ TF ACTIVITY ###############
############################################

### TCGA reference for regulons
data("dorothea_hs_pancancer", package = "dorothea")
regulons_tcga = dorothea_hs_pancancer %>%
  filter(confidence %in% c("A", "B"))

tf_activities_tcga <- run_viper(expression, regulons_tcga, 
                                options =  list(method = "scale", minsize = 4,eset.filter = FALSE, cores = 1, verbose = FALSE))
#write.csv(tf_activities_tcga,"A:/Downloads/Projects/workFromHome/Projects/LungPredict/TF_analysis/TF_activities_TCGA.csv")

### GTEX reference for regulons
data("dorothea_hs", package = "dorothea")
regulons_gtex = dorothea_hs %>%
  filter(confidence %in% c("A", "B"))

tf_activities_gtex <- run_viper(expression, regulons_gtex, 
                                options =  list(method = "scale", minsize = 4,eset.filter = FALSE, cores = 1, verbose = FALSE))

#write.csv(tf_activities_gtex,"A:/Downloads/Projects/workFromHome/Projects/LungPredict/TF_analysis/TF_activities_GTEX.csv")


################################################################################
########## Comparison between results from TCGA and GTEX references  ###########
################################################################################
tf_activities_tcga_df <- as.data.frame(tf_activities_tcga)
tf_activities_gtex_df <- as.data.frame(tf_activities_gtex)
print(rownames(tf_activities_gtex_df))

## filter tf_activities_gtex_df as it has more TF estimated than tf_activities_tcga_df
tf_activities_gtex_df_filtered <- tf_activities_gtex_df[row.names(tf_activities_gtex_df) %in% row.names(tf_activities_tcga_df),]   

tf_activities_tcga_df_filtered <- tf_activities_tcga_df[row.names(tf_activities_tcga_df) %in% row.names(tf_activities_gtex_df_filtered),]   


### compute TF correlations by transposing the matrices (remove the transpose to get patient correlation)
tcga_gtex.patientcor <- cor(tf_activities_tcga_df_filtered, tf_activities_gtex_df_filtered)
?cor
heatmap(tcga_gtex.patientcor)
tcga_gtex.TFcor <- cor(t(tf_activities_tcga_df_filtered), t(tf_activities_gtex_df_filtered))

?heatmap
### compute value differences per patient for each estimated TF
diff <- tf_activities_tcga_df - tf_activities_gtex_df_filtered
diff2 <- diff[rownames(diff) != "SOX2",] %>% as.matrix()
Heatmap(diff2,
         row_names_gp = gpar(fontsize = 7),
         column_names_gp = gpar(fontsize =7),
         )
################################################################################
################################################################################
################################################################################

# drop the sample column for plotting
filtered_annotations <- subset(filtered_annotations, select = -c(sample) )


my_matrix <- as.matrix(tf_activities_tcga)



greyscale = grey.colors(10, rev=T)

# ann_colors = list(
#   Stage = c("IA-1" = greyscale[1], "IA-2" = greyscale[2], "IA-3" = greyscale[3],
#             "IB" = greyscale[4], "IIA"  = greyscale[5], "IIB" = greyscale[6],
#             "IIIA" = greyscale[7], "IIIB" = greyscale[8],"IV-A" = greyscale[9],
#             "IV-B" = greyscale[10]),
#   Mutation = c("EGFR" = "red", "No" = "white", "STK11" = "brown", "TP53" = "purple",
#                "KRAS.G12V" = "deepskyblue", "KRAS.G12D" = "deepskyblue3",
#                "KRAS.G13C" = "deepskyblue4", "KRAS.Q61E" = "blue", "KRAS" = "darkblue",
#                "KRAS.G12C" = "deepskyblue2", "NRAS"="orange", "ERBB2"="pink",
#                "CDKN2A"="darkturquoise",  "MET"="grey", "PTEN"= "aquamarine",
#                "FGFR3"="chocolate", "BRAF"="yellow", "SMAD4"="tan",  "p.E542K"="snow2"),
#   Surgery = c("No" = "darkmagenta", "Yes" = "darkgoldenrod1"),
#   Metastatic = c("Yes"="red", "No"="white"),
#   Sexe = c("M" = "blue", "F" = "pink"),
#   Smoker = c("No"="black", "Yes"="red", "Old"="grey", "Unknown"="white"),
#   Location = c("Lung" = "salmon2", "extra-thoracic" = "seagreen"),
#   Age = c("40-45" = greyscale[1], "45-50" = greyscale[2], "50-55" = greyscale[3],
#   "55-60" = greyscale[4], "60-65"  = greyscale[5], "65-70" = greyscale[6],
#   "70-75" = greyscale[7], "75-80" = greyscale[8],"80-85" = greyscale[9])
# ) 

?rowAnnotation
ann_colors = HeatmapAnnotation( Stage = filtered_annotations$Stage, 
                            Mutation = filtered_annotations$Mutation, 
                            Surgery = filtered_annotations$Surgery,
                            Metastatic = filtered_annotations$Metastatic, 
                            Sexe = filtered_annotations$Sexe, 
                            Smoker = filtered_annotations$Smoker,
                            Location = filtered_annotations$Localisation, 
                            Age = filtered_annotations$Age, 
                            col = list(
                              Stage = c("IA-1" = greyscale[1], "IA-2" = greyscale[2], "IA-3" = greyscale[3],
                                        "IB" = greyscale[4], "IIA"  = greyscale[5], "IIB" = greyscale[6],
                                        "IIIA" = greyscale[7], "IIIB" = greyscale[8],"IV-A" = greyscale[9],
                                        "IV-B" = greyscale[10]),
                              Mutation = c("EGFR" = "red", "No" = "white", "STK11" = "brown", "TP53" = "purple",
                                           "KRAS.G12V" = "deepskyblue", "KRAS.G12D" = "deepskyblue3",
                                           "KRAS.G13C" = "deepskyblue4", "KRAS.Q61E" = "blue", "KRAS" = "darkblue",
                                           "KRAS.G12C" = "deepskyblue2", "NRAS"="orange", "ERBB2"="pink",
                                           "CDKN2A"="darkturquoise",  "MET"="grey", "PTEN"= "aquamarine",
                                           "FGFR3"="chocolate", "BRAF"="yellow", "SMAD4"="tan",  "p.E542K"="snow2"),
                              Surgery = c("No" = "darkmagenta", "Yes" = "darkgoldenrod1"),
                              Metastatic = c("Yes"="red", "No"="white"),
                              Sexe = c("M" = "blue", "F" = "pink"),
                              Smoker = c("No"="black", "Yes"="red", "Old"="salmon", "Unknown"="white"),
                              Location = c("Lung" = "salmon2", "extra-thoracic" = "seagreen"),
                              Age = c("40-45" = greyscale[1], "45-50" = greyscale[2], "50-55" = greyscale[3],
                                      "55-60" = greyscale[4], "60-65"  = greyscale[5], "65-70" = greyscale[6],
                                      "70-75" = greyscale[7], "75-80" = greyscale[8],"80-85" = greyscale[9])
                            ) )
# pheatmap(my_matrix, 
#          main="LP Patient Clustering based on TF activities",
#          annotation_col = filtered_annotations, 
#          annotation_colors = ann_colors,
#          show_rownames = TRUE, 
#          #clustering_distance_rows="euclidean",
#          clustering_method = "complete",
#          scale = "column",
#          cellwidth = 5,
#          cellheight = 5,
#          fontsize_col = 5,
#          fontsize_row = 4)

?Heatmap
Heatmap(my_matrix,
        top_annotation = ann_colors, 
        name = "Patient x TF activities",
        show_row_names = T, 
        show_heatmap_legend = T, 
        row_names_gp = gpar(fontsize = 7),
        column_names_gp = gpar(fontsize =7),
        #          cellwidth = 5,
        #          cellheight = 5,
        #          fontsize_col = 5,
        #          fontsize_row = 4
)


cor.matrix <- cor(my_matrix)

Heatmap(cor.matrix, top_annotation = ann_colors, 
        name = "patient correlation based on TF activities",
        show_row_names = T, show_heatmap_legend = T, 
        row_names_gp = gpar(fontsize = 7),
        column_names_gp = gpar(fontsize =7),
)


cor.matrix.expression <- cor(expression)

# pheatmap(cor.matrix.expression, 
#          main="LP Patient Correlation based on gene expression",
#          annotation_col = filtered_annotations, 
#          annotation_colors = ann_colors,
#          show_rownames = TRUE, 
#          clustering_method = "complete",
#          cellwidth = 5,
#          cellheight = 5,
#          fontsize_col = 5,
#          fontsize_row = 4)

Heatmap(cor.matrix.expression, top_annotation = ann_colors, 
        name = "patient correlation based on expression",
        show_row_names = T, show_heatmap_legend = T, 
        row_names_gp = gpar(fontsize = 7),
        column_names_gp = gpar(fontsize =7),
)



t_matrix <- t(my_matrix)
cor.t_matrix <- cor(t_matrix)
Heatmap(cor.t_matrix, 
        name = "TF correlation based on 64 patients",
        show_row_names = T, show_heatmap_legend = T, 
        row_names_gp = gpar(fontsize = 7),
        column_names_gp = gpar(fontsize =7),
)



#### PCA
res.pca <- prcomp(expression, scale = TRUE)
install.packages("factoextra")
library(factoextra)
fviz_eig(res.pca)
fviz_pca_ind(res.pca,
             col.ind = "cos2", # Color by the quality of representation
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)

library(corbetw2mat)



## Ting
ann_col <- data.frame(as.character(resCut), FFPE_Clinical[Adenocarcinoma,]$Surgery, FFPE_Clinical[Adenocarcinoma,]$diagnoostic_VI_metastasis, FFPE_Clinical[Adenocarcinoma,]$mutation, FFPE_Clinical[Adenocarcinoma,]$Stage, FFPE_Clinical[Adenocarcinoma,]$sexe, FFPE_Clinical[Adenocarcinoma,]$smoker,FFPE_Clinical[Adenocarcinoma,]$localisation_ponction)
rownames(ann_col) <- names(resCut)
colnames(ann_col)=c('Hclust', 'Surgery', 'Metastatic', 'Mutation', 'Stage', 'Sexe', 'Smoker')
greyscale = grey.colors(10, rev=T)
ann_colors = list(
  Hclust = c("1" = "deeppink", "2" = "deepskyblue"),
  Surgery = c("No" = "darkmagenta", "Yes" = "darkgoldenrod1"),
  Metastatic = c("Yes"="red", "No"="blue"),
  Mutation = c("EGFR" = "red", "No" = "white", "STK11" = "brown", "TP53" = "purple",
               "KRAS.G12V" = "deepskyblue", "KRAS.G12D" = "deepskyblue3",
               "KRAS.G13C" = "deepskyblue4", "KRAS.Q61E" = "blue", "KRAS" = "darkblue",
               "KRAS.G12C" = "deepskyblue2", "NRAS"="orange", "ERBB2"="pink",
               "CDKN2A"="darkturquoise",  "MET"="grey", "PTEN"= "aquamarine",
               "FGFR3"="chocolate", "BRAF"="yellow", "SMAD4"="tan",  "p.E542K"="snow2"),
  Stage = c("IA-1" = greyscale[1], "IA-2" = greyscale[2], "IA-3" = greyscale[3],
            "IB" = greyscale[4], "IIA"  = greyscale[5], "IIB" = greyscale[6],
            "IIIA" = greyscale[7], "IIIB" = greyscale[8],"IV-A" = greyscale[9],
            "IV-B" = greyscale[10]),
  Sexe = c("M" = "red", "F" = "blue"),
  Smoker = c("No"="black", "Yes"="red", "Old"="grey"),
  Location = c("Lung" = "salmon2", "extra-thoracic" = "seagreen")
) 
cluster_rows = F, annotation_col = ann_row

res <- pheatmap(t(FFPE_STAR_Normalisation[rownames(variable), Adenocarcinoma]), scale = "row", fontsize_col = 10, show_rownames = T,  show_colnames = F, fontsize_row = 10, main = "Clustering based on gene expression")
resCut = cutree(res$tree_row, 2)
resc1 = names(which(resCut == 1))
resc2 = names(which(resCut == 2))

pheatmap(FFPE_STAR_Normalisation[rownames(variable), Adenocarcinoma], scale = "column", fontsize_col = 10, show_rownames = F, show_colnames = T, cluster_rows = F, annotation_col = ann_col, annotation_colors = ann_colors, fontsize_row = 5, main = "Clustering based on gene expression (Adenocarcinoma)")

















#Julien

annotZ <- rowAnnotation(Mutation = ihc1_l1_network$mutation, 
                        col = list(Mutation = c("other" = "black", "none"="white",
                                                "EGFR" = "red", "EGFR et al" = "orange",
                                                "KRAS" ="blue",  "KRAS et al"="lightblue", 
                                                "STK11" ="purple", "STK11 et al"="pink")))

Heatmap(scale(ihc1_l1_network[,2:5]), column_title = "Cell type estimation across IHC (rescaled)",
        show_column_dend = F, right_annotation = annotZ,
        name = "% of cell", show_row_names = T, show_heatmap_legend = T, row_names_gp = gpar(fontsize = 8))


## Jacobo
col_b <- colorRamp2(c(min(type_moffit["Basal", ]), max(type_moffit["Basal", ])), c("white", "chocolate1")) ## quantitative 
col_collison <- c("NA" = "white", "QM" = "orange", "Classical" = "deepskyblue1") ## category
column_ha <- HeatmapAnnotation(
  Serine = type_serine,
  PHGDH = ccle_sym.enz["PHGDH", ],
  CBS = ccle_sym.enz["CBS", ],
  PSAT1 = ccle_sym.enz["PSAT1", ],
  Classical_enr = type_moffit["Classical", ],
  Basal_enr = type_moffit["Basal", ],
  Collison = type_collison,
  col = list(
    Classical_enr = col_c,
    Basal_enr = col_b,
    Collison = col_collison,
    Serine = col_serine,
    PHGDH = col_PHGDH,
    CBS = col_CBS,
    PSAT1 = col_PSAT1
  )
)
Heatmap(ccle_sym.transclass,
        top_annotation = column_ha, # left_annotation = row_ha, #This one is only for when ccle_sym.transclass is Moffit
        cluster_columns = F
)



# ?viper
# ?msviper
# 
# run_msviper <- function(exprs_m, dorothea, use_aracne, ref, treat, ref_name, treat_name, minsize, ges.filter) {
#   # First we need to generate the phenotype table (AnnotatedDataFrame)
#   conditions <- rep("NA", ncol(exprs_m))
#   conditions[ref] <- ref_name
#   conditions[treat] <- treat_name
#   names(conditions) <- colnames(exprs_m)
#   conditions <- conditions[which(conditions != "NA")]
#   
#   phenotype <- data.frame(condition = factor(conditions))
#   rownames(phenotype) <- names(conditions)
#   
#   phenoData <- new("AnnotatedDataFrame", data = phenotype)
#   
#   exprs_m <- exprs_m[,which(colnames(exprs_m) %in% rownames(phenotype))] %>% as.matrix()
#   
#   # Create Expression set from phenotyble table and expression matrix
#   dset_viper <- ExpressionSet(assayData = exprs_m, phenoData = phenoData)
#   dset_viper$sampleID <- factor(colnames(exprs_m))
#   
#   # Aracne can be used to estimate the mor instead using the -1, 1 from dorothea
#   regulons <- NULL
#   if (use_aracne) {
#     regulons <- dorothea2aracne2viper_regulons(dorothea, dset_viper)
#   } else {
#     regulons <- dorothea2viper_regulons(dorothea)
#   }
#   
#   # We need to create the statistics signature from the conditions
#   signature <- rowTtest(dset_viper, "condition", treat_name, ref_name)
#   statistics_signature <- (qnorm(signature$p.value / 2, lower.tail = FALSE) * sign(signature$statistic))[, 1]
#   # Generate the null model with bootstrapping (1000 iterations)
#   nullmodel <- ttestNull(dset_viper, "condition", treat_name, ref_name, per = 1000, repos = T, verbose = F)
#   # Run msviper using the statistics signature, the regulons converted from dorothea table, the null model the minSize of regulon and the ges.filter
#   mrs <- msviper(ges = statistics_signature, regulon = regulons, nullmodel = nullmodel, minsize = minsize, ges.filter = ges.filter, verbose = F)
#   # Convert the msviper regulons to dorothea
#   dorothea_mrs_regulons <- viper_regulons2dorothea(mrs$regulon) %>%
#     mutate(state = ifelse(mor > 0, "activation", "inhibition"))
#   # Generate a table with the TFs, the regulon size, the NES score, the pval and the pval.fdr
#   mrs_table <- tibble(TF = names(mrs$es$p.value), size = mrs$es$size, nes = mrs$es$nes, pval = mrs$es$p.value, pval.fdr = p.adjust(mrs$es$p.value, method = "fdr")) %>% arrange(pval)
#   
#   list(mrs_table = mrs_table, mrs = mrs, regulons = dorothea_mrs_regulons)
# }

TF_activities_analysis.R 
#required packages
# install.packages("pillar")
# install.packages("ggplot2")
# install.packages("dplyr")
# install.packages("forcats")
# BiocManager::install("viper")
#BiocManager::install("dorothea")
# ??dorothea
#install.packages("devtools")
# library(devtools)
# install_github("saezlab/dorothea") # to get the TCGA (pan cancer) ref for regulons

library("tidyverse")
library("viper")
library("dorothea")
library(pheatmap)
library(RColorBrewer)
library(viridis)
library(ComplexHeatmap)

#################################
######## EXPRESSION DATA ########
#################################

#loading expression file
expression=read.csv("A:/Downloads/Projects/workFromHome/Projects/LungPredict/TF_analysis/LP_FFPE_STAR_RSEM_TPM.txt", sep = "	")
expression=expression%>% column_to_rownames(var='Gene') %>% as.matrix()

#print(colnames(expression))
### remove patient 69 from expression data
# expression <- subset(expression, select = -c(LP.01.69))

#########################################################
###### FILTER PATIENT BASED ON CLINICAL FEATURES  #######
#########################################################

### create annotations clinical feature file
clinical_features <- read.table("A:/Downloads/Projects/workFromHome/Projects/LungPredict/TF_analysis/clinic_data_v2_clean.csv", header = TRUE, sep = ",")

# look for specific mutations
egfr_mutant <- clinical_features$sample[clinical_features$mutation == "EGFR"]
kras_mutant <- clinical_features$sample[grep("KRAS", clinical_features$mutation)]
stk11_mutant <- clinical_features$sample[clinical_features$mutation == "STK11"]
`%ni%` <- Negate(`%in%`)
clinical_features <- clinical_features %>% 
  mutate(KRAS = case_when(
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
  ))


# add age category
clinical_features <- clinical_features %>%
  mutate(age_category = case_when(
    age < 85 & age >= 80 ~ '80-85',
    age < 80 & age >= 75 ~ '75-80',
    age < 75 & age >= 70 ~ '70-75',
    age < 70 & age >= 65 ~ '65-70',
    age < 65 & age >= 60 ~ '60-65',
    age < 60 & age >= 55 ~ '55-60',
    age < 55 & age >= 50 ~ '50-55',
    age < 50 & age >= 45 ~ '45-50',
    age < 45 & age >= 40 ~ '40-45'
  ))


# expBasedClusters <- read.table("A:/Downloads/Projects/workFromHome/Projects/LungPredict/TF_analysis/expression_based_patient_clusters2.txt", header = TRUE, sep = "	")
expBasedClusters <- read.table("A:/Downloads/Projects/workFromHome/Projects/LungPredict/TF_analysis/a.txt", header = TRUE, sep = "	")
colnames(expBasedClusters)[1] <- "sample"

reactomeClusters <- read.table("A:/Downloads/Projects/workFromHome/Projects/LungPredict/TF_analysis/ReactomeClustersAllPatients.csv", header = TRUE, sep = ",")
colnames(reactomeClusters)[1] <- "sample"

deconvCancerClusters <- read.table("A:/Downloads/Projects/workFromHome/Projects/LungPredict/TF_analysis/DeconvCancerClusters.txt", header = TRUE, sep = "	")
colnames(deconvCancerClusters)[1] <- "sample"




annotations <- data.frame(clinical_features$sample, 
                          clinical_features$Stage, 
                          clinical_features$Surgery, 
                          clinical_features$diagnostic,
                          clinical_features$diagnoostic_VI_metastasis, 
                          clinical_features$sexe, 
                          clinical_features$smoker, 
                          clinical_features$localisation_ponction,
                          clinical_features$age_category,
                          clinical_features$KRAS,
                          clinical_features$EGFR,
                          clinical_features$STK11)

# create the annotation DF
rownames(annotations) <- annotations$clinical_features.sample


colnames(annotations) <- c("sample", "Stage", "Surgery", "diagnostic", "Metastatic", "Sexe", "Smoker", "Localisation", 
                           "Age", "KRAS", "EGFR", "STK11")

############################################################################  
################ REMOVE UNDESIRED PATIENTS ################################  
################ AND FILTER THE ANNOTATIONS AND EXPRESSION DATA FOR OUR DATA
############################################################################ 

### remove patient 56B from annotation data
filtered_annotations <- annotations[annotations$sample != 'LP.01.56B',]
### remove patient 69 from annotation data
filtered_annotations <- filtered_annotations[filtered_annotations$sample != 'LP.01.69',]
### remove extra-thoracic samples from annotation data
filtered_annotations <- filtered_annotations[filtered_annotations$Localisation != "extra-thoracic",]
### keep only the adenocarcinoma ??
#### filtered_annotations <- filtered_annotations[filtered_annotations$diagnostic == "Adenocarcinoma",]


#Reduce(function(x, y) merge(x, y, by = "sample"), list(df1, df2, df3))

?right_join

outputAnnotations <-  right_join(expBasedClusters, filtered_annotations, by = "sample")
colnames(outputAnnotations)[2] <- "GE_Cluster"


outputAnnotations <-  right_join(reactomeClusters, outputAnnotations, by = "sample")
colnames(outputAnnotations)[2] <- "reactomeCluster"

outputAnnotations2 <-  right_join(deconvCancerClusters[,1:2], outputAnnotations, by = "sample")
outputAnnotations3 <-  right_join(deconvCancerClusters[,c(1,3)], outputAnnotations2, by = "sample")

outputAnnotations4 <- outputAnnotations3[,-12]



write.table(outputAnnotations4, file= "A:/Downloads/Projects/workFromHome/Projects/LungPredict/TF_analysis/full_annotations_with_clusters_corrected1.txt", sep="	", quote=F, row.names=FALSE)






### keep only samples that are also in expression data
filtered_annotations <- filtered_annotations[filtered_annotations$sample %in% colnames(expression),]



### filter expression data columns for samples that are in the filtered annotations rows
expression <- expression[,colnames(expression) %in% rownames(filtered_annotations)]


# create the annotation
greyscale <- grey.colors(10, rev = T)
ann_colors <- HeatmapAnnotation(
  Stage = filtered_annotations$Stage, Surgery = filtered_annotations$Surgery,
  Age = filtered_annotations$Age,
  Metastatic = filtered_annotations$Metastatic, 
  Sexe = filtered_annotations$Sexe, 
  Smoker = filtered_annotations$Smoker,
  Location = filtered_annotations$Localisation, 
  Age = filtered_annotations$Age, 
  KRAS = filtered_annotations$KRAS,
  EGFR = filtered_annotations$EGFR,
  STK11 = filtered_annotations$STK11,
  col = list(
    Stage = c(
      "IA-1" = greyscale[1], "IA-2" = greyscale[2], "IA-3" = greyscale[3],
      "IB" = greyscale[4], "IIA" = greyscale[5], "IIB" = greyscale[6],
      "IIIA" = greyscale[7], "IIIB" = greyscale[8], "IV-A" = greyscale[9],
      "IV-B" = greyscale[10]
    ),
    Surgery = c("No" = "darkmagenta", "Yes" = "darkgoldenrod1"),
    Metastatic = c("Yes" = "red", "No" = "white"),
    Sexe = c("M" = "blue", "F" = "pink"),
    Smoker = c("No" = "black", "Yes" = "red", "Old" = "salmon", "Unknown" = "white"),
    Location = c("Lung" = "salmon2", "extra-thoracic" = "seagreen"),
    Age = c(
      "40-45" = greyscale[1], "45-50" = greyscale[2], "50-55" = greyscale[3],
      "55-60" = greyscale[4], "60-65" = greyscale[5], "65-70" = greyscale[6],
      "70-75" = greyscale[7], "75-80" = greyscale[8], "80-85" = greyscale[9]
    ),
    KRAS = c("yes" = "blue", "no" = "white"),
    EGFR = c("yes" = "red", "no" = "white"),
    STK11 = c("yes" = "green", "no" = "white")
  )
)

############################################
## GENE EXPRESSION BASED CLUSTERING ##
############################################

cor.matrix.expression <- cor(expression)

# myPheatmap <- pheatmap(cor.matrix.expression,
#          main="LP Patient Correlation based on gene expression",
#          annotation_col = filtered_annotations,
#          annotation_colors = ann_colors,
#          show_rownames = TRUE,
#          clustering_method = "complete",
#          cellwidth = 5,
#          cellheight = 5,
#          fontsize_col = 5,
#          fontsize_row = 4)

# drop the sample column for plotting
filtered_annotations <- subset(filtered_annotations, select = -c(sample) )

myExpBasedPatientCorrHeatmap <- Heatmap(cor.matrix.expression, top_annotation = ann_colors, 
        name = "patient correlation based on expression",
        show_row_names = T, show_heatmap_legend = T, 
        row_names_gp = gpar(fontsize = 8),
        column_names_gp = gpar(fontsize = 8)
)

print(myExpBasedPatientCorrHeatmap)

### Extract clusters from the heatmap
HM = Heatmap(cor.matrix.expression, top_annotation = ann_colors, 
             name = "patient correlation based on expression",
             show_row_names = T, show_heatmap_legend = T, 
             row_names_gp = gpar(fontsize = 8),
             column_names_gp = gpar(fontsize = 8), 
             km = 2)
?Heatmap

hc <- hclust(dist(cor.matrix.expression)) 
resCut <- cutree(hc, k = 3)
resCut.dataframe <- as.data.frame(resCut)
# resCut.dataframe$sample <- rownames(resCut.dataframe)

write.table(data.frame("sample"=rownames(resCut.dataframe), resCut.dataframe), "A:/Downloads/Projects/workFromHome/Projects/LungPredict/TF_analysis/a.txt", row.names=FALSE, quote=F, sep="	")



write.table(resCut.dataframe, file= "A:/Downloads/Projects/workFromHome/Projects/LungPredict/TF_analysis/expression_based_patient_clusters_resCutCorrected.txt", sep="	", quote=F, row.names=TRUE)





HM = draw(HM)
r.dend <- row_dend(HM)  #Extract row dendrogram
rcl.list <- row_order(HM)  #Extract clusters (output is a list)
lapply(rcl.list, function(x) length(x))  #check/confirm size clusters

# loop to extract patient_id for each cluster.
for (i in 1:length(row_order(HM))){
  if (i == 1) {
    clu <- t(t(row.names(cor.matrix.expression[row_order(HM)[[i]],])))
    out <- cbind(clu, paste("Cluster", i, sep=""))
    colnames(out) <- c("PatientID", "Cluster")
    } else {
      clu <- t(t(row.names(cor.matrix.expression[row_order(HM)[[i]],])))
      clu <- cbind(clu, paste("Cluster", i, sep=""))
      out <- rbind(out, clu)
      }
}

#check
out 

#export
write.table(out, file= "A:/Downloads/Projects/workFromHome/Projects/LungPredict/TF_analysis/expression_based_patient_clusters2.txt", sep="	", quote=F, row.names=FALSE)

############################################
################ TF ACTIVITY ###############
############################################

### TCGA reference for regulons
data("dorothea_hs_pancancer", package = "dorothea")
regulons_tcga = dorothea_hs_pancancer %>%
  filter(confidence %in% c("A", "B"))

tf_activities_tcga <- run_viper(expression, regulons_tcga, 
                                options =  list(method = "scale", minsize = 4,eset.filter = FALSE, cores = 1, verbose = FALSE))
#write.csv(tf_activities_tcga,"A:/Downloads/Projects/workFromHome/Projects/LungPredict/TF_analysis/TF_activities_TCGA.csv")

### GTEX reference for regulons
data("dorothea_hs", package = "dorothea")
regulons_gtex = dorothea_hs %>%
  filter(confidence %in% c("A", "B"))

tf_activities_gtex <- run_viper(expression, regulons_gtex, 
                                options =  list(method = "scale", minsize = 4,eset.filter = FALSE, cores = 1, verbose = FALSE))

#write.csv(tf_activities_gtex,"A:/Downloads/Projects/workFromHome/Projects/LungPredict/TF_analysis/TF_activities_GTEX.csv")


################################################################################
########## Comparison between results from TCGA and GTEX references  ###########
################################################################################
tf_activities_tcga_df <- as.data.frame(tf_activities_tcga)
tf_activities_gtex_df <- as.data.frame(tf_activities_gtex)
print(rownames(tf_activities_gtex_df))

## filter tf_activities_gtex_df as it has more TF estimated than tf_activities_tcga_df
tf_activities_gtex_df_filtered <- tf_activities_gtex_df[row.names(tf_activities_gtex_df) %in% row.names(tf_activities_tcga_df),]   

tf_activities_tcga_df_filtered <- tf_activities_tcga_df[row.names(tf_activities_tcga_df) %in% row.names(tf_activities_gtex_df_filtered),]   


### compute TF correlations by transposing the matrices (remove the transpose to get patient correlation)
tcga_gtex.patientcor <- cor(tf_activities_tcga_df_filtered, tf_activities_gtex_df_filtered)
?cor
heatmap(tcga_gtex.patientcor)
tcga_gtex.TFcor <- cor(t(tf_activities_tcga_df_filtered), t(tf_activities_gtex_df_filtered))

?heatmap
### compute value differences per patient for each estimated TF
diff <- tf_activities_tcga_df - tf_activities_gtex_df_filtered
diff2 <- diff[rownames(diff) != "SOX2",] %>% as.matrix()
Heatmap(diff2,
         row_names_gp = gpar(fontsize = 7),
         column_names_gp = gpar(fontsize =7),
         )
################################################################################
################################################################################
################################################################################

# drop the sample column for plotting
filtered_annotations <- subset(filtered_annotations, select = -c(sample) )


my_matrix_TF <- as.matrix(tf_activities_tcga)


################################################################################
######### REMOVE TF THAT DONT VARY A LOT #######################################
################################################################################

TFstd <- apply(my_matrix_TF, 1, sd)
# print(TFstd)
hist(TFstd, breaks=106)

### Remove all TF that have std < 2

my_matrix_TF_filtered <- my_matrix_TF[TFstd > 2,]

################################################################################
##### UPDATE THE ANNOTATIONS WITH EXPRESSION BASED CLUSTERS #####
### AND THE DECONV CLUSTERS
#### AND THE REACTOME CLUSTERS
################################################################################

expBasedClusters <- read.table("A:/Downloads/Projects/workFromHome/Projects/LungPredict/TF_analysis/expression_based_patient_clusters.txt", header = TRUE, sep = "	")
colnames(expBasedClusters)[1] <- "sample"

reactomeClusters <- read.table("A:/Downloads/Projects/workFromHome/Projects/LungPredict/TF_analysis/Reactomes_cluster.txt", header = TRUE, sep = "	")
colnames(expBasedClusters)[1] <- "sample"

mergedAnnotations <- inner_join(expBasedClusters, annotations, by = "sample")
# ?inner_join

ann_colors_merged <- HeatmapAnnotation(
  Stage = mergedAnnotations$Stage, Surgery = mergedAnnotations$Surgery,
  Age = mergedAnnotations$Age,
  Metastatic = mergedAnnotations$Metastatic, 
  Sexe = mergedAnnotations$Sexe, 
  Smoker = mergedAnnotations$Smoker,
  Location = mergedAnnotations$Localisation, 
  Age = mergedAnnotations$Age, 
  KRAS = mergedAnnotations$KRAS,
  EGFR = mergedAnnotations$EGFR,
  STK11 = mergedAnnotations$STK11,
  Cluster = mergedAnnotations$Cluster,
  col = list(
    Stage = c(
      "IA-1" = greyscale[1], "IA-2" = greyscale[2], "IA-3" = greyscale[3],
      "IB" = greyscale[4], "IIA" = greyscale[5], "IIB" = greyscale[6],
      "IIIA" = greyscale[7], "IIIB" = greyscale[8], "IV-A" = greyscale[9],
      "IV-B" = greyscale[10]
    ),
    Surgery = c("No" = "darkmagenta", "Yes" = "darkgoldenrod1"),
    Metastatic = c("Yes" = "red", "No" = "white"),
    Sexe = c("M" = "blue", "F" = "pink"),
    Smoker = c("No" = "black", "Yes" = "red", "Old" = "salmon", "Unknown" = "white"),
    Location = c("Lung" = "salmon2", "extra-thoracic" = "seagreen"),
    Age = c(
      "40-45" = greyscale[1], "45-50" = greyscale[2], "50-55" = greyscale[3],
      "55-60" = greyscale[4], "60-65" = greyscale[5], "65-70" = greyscale[6],
      "70-75" = greyscale[7], "75-80" = greyscale[8], "80-85" = greyscale[9]
    ),
    KRAS = c("yes" = "blue", "no" = "white"),
    EGFR = c("yes" = "red", "no" = "white"),
    STK11 = c("yes" = "green", "no" = "white"),
    Cluster = c("Cluster1" = "pink", "Cluster2" = "blue")
  )
)

Heatmap(my_matrix_TF_filtered,
        top_annotation = ann_colors_merged, 
        name = "Patient x TF activities",
        show_row_names = T, 
        show_heatmap_legend = T, 
        row_names_gp = gpar(fontsize = 7),
        column_names_gp = gpar(fontsize =7),
        #          cellwidth = 5,
        #          cellheight = 5,
        #          fontsize_col = 5,
        #          fontsize_row = 4
)

################################################################################
###### TF-BASED PATIENT CORRELATION MATRIX #####
################################################################################

cor.matrix.TF.filtered <- cor(my_matrix_TF_filtered)

Heatmap(cor.matrix.TF.filtered, top_annotation = ann_colors_merged, 
        name = "patient correlation based on TF activities",
        show_row_names = T, show_heatmap_legend = T, 
        row_names_gp = gpar(fontsize = 7),
        column_names_gp = gpar(fontsize =7),
)




t_matrix <- t(my_matrix)
cor.t_matrix <- cor(t_matrix)
Heatmap(cor.t_matrix, 
        name = "TF correlation based on 64 patients",
        show_row_names = T, show_heatmap_legend = T, 
        row_names_gp = gpar(fontsize = 7),
        column_names_gp = gpar(fontsize =7),
)



#### PCA
res.pca <- prcomp(expression, scale = TRUE)
install.packages("factoextra")
library(factoextra)
fviz_eig(res.pca)
fviz_pca_ind(res.pca,
             col.ind = "cos2", # Color by the quality of representation
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)

library(corbetw2mat)



## Ting
ann_col <- data.frame(as.character(resCut), FFPE_Clinical[Adenocarcinoma,]$Surgery, FFPE_Clinical[Adenocarcinoma,]$diagnoostic_VI_metastasis, FFPE_Clinical[Adenocarcinoma,]$mutation, FFPE_Clinical[Adenocarcinoma,]$Stage, FFPE_Clinical[Adenocarcinoma,]$sexe, FFPE_Clinical[Adenocarcinoma,]$smoker,FFPE_Clinical[Adenocarcinoma,]$localisation_ponction)
rownames(ann_col) <- names(resCut)
colnames(ann_col)=c('Hclust', 'Surgery', 'Metastatic', 'Mutation', 'Stage', 'Sexe', 'Smoker')
greyscale = grey.colors(10, rev=T)
ann_colors = list(
  Hclust = c("1" = "deeppink", "2" = "deepskyblue"),
  Surgery = c("No" = "darkmagenta", "Yes" = "darkgoldenrod1"),
  Metastatic = c("Yes"="red", "No"="blue"),
  Mutation = c("EGFR" = "red", "No" = "white", "STK11" = "brown", "TP53" = "purple",
               "KRAS.G12V" = "deepskyblue", "KRAS.G12D" = "deepskyblue3",
               "KRAS.G13C" = "deepskyblue4", "KRAS.Q61E" = "blue", "KRAS" = "darkblue",
               "KRAS.G12C" = "deepskyblue2", "NRAS"="orange", "ERBB2"="pink",
               "CDKN2A"="darkturquoise",  "MET"="grey", "PTEN"= "aquamarine",
               "FGFR3"="chocolate", "BRAF"="yellow", "SMAD4"="tan",  "p.E542K"="snow2"),
  Stage = c("IA-1" = greyscale[1], "IA-2" = greyscale[2], "IA-3" = greyscale[3],
            "IB" = greyscale[4], "IIA"  = greyscale[5], "IIB" = greyscale[6],
            "IIIA" = greyscale[7], "IIIB" = greyscale[8],"IV-A" = greyscale[9],
            "IV-B" = greyscale[10]),
  Sexe = c("M" = "red", "F" = "blue"),
  Smoker = c("No"="black", "Yes"="red", "Old"="grey"),
  Location = c("Lung" = "salmon2", "extra-thoracic" = "seagreen")
) 
cluster_rows = F, annotation_col = ann_row

res <- pheatmap(t(FFPE_STAR_Normalisation[rownames(variable), Adenocarcinoma]), scale = "row", fontsize_col = 10, show_rownames = T,  show_colnames = F, fontsize_row = 10, main = "Clustering based on gene expression")
resCut = cutree(res$tree_row, 2)
resc1 = names(which(resCut == 1))
resc2 = names(which(resCut == 2))

pheatmap(FFPE_STAR_Normalisation[rownames(variable), Adenocarcinoma], scale = "column", fontsize_col = 10, show_rownames = F, show_colnames = T, cluster_rows = F, annotation_col = ann_col, annotation_colors = ann_colors, fontsize_row = 5, main = "Clustering based on gene expression (Adenocarcinoma)")

















#Julien

annotZ <- rowAnnotation(Mutation = ihc1_l1_network$mutation, 
                        col = list(Mutation = c("other" = "black", "none"="white",
                                                "EGFR" = "red", "EGFR et al" = "orange",
                                                "KRAS" ="blue",  "KRAS et al"="lightblue", 
                                                "STK11" ="purple", "STK11 et al"="pink")))

Heatmap(scale(ihc1_l1_network[,2:5]), column_title = "Cell type estimation across IHC (rescaled)",
        show_column_dend = F, right_annotation = annotZ,
        name = "% of cell", show_row_names = T, show_heatmap_legend = T, row_names_gp = gpar(fontsize = 8))


## Jacobo
col_b <- colorRamp2(c(min(type_moffit["Basal", ]), max(type_moffit["Basal", ])), c("white", "chocolate1")) ## quantitative 
col_collison <- c("NA" = "white", "QM" = "orange", "Classical" = "deepskyblue1") ## category
column_ha <- HeatmapAnnotation(
  Serine = type_serine,
  PHGDH = ccle_sym.enz["PHGDH", ],
  CBS = ccle_sym.enz["CBS", ],
  PSAT1 = ccle_sym.enz["PSAT1", ],
  Classical_enr = type_moffit["Classical", ],
  Basal_enr = type_moffit["Basal", ],
  Collison = type_collison,
  col = list(
    Classical_enr = col_c,
    Basal_enr = col_b,
    Collison = col_collison,
    Serine = col_serine,
    PHGDH = col_PHGDH,
    CBS = col_CBS,
    PSAT1 = col_PSAT1
  )
)
Heatmap(ccle_sym.transclass,
        top_annotation = column_ha, # left_annotation = row_ha, #This one is only for when ccle_sym.transclass is Moffit
        cluster_columns = F
)



# ?viper
# ?msviper
# 
# run_msviper <- function(exprs_m, dorothea, use_aracne, ref, treat, ref_name, treat_name, minsize, ges.filter) {
#   # First we need to generate the phenotype table (AnnotatedDataFrame)
#   conditions <- rep("NA", ncol(exprs_m))
#   conditions[ref] <- ref_name
#   conditions[treat] <- treat_name
#   names(conditions) <- colnames(exprs_m)
#   conditions <- conditions[which(conditions != "NA")]
#   
#   phenotype <- data.frame(condition = factor(conditions))
#   rownames(phenotype) <- names(conditions)
#   
#   phenoData <- new("AnnotatedDataFrame", data = phenotype)
#   
#   exprs_m <- exprs_m[,which(colnames(exprs_m) %in% rownames(phenotype))] %>% as.matrix()
#   
#   # Create Expression set from phenotyble table and expression matrix
#   dset_viper <- ExpressionSet(assayData = exprs_m, phenoData = phenoData)
#   dset_viper$sampleID <- factor(colnames(exprs_m))
#   
#   # Aracne can be used to estimate the mor instead using the -1, 1 from dorothea
#   regulons <- NULL
#   if (use_aracne) {
#     regulons <- dorothea2aracne2viper_regulons(dorothea, dset_viper)
#   } else {
#     regulons <- dorothea2viper_regulons(dorothea)
#   }
#   
#   # We need to create the statistics signature from the conditions
#   signature <- rowTtest(dset_viper, "condition", treat_name, ref_name)
#   statistics_signature <- (qnorm(signature$p.value / 2, lower.tail = FALSE) * sign(signature$statistic))[, 1]
#   # Generate the null model with bootstrapping (1000 iterations)
#   nullmodel <- ttestNull(dset_viper, "condition", treat_name, ref_name, per = 1000, repos = T, verbose = F)
#   # Run msviper using the statistics signature, the regulons converted from dorothea table, the null model the minSize of regulon and the ges.filter
#   mrs <- msviper(ges = statistics_signature, regulon = regulons, nullmodel = nullmodel, minsize = minsize, ges.filter = ges.filter, verbose = F)
#   # Convert the msviper regulons to dorothea
#   dorothea_mrs_regulons <- viper_regulons2dorothea(mrs$regulon) %>%
#     mutate(state = ifelse(mor > 0, "activation", "inhibition"))
#   # Generate a table with the TFs, the regulon size, the NES score, the pval and the pval.fdr
#   mrs_table <- tibble(TF = names(mrs$es$p.value), size = mrs$es$size, nes = mrs$es$nes, pval = mrs$es$p.value, pval.fdr = p.adjust(mrs$es$p.value, method = "fdr")) %>% arrange(pval)
#   
#   list(mrs_table = mrs_table, mrs = mrs, regulons = dorothea_mrs_regulons)
# }

lungpredict_Rscripts.R 

Gene expression and TF activities plots (Nina).R 
annotations <- read.table("A:/Downloads/Projects/workFromHome/Projects/LungPredict/TF_analysis/full_annotations_with_clusters_corrected1.txt", sep="	", header=T)


#loading expression file
expression=read.csv("A:/Downloads/Projects/workFromHome/Projects/LungPredict/TF_analysis/LP_FFPE_STAR_RSEM_TPM.txt", sep = "	")
expression=expression%>% column_to_rownames(var='Gene') %>% as.matrix()

selector <- intersect(annotations$sample, colnames(expression))

# filter expression for the annotated patients
expression  <- expression[,selector]

### keep only samples that are also in expression
rownames(annotations) <- annotations$sample

# filtered_annotations <- annotations[annotations$sample %in% colnames(expression),]
filtered_annotations <- annotations[selector,]
annotations <- filtered_annotations

# create the annotation
greyscale <- grey.colors(10, rev = T)

ann_colors <- HeatmapAnnotation(
  Stage = annotations$Stage, 
  Surgery = annotations$Surgery,
  Diagnostic = annotations$diagnostic,
  Metastatic = annotations$Metastatic,
  Sex = annotations$Sexe,
  Smoker = annotations$Smoker,
  Age = annotations$Age,
  KRAS = annotations$KRAS,
  EGFR = annotations$EGFR,
  STK11 = annotations$STK11,
  Cancer_content = annotations$Cancer.clust,
  Deconv_clust = annotations$Dclust,
  Reactome_clust = annotations$reactomeCluster,
  Gene_clust = annotations$GE_Cluster,
  
  col = list(
    Stage = c(
      "IA-1" = greyscale[1], "IA-2" = greyscale[2], "IA-3" = greyscale[3],
      "IB" = greyscale[4], "IIA" = greyscale[5], "IIB" = greyscale[6],
      "IIIA" = greyscale[7], "IIIB" = greyscale[8], "IV-A" = greyscale[9],
      "IV-B" = greyscale[10]
    ),
    KRAS = c("yes" = "blue", "no" = "white"),
    EGFR = c("yes" = "red", "no" = "white"),
    STK11 = c("yes" = "green", "no" = "white"),
    Surgery = c("No" = "darkmagenta", "Yes" = "darkgoldenrod1"),
    Diagnostic = c("Adenocarcinoma" = "blue", "epidermoid carcinomas" = "pink", "Neuroendocrine tumors" = "black", "Other" = "white"),
    Metastatic = c("Yes" = "red", "No" = "white"),
    Sex = c("M" = "blue", "F" = "pink"),
    Smoker = c("No" = "black", "Yes" = "red", "Old" = "salmon", "Unknown" = "white"),
    Age = c(
      "40-45" = greyscale[1], "45-50" = greyscale[2], "50-55" = greyscale[3],
      "55-60" = greyscale[4], "60-65" = greyscale[5], "65-70" = greyscale[6],
      "70-75" = greyscale[7], "75-80" = greyscale[8], "80-85" = greyscale[9]
    ),
    Cancer_content = c("high" = "lightskyblue", "low" = "lightyellow"),
    Deconv_clust = c("c1" = "antiquewhite", "c2" = "antiquewhite4" ),
    Reactome_clust = c( "1" = "blue", "2" = "pink") ,
    Gene_clust = c("1" = "chocolate1", "2" = "chocolate4", "3" = "chocolate3")
  )
)

# #plot the matrix gene x patient
# ?Heatmap
Heatmap(expression,
        cluster_rows = FALSE,
        cluster_row_slices = FALSE,
        show_row_dend = FALSE,
        top_annotation = ann_colors,
        name = "Patient Correlation based on Expression",
        show_row_names = F,
        show_heatmap_legend = T,
        row_names_gp = gpar(fontsize = 7),
        column_names_gp = gpar(fontsize =7),
)

# #plot the matrix patient x patient correlation based on gene expression

library(Hmisc)
cor.matrix.expression2 <- rcorr(as.matrix(expression))

cor.matrix.expression2$r[cor.matrix.expression2$P > 0.05] <- 0
colors <- circlize::colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))

Heatmap(cor.matrix.expression2$r, col = colors,
        top_annotation = ann_colors, 
        name = "GE-PatientCorr",
        show_row_names = T, 
        show_heatmap_legend = T, 
        row_names_gp = gpar(fontsize = 7),
        column_names_gp = gpar(fontsize =7),
)


############################################
################ TF ACTIVITY ###############
############################################

### TCGA reference for regulons
data("dorothea_hs_pancancer", package = "dorothea")
regulons_tcga = dorothea_hs_pancancer %>%
  filter(confidence %in% c("A", "B"))

tf_activities_tcga <- run_viper(expression, regulons_tcga, 
                                options =  list(method = "scale", minsize = 4,eset.filter = FALSE, cores = 1, verbose = FALSE))



## plot patient correlation based on TF activities
cor.matrix.TF.filtered <- cor(tf_activities_tcga)

Heatmap(cor.matrix.TF.filtered, 
        top_annotation = ann_colors, 
        name = "patient correlation based on TF activities",
        show_row_names = T, show_heatmap_legend = T, 
        row_names_gp = gpar(fontsize = 7),
        column_names_gp = gpar(fontsize =7),
)

## plot TF x patient
Heatmap(tf_activities_tcga, 
        top_annotation = ann_colors, 
        name = "patient x TF activities not filtered",
        show_row_names = T, show_heatmap_legend = T, 
        row_names_gp = gpar(fontsize = 7),
        column_names_gp = gpar(fontsize =7),
)

## keep TF that vary a lot only
TFstd <- apply(tf_activities_tcga, 1, sd)
# print(TFstd)
hist(TFstd, breaks=106)

### Remove all TF that have std < 2
my_matrix_TF_filtered <- tf_activities_tcga[TFstd > 2,]


Heatmap(my_matrix_TF_filtered, 
        top_annotation = ann_colors, 
        name = "patient x TF activities not filtered",
        show_row_names = T, show_heatmap_legend = T, 
        row_names_gp = gpar(fontsize = 7),
        column_names_gp = gpar(fontsize =7),
)

Heatmaps_annotations.R 
annotations <- read.table("A:/Downloads/Projects/workFromHome/Projects/LungPredict/TF_analysis/full_annotations_with_clusters_corrected1.txt", sep="	", header=T)



#loading expression file
expression=read.csv("A:/Downloads/Projects/workFromHome/Projects/LungPredict/TF_analysis/LP_FFPE_STAR_RSEM_TPM.txt", sep = "	")
expression=expression%>% column_to_rownames(var='Gene') %>% as.matrix()

selector <- intersect(annotations$sample, colnames(expression))

# filter expression for the annotated patients
expression  <- expression[,selector]

### keep only samples that are also in expression
rownames(annotations) <- annotations$sample

# filtered_annotations <- annotations[annotations$sample %in% colnames(expression),]
filtered_annotations <- annotations[selector,]
annotations <- filtered_annotations

# create the annotation
greyscale <- grey.colors(10, rev = T)

ann_colors <- HeatmapAnnotation(
  Stage = annotations$Stage, 
  Surgery = annotations$Surgery,
  Diagnostic = annotations$diagnostic,
  Metastatic = annotations$Metastatic,
  Sex = annotations$Sexe,
  Smoker = annotations$Smoker,
  Age = annotations$Age,
  KRAS = annotations$KRAS,
  EGFR = annotations$EGFR,
  STK11 = annotations$STK11,
  Cancer_content = annotations$Cancer.clust,
  Deconv_clust = annotations$Dclust,
  Reactome_clust = annotations$reactomeCluster,
  Gene_clust = annotations$GE_Cluster,
  
  col = list(
    Stage = c(
      "IA-1" = greyscale[1], "IA-2" = greyscale[2], "IA-3" = greyscale[3],
      "IB" = greyscale[4], "IIA" = greyscale[5], "IIB" = greyscale[6],
      "IIIA" = greyscale[7], "IIIB" = greyscale[8], "IV-A" = greyscale[9],
      "IV-B" = greyscale[10]
    ),
    KRAS = c("yes" = "blue", "no" = "white"),
    EGFR = c("yes" = "red", "no" = "white"),
    STK11 = c("yes" = "green", "no" = "white"),
    Surgery = c("No" = "darkmagenta", "Yes" = "darkgoldenrod1"),
    Diagnostic = c("Adenocarcinoma" = "blue", "epidermoid carcinomas" = "pink", "Neuroendocrine tumors" = "black", "Other" = "white"),
    Metastatic = c("Yes" = "red", "No" = "white"),
    Sex = c("M" = "blue", "F" = "pink"),
    Smoker = c("No" = "black", "Yes" = "red", "Old" = "salmon", "Unknown" = "white"),
    Age = c(
      "40-45" = greyscale[1], "45-50" = greyscale[2], "50-55" = greyscale[3],
      "55-60" = greyscale[4], "60-65" = greyscale[5], "65-70" = greyscale[6],
      "70-75" = greyscale[7], "75-80" = greyscale[8], "80-85" = greyscale[9]
    ),
    Cancer_content = c("high" = "lightskyblue", "low" = "lightyellow"),
    Deconv_clust = c("c1" = "antiquewhite", "c2" = "antiquewhite4" ),
    Reactome_clust = c( "1" = "blue", "2" = "pink") ,
    Gene_clust = c("1" = "chocolate1", "2" = "chocolate4", "3" = "chocolate3")
  )
)

# #plot the matrix gene x patient
# ?Heatmap
Heatmap(expression,
        cluster_rows = FALSE,
        cluster_row_slices = FALSE,
        show_row_dend = FALSE,
        top_annotation = ann_colors,
        name = "Patient Correlation based on Expression",
        show_row_names = F,
        show_heatmap_legend = T,
        row_names_gp = gpar(fontsize = 7),
        column_names_gp = gpar(fontsize =7),
)

# #plot the matrix patient x patient correlation based on gene expression
# cor.matrix.expression <- cor(expression)

library(Hmisc)
cor.matrix.expression2 <- rcorr(as.matrix(expression))

cor.matrix.expression2$r[cor.matrix.expression2$P > 0.05] <- 0
colors <- circlize::colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))

Heatmap(cor.matrix.expression2$r, col = colors,
        top_annotation = ann_colors, 
        name = "GE-PatientCorr",
        show_row_names = T, 
        show_heatmap_legend = T, 
        row_names_gp = gpar(fontsize = 7),
        column_names_gp = gpar(fontsize =7),
)



# my_matrix <- cor.matrix.expression
# 
# my_matrix_filtered <- my_matrix[selector , selector]
# 
# 
# 
# 
# # annotations <- annotations[,selector]
# # annotations <- annotations[,-1]
# # annotations <- subset(annotations, select = -c(sample) )
# 
# 
# 
# 
# 
# 
# Heatmap(my_matrix_filtered,
#         top_annotation = ann_colors, 
#         name = "Patient Correlation based on Expression",
#         show_row_names = T, 
#         show_heatmap_legend = T, 
#         row_names_gp = gpar(fontsize = 7),
#         column_names_gp = gpar(fontsize =7),
# )


############################################
################ TF ACTIVITY ###############
############################################

### TCGA reference for regulons
data("dorothea_hs_pancancer", package = "dorothea")
regulons_tcga = dorothea_hs_pancancer %>%
  filter(confidence %in% c("A", "B"))

tf_activities_tcga <- run_viper(expression, regulons_tcga, 
                                options =  list(method = "scale", minsize = 4,eset.filter = FALSE, cores = 1, verbose = FALSE))



## plot patient correlation based on TF activities
cor.matrix.TF.filtered <- cor(tf_activities_tcga)

Heatmap(cor.matrix.TF.filtered, 
        top_annotation = ann_colors, 
        name = "patient correlation based on TF activities",
        show_row_names = T, show_heatmap_legend = T, 
        row_names_gp = gpar(fontsize = 7),
        column_names_gp = gpar(fontsize =7),
)

## plot TF x patient
Heatmap(tf_activities_tcga, 
        top_annotation = ann_colors, 
        name = "patient x TF activities not filtered",
        show_row_names = T, show_heatmap_legend = T, 
        row_names_gp = gpar(fontsize = 7),
        column_names_gp = gpar(fontsize =7),
)

## keep TF that vary a lot only
TFstd <- apply(tf_activities_tcga, 1, sd)
# print(TFstd)
hist(TFstd, breaks=106)

### Remove all TF that have std < 2
my_matrix_TF_filtered <- tf_activities_tcga[TFstd > 2,]


Heatmap(my_matrix_TF_filtered, 
        top_annotation = ann_colors, 
        name = "patient x TF activities not filtered",
        show_row_names = T, show_heatmap_legend = T, 
        row_names_gp = gpar(fontsize = 7),
        column_names_gp = gpar(fontsize =7),
)










?Heatmap








# then your heatmap
Heatmap(whatever_data_you_use,
        top_annotation = ann_colors, # annotations just created
        name = "Name on the main legend",
        show_heatmap_legend = T, 
        row_names_gp = gpar(fontsize = 8), column_names_gp = gpar(fontsize = 8) # set size of row and columns names
)

TF_activities.R 
#required packages
# install.packages("pillar")
# install.packages("ggplot2")
# install.packages("dplyr")
# install.packages("forcats")
# BiocManager::install("viper")
#BiocManager::install("dorothea")
# ??dorothea
#install.packages("devtools")
# library(devtools)
# install_github("saezlab/dorothea") # to get the TCGA (pan cancer) ref for regulons

library("tidyverse")
library("viper")
library("dorothea")
library(pheatmap)
library(RColorBrewer)
library(viridis)
library(ComplexHeatmap)

#################################
######## EXPRESSION DATA ########
#################################

#loading expression file
expression=read.csv("A:/Downloads/Projects/workFromHome/Projects/LungPredict/TF_analysis/LP_FFPE_STAR_RSEM_TPM.txt", sep = "	")
expression=expression%>% column_to_rownames(var='Gene') %>% as.matrix()
print(colnames(expression))
### remove patient 69 from expression data
# expression <- subset(expression, select = -c(LP.01.69))

#########################################################
###### FILTER PATIENT BASED ON CLINICAL FEATURES  #######
#########################################################

### create annotations
clinical_features <- read.table("A:/Downloads/Projects/workFromHome/Projects/LungPredict/TF_analysis/clinic_data_v2_clean.csv", header = TRUE, sep = ",")
# add age category
clinical_features <- clinical_features %>%
  mutate(age_category = case_when(
    age < 85 & age >= 80 ~ '80-85',
    age < 80 & age >= 75 ~ '75-80',
    age < 75 & age >= 70 ~ '70-75',
    age < 70 & age >= 65 ~ '65-70',
    age < 65 & age >= 60 ~ '60-65',
    age < 60 & age >= 55 ~ '55-60',
    age < 55 & age >= 50 ~ '50-55',
    age < 50 & age >= 45 ~ '45-50',
    age < 45 & age >= 40 ~ '40-45'
  ))

annotations <- data.frame(clinical_features$sample, 
                          clinical_features$Stage, 
                          clinical_features$mutation, 
                          clinical_features$Surgery, 
                          clinical_features$diagnoostic_VI_metastasis, 
                          clinical_features$sexe, 
                          clinical_features$smoker, 
                          clinical_features$localisation_ponction,
                          clinical_features$age_category)
rownames(annotations) <- annotations$clinical_features.sample


colnames(annotations) <- c('sample', 'Stage', 'Mutation', 'Surgery', 'Metastatic', 'Sexe', 'Smoker', 'Localisation', 'Age')


### remove patient 56B from annotation data
filtered_annotations <- annotations[annotations$sample != 'LP.01.56B',]
### remove patient 69 from annotation data
filtered_annotations <- filtered_annotations[filtered_annotations$sample != 'LP.01.69',]
### remove extra-thoracic samples from annotation data
filtered_annotations <- filtered_annotations[filtered_annotations$Localisation != "extra-thoracic",]
### keep only samples that are also in expression data
filtered_annotations <- filtered_annotations[filtered_annotations$sample %in% colnames(expression),]
### filter expression data columns for samples that are in the filtered annotations rows
expression <- expression[,colnames(expression) %in% rownames(filtered_annotations)]


############################################
################ TF ACTIVITY ###############
############################################

### TCGA reference for regulons
data("dorothea_hs_pancancer", package = "dorothea")
regulons_tcga = dorothea_hs_pancancer %>%
  filter(confidence %in% c("A", "B"))

tf_activities_tcga <- run_viper(expression, regulons_tcga, 
                                options =  list(method = "scale", minsize = 4,eset.filter = FALSE, cores = 1, verbose = FALSE))
#write.csv(tf_activities_tcga,"A:/Downloads/Projects/workFromHome/Projects/LungPredict/TF_analysis/TF_activities_TCGA.csv")

### GTEX reference for regulons
data("dorothea_hs", package = "dorothea")
regulons_gtex = dorothea_hs %>%
  filter(confidence %in% c("A", "B"))

tf_activities_gtex <- run_viper(expression, regulons_gtex, 
                                options =  list(method = "scale", minsize = 4,eset.filter = FALSE, cores = 1, verbose = FALSE))

#write.csv(tf_activities_gtex,"A:/Downloads/Projects/workFromHome/Projects/LungPredict/TF_analysis/TF_activities_GTEX.csv")


################################################################################
########## Comparison between results from TCGA and GTEX references  ###########
################################################################################
tf_activities_tcga_df <- as.data.frame(tf_activities_tcga)
tf_activities_gtex_df <- as.data.frame(tf_activities_gtex)
print(rownames(tf_activities_gtex_df))

## filter tf_activities_gtex_df as it has more TF estimated than tf_activities_tcga_df
tf_activities_gtex_df_filtered <- tf_activities_gtex_df[row.names(tf_activities_gtex_df) %in% row.names(tf_activities_tcga_df),]   

tf_activities_tcga_df_filtered <- tf_activities_tcga_df[row.names(tf_activities_tcga_df) %in% row.names(tf_activities_gtex_df_filtered),]   


### compute TF correlations by transposing the matrices (remove the transpose to get patient correlation)
tcga_gtex.patientcor <- cor(tf_activities_tcga_df_filtered, tf_activities_gtex_df_filtered)
?cor
heatmap(tcga_gtex.patientcor)
tcga_gtex.TFcor <- cor(t(tf_activities_tcga_df_filtered), t(tf_activities_gtex_df_filtered))

?heatmap
### compute value differences per patient for each estimated TF
diff <- tf_activities_tcga_df - tf_activities_gtex_df_filtered
diff2 <- diff[rownames(diff) != "SOX2",] %>% as.matrix()
Heatmap(diff2,
         row_names_gp = gpar(fontsize = 7),
         column_names_gp = gpar(fontsize =7),
         )
################################################################################
################################################################################
################################################################################

# drop the sample column for plotting
filtered_annotations <- subset(filtered_annotations, select = -c(sample) )


my_matrix <- as.matrix(tf_activities_tcga)



greyscale = grey.colors(10, rev=T)

# ann_colors = list(
#   Stage = c("IA-1" = greyscale[1], "IA-2" = greyscale[2], "IA-3" = greyscale[3],
#             "IB" = greyscale[4], "IIA"  = greyscale[5], "IIB" = greyscale[6],
#             "IIIA" = greyscale[7], "IIIB" = greyscale[8],"IV-A" = greyscale[9],
#             "IV-B" = greyscale[10]),
#   Mutation = c("EGFR" = "red", "No" = "white", "STK11" = "brown", "TP53" = "purple",
#                "KRAS.G12V" = "deepskyblue", "KRAS.G12D" = "deepskyblue3",
#                "KRAS.G13C" = "deepskyblue4", "KRAS.Q61E" = "blue", "KRAS" = "darkblue",
#                "KRAS.G12C" = "deepskyblue2", "NRAS"="orange", "ERBB2"="pink",
#                "CDKN2A"="darkturquoise",  "MET"="grey", "PTEN"= "aquamarine",
#                "FGFR3"="chocolate", "BRAF"="yellow", "SMAD4"="tan",  "p.E542K"="snow2"),
#   Surgery = c("No" = "darkmagenta", "Yes" = "darkgoldenrod1"),
#   Metastatic = c("Yes"="red", "No"="white"),
#   Sexe = c("M" = "blue", "F" = "pink"),
#   Smoker = c("No"="black", "Yes"="red", "Old"="grey", "Unknown"="white"),
#   Location = c("Lung" = "salmon2", "extra-thoracic" = "seagreen"),
#   Age = c("40-45" = greyscale[1], "45-50" = greyscale[2], "50-55" = greyscale[3],
#   "55-60" = greyscale[4], "60-65"  = greyscale[5], "65-70" = greyscale[6],
#   "70-75" = greyscale[7], "75-80" = greyscale[8],"80-85" = greyscale[9])
# ) 

?rowAnnotation
ann_colors = HeatmapAnnotation( Stage = filtered_annotations$Stage, 
                            Mutation = filtered_annotations$Mutation, 
                            Surgery = filtered_annotations$Surgery,
                            Metastatic = filtered_annotations$Metastatic, 
                            Sexe = filtered_annotations$Sexe, 
                            Smoker = filtered_annotations$Smoker,
                            Location = filtered_annotations$Localisation, 
                            Age = filtered_annotations$Age, 
                            col = list(
                              Stage = c("IA-1" = greyscale[1], "IA-2" = greyscale[2], "IA-3" = greyscale[3],
                                        "IB" = greyscale[4], "IIA"  = greyscale[5], "IIB" = greyscale[6],
                                        "IIIA" = greyscale[7], "IIIB" = greyscale[8],"IV-A" = greyscale[9],
                                        "IV-B" = greyscale[10]),
                              Mutation = c("EGFR" = "red", "No" = "white", "STK11" = "brown", "TP53" = "purple",
                                           "KRAS.G12V" = "deepskyblue", "KRAS.G12D" = "deepskyblue3",
                                           "KRAS.G13C" = "deepskyblue4", "KRAS.Q61E" = "blue", "KRAS" = "darkblue",
                                           "KRAS.G12C" = "deepskyblue2", "NRAS"="orange", "ERBB2"="pink",
                                           "CDKN2A"="darkturquoise",  "MET"="grey", "PTEN"= "aquamarine",
                                           "FGFR3"="chocolate", "BRAF"="yellow", "SMAD4"="tan",  "p.E542K"="snow2"),
                              Surgery = c("No" = "darkmagenta", "Yes" = "darkgoldenrod1"),
                              Metastatic = c("Yes"="red", "No"="white"),
                              Sexe = c("M" = "blue", "F" = "pink"),
                              Smoker = c("No"="black", "Yes"="red", "Old"="salmon", "Unknown"="white"),
                              Location = c("Lung" = "salmon2", "extra-thoracic" = "seagreen"),
                              Age = c("40-45" = greyscale[1], "45-50" = greyscale[2], "50-55" = greyscale[3],
                                      "55-60" = greyscale[4], "60-65"  = greyscale[5], "65-70" = greyscale[6],
                                      "70-75" = greyscale[7], "75-80" = greyscale[8],"80-85" = greyscale[9])
                            ) )
# pheatmap(my_matrix, 
#          main="LP Patient Clustering based on TF activities",
#          annotation_col = filtered_annotations, 
#          annotation_colors = ann_colors,
#          show_rownames = TRUE, 
#          #clustering_distance_rows="euclidean",
#          clustering_method = "complete",
#          scale = "column",
#          cellwidth = 5,
#          cellheight = 5,
#          fontsize_col = 5,
#          fontsize_row = 4)

?Heatmap
Heatmap(my_matrix,
        top_annotation = ann_colors, 
        name = "Patient x TF activities",
        show_row_names = T, 
        show_heatmap_legend = T, 
        row_names_gp = gpar(fontsize = 7),
        column_names_gp = gpar(fontsize =7),
        #          cellwidth = 5,
        #          cellheight = 5,
        #          fontsize_col = 5,
        #          fontsize_row = 4
)


cor.matrix <- cor(my_matrix)

Heatmap(cor.matrix, top_annotation = ann_colors, 
        name = "patient correlation based on TF activities",
        show_row_names = T, show_heatmap_legend = T, 
        row_names_gp = gpar(fontsize = 7),
        column_names_gp = gpar(fontsize =7),
)


cor.matrix.expression <- cor(expression)

# pheatmap(cor.matrix.expression, 
#          main="LP Patient Correlation based on gene expression",
#          annotation_col = filtered_annotations, 
#          annotation_colors = ann_colors,
#          show_rownames = TRUE, 
#          clustering_method = "complete",
#          cellwidth = 5,
#          cellheight = 5,
#          fontsize_col = 5,
#          fontsize_row = 4)

Heatmap(cor.matrix.expression, top_annotation = ann_colors, 
        name = "patient correlation based on expression",
        show_row_names = T, show_heatmap_legend = T, 
        row_names_gp = gpar(fontsize = 7),
        column_names_gp = gpar(fontsize =7),
)



t_matrix <- t(my_matrix)
cor.t_matrix <- cor(t_matrix)
Heatmap(cor.t_matrix, 
        name = "TF correlation based on 64 patients",
        show_row_names = T, show_heatmap_legend = T, 
        row_names_gp = gpar(fontsize = 7),
        column_names_gp = gpar(fontsize =7),
)



#### PCA
res.pca <- prcomp(expression, scale = TRUE)
install.packages("factoextra")
library(factoextra)
fviz_eig(res.pca)
fviz_pca_ind(res.pca,
             col.ind = "cos2", # Color by the quality of representation
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)

library(corbetw2mat)



## Ting
ann_col <- data.frame(as.character(resCut), FFPE_Clinical[Adenocarcinoma,]$Surgery, FFPE_Clinical[Adenocarcinoma,]$diagnoostic_VI_metastasis, FFPE_Clinical[Adenocarcinoma,]$mutation, FFPE_Clinical[Adenocarcinoma,]$Stage, FFPE_Clinical[Adenocarcinoma,]$sexe, FFPE_Clinical[Adenocarcinoma,]$smoker,FFPE_Clinical[Adenocarcinoma,]$localisation_ponction)
rownames(ann_col) <- names(resCut)
colnames(ann_col)=c('Hclust', 'Surgery', 'Metastatic', 'Mutation', 'Stage', 'Sexe', 'Smoker')
greyscale = grey.colors(10, rev=T)
ann_colors = list(
  Hclust = c("1" = "deeppink", "2" = "deepskyblue"),
  Surgery = c("No" = "darkmagenta", "Yes" = "darkgoldenrod1"),
  Metastatic = c("Yes"="red", "No"="blue"),
  Mutation = c("EGFR" = "red", "No" = "white", "STK11" = "brown", "TP53" = "purple",
               "KRAS.G12V" = "deepskyblue", "KRAS.G12D" = "deepskyblue3",
               "KRAS.G13C" = "deepskyblue4", "KRAS.Q61E" = "blue", "KRAS" = "darkblue",
               "KRAS.G12C" = "deepskyblue2", "NRAS"="orange", "ERBB2"="pink",
               "CDKN2A"="darkturquoise",  "MET"="grey", "PTEN"= "aquamarine",
               "FGFR3"="chocolate", "BRAF"="yellow", "SMAD4"="tan",  "p.E542K"="snow2"),
  Stage = c("IA-1" = greyscale[1], "IA-2" = greyscale[2], "IA-3" = greyscale[3],
            "IB" = greyscale[4], "IIA"  = greyscale[5], "IIB" = greyscale[6],
            "IIIA" = greyscale[7], "IIIB" = greyscale[8],"IV-A" = greyscale[9],
            "IV-B" = greyscale[10]),
  Sexe = c("M" = "red", "F" = "blue"),
  Smoker = c("No"="black", "Yes"="red", "Old"="grey"),
  Location = c("Lung" = "salmon2", "extra-thoracic" = "seagreen")
) 
cluster_rows = F, annotation_col = ann_row

res <- pheatmap(t(FFPE_STAR_Normalisation[rownames(variable), Adenocarcinoma]), scale = "row", fontsize_col = 10, show_rownames = T,  show_colnames = F, fontsize_row = 10, main = "Clustering based on gene expression")
resCut = cutree(res$tree_row, 2)
resc1 = names(which(resCut == 1))
resc2 = names(which(resCut == 2))

pheatmap(FFPE_STAR_Normalisation[rownames(variable), Adenocarcinoma], scale = "column", fontsize_col = 10, show_rownames = F, show_colnames = T, cluster_rows = F, annotation_col = ann_col, annotation_colors = ann_colors, fontsize_row = 5, main = "Clustering based on gene expression (Adenocarcinoma)")

















#Julien

annotZ <- rowAnnotation(Mutation = ihc1_l1_network$mutation, 
                        col = list(Mutation = c("other" = "black", "none"="white",
                                                "EGFR" = "red", "EGFR et al" = "orange",
                                                "KRAS" ="blue",  "KRAS et al"="lightblue", 
                                                "STK11" ="purple", "STK11 et al"="pink")))

Heatmap(scale(ihc1_l1_network[,2:5]), column_title = "Cell type estimation across IHC (rescaled)",
        show_column_dend = F, right_annotation = annotZ,
        name = "% of cell", show_row_names = T, show_heatmap_legend = T, row_names_gp = gpar(fontsize = 8))


## Jacobo
col_b <- colorRamp2(c(min(type_moffit["Basal", ]), max(type_moffit["Basal", ])), c("white", "chocolate1")) ## quantitative 
col_collison <- c("NA" = "white", "QM" = "orange", "Classical" = "deepskyblue1") ## category
column_ha <- HeatmapAnnotation(
  Serine = type_serine,
  PHGDH = ccle_sym.enz["PHGDH", ],
  CBS = ccle_sym.enz["CBS", ],
  PSAT1 = ccle_sym.enz["PSAT1", ],
  Classical_enr = type_moffit["Classical", ],
  Basal_enr = type_moffit["Basal", ],
  Collison = type_collison,
  col = list(
    Classical_enr = col_c,
    Basal_enr = col_b,
    Collison = col_collison,
    Serine = col_serine,
    PHGDH = col_PHGDH,
    CBS = col_CBS,
    PSAT1 = col_PSAT1
  )
)
Heatmap(ccle_sym.transclass,
        top_annotation = column_ha, # left_annotation = row_ha, #This one is only for when ccle_sym.transclass is Moffit
        cluster_columns = F
)



# ?viper
# ?msviper
# 
# run_msviper <- function(exprs_m, dorothea, use_aracne, ref, treat, ref_name, treat_name, minsize, ges.filter) {
#   # First we need to generate the phenotype table (AnnotatedDataFrame)
#   conditions <- rep("NA", ncol(exprs_m))
#   conditions[ref] <- ref_name
#   conditions[treat] <- treat_name
#   names(conditions) <- colnames(exprs_m)
#   conditions <- conditions[which(conditions != "NA")]
#   
#   phenotype <- data.frame(condition = factor(conditions))
#   rownames(phenotype) <- names(conditions)
#   
#   phenoData <- new("AnnotatedDataFrame", data = phenotype)
#   
#   exprs_m <- exprs_m[,which(colnames(exprs_m) %in% rownames(phenotype))] %>% as.matrix()
#   
#   # Create Expression set from phenotyble table and expression matrix
#   dset_viper <- ExpressionSet(assayData = exprs_m, phenoData = phenoData)
#   dset_viper$sampleID <- factor(colnames(exprs_m))
#   
#   # Aracne can be used to estimate the mor instead using the -1, 1 from dorothea
#   regulons <- NULL
#   if (use_aracne) {
#     regulons <- dorothea2aracne2viper_regulons(dorothea, dset_viper)
#   } else {
#     regulons <- dorothea2viper_regulons(dorothea)
#   }
#   
#   # We need to create the statistics signature from the conditions
#   signature <- rowTtest(dset_viper, "condition", treat_name, ref_name)
#   statistics_signature <- (qnorm(signature$p.value / 2, lower.tail = FALSE) * sign(signature$statistic))[, 1]
#   # Generate the null model with bootstrapping (1000 iterations)
#   nullmodel <- ttestNull(dset_viper, "condition", treat_name, ref_name, per = 1000, repos = T, verbose = F)
#   # Run msviper using the statistics signature, the regulons converted from dorothea table, the null model the minSize of regulon and the ges.filter
#   mrs <- msviper(ges = statistics_signature, regulon = regulons, nullmodel = nullmodel, minsize = minsize, ges.filter = ges.filter, verbose = F)
#   # Convert the msviper regulons to dorothea
#   dorothea_mrs_regulons <- viper_regulons2dorothea(mrs$regulon) %>%
#     mutate(state = ifelse(mor > 0, "activation", "inhibition"))
#   # Generate a table with the TFs, the regulon size, the NES score, the pval and the pval.fdr
#   mrs_table <- tibble(TF = names(mrs$es$p.value), size = mrs$es$size, nes = mrs$es$nes, pval = mrs$es$p.value, pval.fdr = p.adjust(mrs$es$p.value, method = "fdr")) %>% arrange(pval)
#   
#   list(mrs_table = mrs_table, mrs = mrs, regulons = dorothea_mrs_regulons)
# }

TF_activities_analysis.R 
#required packages
# install.packages("pillar")
# install.packages("ggplot2")
# install.packages("dplyr")
# install.packages("forcats")
# BiocManager::install("viper")
#BiocManager::install("dorothea")
# ??dorothea
#install.packages("devtools")
# library(devtools)
# install_github("saezlab/dorothea") # to get the TCGA (pan cancer) ref for regulons

library("tidyverse")
library("viper")
library("dorothea")
library(pheatmap)
library(RColorBrewer)
library(viridis)
library(ComplexHeatmap)

#################################
######## EXPRESSION DATA ########
#################################

#loading expression file
expression=read.csv("A:/Downloads/Projects/workFromHome/Projects/LungPredict/TF_analysis/LP_FFPE_STAR_RSEM_TPM.txt", sep = "	")
expression=expression%>% column_to_rownames(var='Gene') %>% as.matrix()

#print(colnames(expression))
### remove patient 69 from expression data
# expression <- subset(expression, select = -c(LP.01.69))

#########################################################
###### FILTER PATIENT BASED ON CLINICAL FEATURES  #######
#########################################################

### create annotations clinical feature file
clinical_features <- read.table("A:/Downloads/Projects/workFromHome/Projects/LungPredict/TF_analysis/clinic_data_v2_clean.csv", header = TRUE, sep = ",")

# look for specific mutations
egfr_mutant <- clinical_features$sample[clinical_features$mutation == "EGFR"]
kras_mutant <- clinical_features$sample[grep("KRAS", clinical_features$mutation)]
stk11_mutant <- clinical_features$sample[clinical_features$mutation == "STK11"]
`%ni%` <- Negate(`%in%`)
clinical_features <- clinical_features %>% 
  mutate(KRAS = case_when(
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
  ))


# add age category
clinical_features <- clinical_features %>%
  mutate(age_category = case_when(
    age < 85 & age >= 80 ~ '80-85',
    age < 80 & age >= 75 ~ '75-80',
    age < 75 & age >= 70 ~ '70-75',
    age < 70 & age >= 65 ~ '65-70',
    age < 65 & age >= 60 ~ '60-65',
    age < 60 & age >= 55 ~ '55-60',
    age < 55 & age >= 50 ~ '50-55',
    age < 50 & age >= 45 ~ '45-50',
    age < 45 & age >= 40 ~ '40-45'
  ))


# expBasedClusters <- read.table("A:/Downloads/Projects/workFromHome/Projects/LungPredict/TF_analysis/expression_based_patient_clusters2.txt", header = TRUE, sep = "	")
expBasedClusters <- read.table("A:/Downloads/Projects/workFromHome/Projects/LungPredict/TF_analysis/a.txt", header = TRUE, sep = "	")
colnames(expBasedClusters)[1] <- "sample"

reactomeClusters <- read.table("A:/Downloads/Projects/workFromHome/Projects/LungPredict/TF_analysis/ReactomeClustersAllPatients.csv", header = TRUE, sep = ",")
colnames(reactomeClusters)[1] <- "sample"

deconvCancerClusters <- read.table("A:/Downloads/Projects/workFromHome/Projects/LungPredict/TF_analysis/DeconvCancerClusters.txt", header = TRUE, sep = "	")
colnames(deconvCancerClusters)[1] <- "sample"




annotations <- data.frame(clinical_features$sample, 
                          clinical_features$Stage, 
                          clinical_features$Surgery, 
                          clinical_features$diagnostic,
                          clinical_features$diagnoostic_VI_metastasis, 
                          clinical_features$sexe, 
                          clinical_features$smoker, 
                          clinical_features$localisation_ponction,
                          clinical_features$age_category,
                          clinical_features$KRAS,
                          clinical_features$EGFR,
                          clinical_features$STK11)

# create the annotation DF
rownames(annotations) <- annotations$clinical_features.sample


colnames(annotations) <- c("sample", "Stage", "Surgery", "diagnostic", "Metastatic", "Sexe", "Smoker", "Localisation", 
                           "Age", "KRAS", "EGFR", "STK11")

############################################################################  
################ REMOVE UNDESIRED PATIENTS ################################  
################ AND FILTER THE ANNOTATIONS AND EXPRESSION DATA FOR OUR DATA
############################################################################ 

### remove patient 56B from annotation data
filtered_annotations <- annotations[annotations$sample != 'LP.01.56B',]
### remove patient 69 from annotation data
filtered_annotations <- filtered_annotations[filtered_annotations$sample != 'LP.01.69',]
### remove extra-thoracic samples from annotation data
filtered_annotations <- filtered_annotations[filtered_annotations$Localisation != "extra-thoracic",]
### keep only the adenocarcinoma ??
#### filtered_annotations <- filtered_annotations[filtered_annotations$diagnostic == "Adenocarcinoma",]


#Reduce(function(x, y) merge(x, y, by = "sample"), list(df1, df2, df3))

?right_join

outputAnnotations <-  right_join(expBasedClusters, filtered_annotations, by = "sample")
colnames(outputAnnotations)[2] <- "GE_Cluster"


outputAnnotations <-  right_join(reactomeClusters, outputAnnotations, by = "sample")
colnames(outputAnnotations)[2] <- "reactomeCluster"

outputAnnotations2 <-  right_join(deconvCancerClusters[,1:2], outputAnnotations, by = "sample")
outputAnnotations3 <-  right_join(deconvCancerClusters[,c(1,3)], outputAnnotations2, by = "sample")

outputAnnotations4 <- outputAnnotations3[,-12]



write.table(outputAnnotations4, file= "A:/Downloads/Projects/workFromHome/Projects/LungPredict/TF_analysis/full_annotations_with_clusters_corrected1.txt", sep="	", quote=F, row.names=FALSE)






### keep only samples that are also in expression data
filtered_annotations <- filtered_annotations[filtered_annotations$sample %in% colnames(expression),]



### filter expression data columns for samples that are in the filtered annotations rows
expression <- expression[,colnames(expression) %in% rownames(filtered_annotations)]


# create the annotation
greyscale <- grey.colors(10, rev = T)
ann_colors <- HeatmapAnnotation(
  Stage = filtered_annotations$Stage, Surgery = filtered_annotations$Surgery,
  Age = filtered_annotations$Age,
  Metastatic = filtered_annotations$Metastatic, 
  Sexe = filtered_annotations$Sexe, 
  Smoker = filtered_annotations$Smoker,
  Location = filtered_annotations$Localisation, 
  Age = filtered_annotations$Age, 
  KRAS = filtered_annotations$KRAS,
  EGFR = filtered_annotations$EGFR,
  STK11 = filtered_annotations$STK11,
  col = list(
    Stage = c(
      "IA-1" = greyscale[1], "IA-2" = greyscale[2], "IA-3" = greyscale[3],
      "IB" = greyscale[4], "IIA" = greyscale[5], "IIB" = greyscale[6],
      "IIIA" = greyscale[7], "IIIB" = greyscale[8], "IV-A" = greyscale[9],
      "IV-B" = greyscale[10]
    ),
    Surgery = c("No" = "darkmagenta", "Yes" = "darkgoldenrod1"),
    Metastatic = c("Yes" = "red", "No" = "white"),
    Sexe = c("M" = "blue", "F" = "pink"),
    Smoker = c("No" = "black", "Yes" = "red", "Old" = "salmon", "Unknown" = "white"),
    Location = c("Lung" = "salmon2", "extra-thoracic" = "seagreen"),
    Age = c(
      "40-45" = greyscale[1], "45-50" = greyscale[2], "50-55" = greyscale[3],
      "55-60" = greyscale[4], "60-65" = greyscale[5], "65-70" = greyscale[6],
      "70-75" = greyscale[7], "75-80" = greyscale[8], "80-85" = greyscale[9]
    ),
    KRAS = c("yes" = "blue", "no" = "white"),
    EGFR = c("yes" = "red", "no" = "white"),
    STK11 = c("yes" = "green", "no" = "white")
  )
)

############################################
## GENE EXPRESSION BASED CLUSTERING ##
############################################

cor.matrix.expression <- cor(expression)

# myPheatmap <- pheatmap(cor.matrix.expression,
#          main="LP Patient Correlation based on gene expression",
#          annotation_col = filtered_annotations,
#          annotation_colors = ann_colors,
#          show_rownames = TRUE,
#          clustering_method = "complete",
#          cellwidth = 5,
#          cellheight = 5,
#          fontsize_col = 5,
#          fontsize_row = 4)

# drop the sample column for plotting
filtered_annotations <- subset(filtered_annotations, select = -c(sample) )

myExpBasedPatientCorrHeatmap <- Heatmap(cor.matrix.expression, top_annotation = ann_colors, 
        name = "patient correlation based on expression",
        show_row_names = T, show_heatmap_legend = T, 
        row_names_gp = gpar(fontsize = 8),
        column_names_gp = gpar(fontsize = 8)
)

print(myExpBasedPatientCorrHeatmap)

### Extract clusters from the heatmap
HM = Heatmap(cor.matrix.expression, top_annotation = ann_colors, 
             name = "patient correlation based on expression",
             show_row_names = T, show_heatmap_legend = T, 
             row_names_gp = gpar(fontsize = 8),
             column_names_gp = gpar(fontsize = 8), 
             km = 2)
?Heatmap

hc <- hclust(dist(cor.matrix.expression)) 
resCut <- cutree(hc, k = 3)
resCut.dataframe <- as.data.frame(resCut)
# resCut.dataframe$sample <- rownames(resCut.dataframe)

write.table(data.frame("sample"=rownames(resCut.dataframe), resCut.dataframe), "A:/Downloads/Projects/workFromHome/Projects/LungPredict/TF_analysis/a.txt", row.names=FALSE, quote=F, sep="	")



write.table(resCut.dataframe, file= "A:/Downloads/Projects/workFromHome/Projects/LungPredict/TF_analysis/expression_based_patient_clusters_resCutCorrected.txt", sep="	", quote=F, row.names=TRUE)





HM = draw(HM)
r.dend <- row_dend(HM)  #Extract row dendrogram
rcl.list <- row_order(HM)  #Extract clusters (output is a list)
lapply(rcl.list, function(x) length(x))  #check/confirm size clusters

# loop to extract patient_id for each cluster.
for (i in 1:length(row_order(HM))){
  if (i == 1) {
    clu <- t(t(row.names(cor.matrix.expression[row_order(HM)[[i]],])))
    out <- cbind(clu, paste("Cluster", i, sep=""))
    colnames(out) <- c("PatientID", "Cluster")
    } else {
      clu <- t(t(row.names(cor.matrix.expression[row_order(HM)[[i]],])))
      clu <- cbind(clu, paste("Cluster", i, sep=""))
      out <- rbind(out, clu)
      }
}

#check
out 

#export
write.table(out, file= "A:/Downloads/Projects/workFromHome/Projects/LungPredict/TF_analysis/expression_based_patient_clusters2.txt", sep="	", quote=F, row.names=FALSE)

############################################
################ TF ACTIVITY ###############
############################################

### TCGA reference for regulons
data("dorothea_hs_pancancer", package = "dorothea")
regulons_tcga = dorothea_hs_pancancer %>%
  filter(confidence %in% c("A", "B"))

tf_activities_tcga <- run_viper(expression, regulons_tcga, 
                                options =  list(method = "scale", minsize = 4,eset.filter = FALSE, cores = 1, verbose = FALSE))
#write.csv(tf_activities_tcga,"A:/Downloads/Projects/workFromHome/Projects/LungPredict/TF_analysis/TF_activities_TCGA.csv")

### GTEX reference for regulons
data("dorothea_hs", package = "dorothea")
regulons_gtex = dorothea_hs %>%
  filter(confidence %in% c("A", "B"))

tf_activities_gtex <- run_viper(expression, regulons_gtex, 
                                options =  list(method = "scale", minsize = 4,eset.filter = FALSE, cores = 1, verbose = FALSE))

#write.csv(tf_activities_gtex,"A:/Downloads/Projects/workFromHome/Projects/LungPredict/TF_analysis/TF_activities_GTEX.csv")


################################################################################
########## Comparison between results from TCGA and GTEX references  ###########
################################################################################
tf_activities_tcga_df <- as.data.frame(tf_activities_tcga)
tf_activities_gtex_df <- as.data.frame(tf_activities_gtex)
print(rownames(tf_activities_gtex_df))

## filter tf_activities_gtex_df as it has more TF estimated than tf_activities_tcga_df
tf_activities_gtex_df_filtered <- tf_activities_gtex_df[row.names(tf_activities_gtex_df) %in% row.names(tf_activities_tcga_df),]   

tf_activities_tcga_df_filtered <- tf_activities_tcga_df[row.names(tf_activities_tcga_df) %in% row.names(tf_activities_gtex_df_filtered),]   


### compute TF correlations by transposing the matrices (remove the transpose to get patient correlation)
tcga_gtex.patientcor <- cor(tf_activities_tcga_df_filtered, tf_activities_gtex_df_filtered)
?cor
heatmap(tcga_gtex.patientcor)
tcga_gtex.TFcor <- cor(t(tf_activities_tcga_df_filtered), t(tf_activities_gtex_df_filtered))

?heatmap
### compute value differences per patient for each estimated TF
diff <- tf_activities_tcga_df - tf_activities_gtex_df_filtered
diff2 <- diff[rownames(diff) != "SOX2",] %>% as.matrix()
Heatmap(diff2,
         row_names_gp = gpar(fontsize = 7),
         column_names_gp = gpar(fontsize =7),
         )
################################################################################
################################################################################
################################################################################

# drop the sample column for plotting
filtered_annotations <- subset(filtered_annotations, select = -c(sample) )


my_matrix_TF <- as.matrix(tf_activities_tcga)


################################################################################
######### REMOVE TF THAT DONT VARY A LOT #######################################
################################################################################

TFstd <- apply(my_matrix_TF, 1, sd)
# print(TFstd)
hist(TFstd, breaks=106)

### Remove all TF that have std < 2

my_matrix_TF_filtered <- my_matrix_TF[TFstd > 2,]

################################################################################
##### UPDATE THE ANNOTATIONS WITH EXPRESSION BASED CLUSTERS #####
### AND THE DECONV CLUSTERS
#### AND THE REACTOME CLUSTERS
################################################################################

expBasedClusters <- read.table("A:/Downloads/Projects/workFromHome/Projects/LungPredict/TF_analysis/expression_based_patient_clusters.txt", header = TRUE, sep = "	")
colnames(expBasedClusters)[1] <- "sample"

reactomeClusters <- read.table("A:/Downloads/Projects/workFromHome/Projects/LungPredict/TF_analysis/Reactomes_cluster.txt", header = TRUE, sep = "	")
colnames(expBasedClusters)[1] <- "sample"

mergedAnnotations <- inner_join(expBasedClusters, annotations, by = "sample")
# ?inner_join

ann_colors_merged <- HeatmapAnnotation(
  Stage = mergedAnnotations$Stage, Surgery = mergedAnnotations$Surgery,
  Age = mergedAnnotations$Age,
  Metastatic = mergedAnnotations$Metastatic, 
  Sexe = mergedAnnotations$Sexe, 
  Smoker = mergedAnnotations$Smoker,
  Location = mergedAnnotations$Localisation, 
  Age = mergedAnnotations$Age, 
  KRAS = mergedAnnotations$KRAS,
  EGFR = mergedAnnotations$EGFR,
  STK11 = mergedAnnotations$STK11,
  Cluster = mergedAnnotations$Cluster,
  col = list(
    Stage = c(
      "IA-1" = greyscale[1], "IA-2" = greyscale[2], "IA-3" = greyscale[3],
      "IB" = greyscale[4], "IIA" = greyscale[5], "IIB" = greyscale[6],
      "IIIA" = greyscale[7], "IIIB" = greyscale[8], "IV-A" = greyscale[9],
      "IV-B" = greyscale[10]
    ),
    Surgery = c("No" = "darkmagenta", "Yes" = "darkgoldenrod1"),
    Metastatic = c("Yes" = "red", "No" = "white"),
    Sexe = c("M" = "blue", "F" = "pink"),
    Smoker = c("No" = "black", "Yes" = "red", "Old" = "salmon", "Unknown" = "white"),
    Location = c("Lung" = "salmon2", "extra-thoracic" = "seagreen"),
    Age = c(
      "40-45" = greyscale[1], "45-50" = greyscale[2], "50-55" = greyscale[3],
      "55-60" = greyscale[4], "60-65" = greyscale[5], "65-70" = greyscale[6],
      "70-75" = greyscale[7], "75-80" = greyscale[8], "80-85" = greyscale[9]
    ),
    KRAS = c("yes" = "blue", "no" = "white"),
    EGFR = c("yes" = "red", "no" = "white"),
    STK11 = c("yes" = "green", "no" = "white"),
    Cluster = c("Cluster1" = "pink", "Cluster2" = "blue")
  )
)

Heatmap(my_matrix_TF_filtered,
        top_annotation = ann_colors_merged, 
        name = "Patient x TF activities",
        show_row_names = T, 
        show_heatmap_legend = T, 
        row_names_gp = gpar(fontsize = 7),
        column_names_gp = gpar(fontsize =7),
        #          cellwidth = 5,
        #          cellheight = 5,
        #          fontsize_col = 5,
        #          fontsize_row = 4
)

################################################################################
###### TF-BASED PATIENT CORRELATION MATRIX #####
################################################################################

cor.matrix.TF.filtered <- cor(my_matrix_TF_filtered)

Heatmap(cor.matrix.TF.filtered, top_annotation = ann_colors_merged, 
        name = "patient correlation based on TF activities",
        show_row_names = T, show_heatmap_legend = T, 
        row_names_gp = gpar(fontsize = 7),
        column_names_gp = gpar(fontsize =7),
)




t_matrix <- t(my_matrix)
cor.t_matrix <- cor(t_matrix)
Heatmap(cor.t_matrix, 
        name = "TF correlation based on 64 patients",
        show_row_names = T, show_heatmap_legend = T, 
        row_names_gp = gpar(fontsize = 7),
        column_names_gp = gpar(fontsize =7),
)



#### PCA
res.pca <- prcomp(expression, scale = TRUE)
install.packages("factoextra")
library(factoextra)
fviz_eig(res.pca)
fviz_pca_ind(res.pca,
             col.ind = "cos2", # Color by the quality of representation
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)

library(corbetw2mat)



## Ting
ann_col <- data.frame(as.character(resCut), FFPE_Clinical[Adenocarcinoma,]$Surgery, FFPE_Clinical[Adenocarcinoma,]$diagnoostic_VI_metastasis, FFPE_Clinical[Adenocarcinoma,]$mutation, FFPE_Clinical[Adenocarcinoma,]$Stage, FFPE_Clinical[Adenocarcinoma,]$sexe, FFPE_Clinical[Adenocarcinoma,]$smoker,FFPE_Clinical[Adenocarcinoma,]$localisation_ponction)
rownames(ann_col) <- names(resCut)
colnames(ann_col)=c('Hclust', 'Surgery', 'Metastatic', 'Mutation', 'Stage', 'Sexe', 'Smoker')
greyscale = grey.colors(10, rev=T)
ann_colors = list(
  Hclust = c("1" = "deeppink", "2" = "deepskyblue"),
  Surgery = c("No" = "darkmagenta", "Yes" = "darkgoldenrod1"),
  Metastatic = c("Yes"="red", "No"="blue"),
  Mutation = c("EGFR" = "red", "No" = "white", "STK11" = "brown", "TP53" = "purple",
               "KRAS.G12V" = "deepskyblue", "KRAS.G12D" = "deepskyblue3",
               "KRAS.G13C" = "deepskyblue4", "KRAS.Q61E" = "blue", "KRAS" = "darkblue",
               "KRAS.G12C" = "deepskyblue2", "NRAS"="orange", "ERBB2"="pink",
               "CDKN2A"="darkturquoise",  "MET"="grey", "PTEN"= "aquamarine",
               "FGFR3"="chocolate", "BRAF"="yellow", "SMAD4"="tan",  "p.E542K"="snow2"),
  Stage = c("IA-1" = greyscale[1], "IA-2" = greyscale[2], "IA-3" = greyscale[3],
            "IB" = greyscale[4], "IIA"  = greyscale[5], "IIB" = greyscale[6],
            "IIIA" = greyscale[7], "IIIB" = greyscale[8],"IV-A" = greyscale[9],
            "IV-B" = greyscale[10]),
  Sexe = c("M" = "red", "F" = "blue"),
  Smoker = c("No"="black", "Yes"="red", "Old"="grey"),
  Location = c("Lung" = "salmon2", "extra-thoracic" = "seagreen")
) 
cluster_rows = F, annotation_col = ann_row

res <- pheatmap(t(FFPE_STAR_Normalisation[rownames(variable), Adenocarcinoma]), scale = "row", fontsize_col = 10, show_rownames = T,  show_colnames = F, fontsize_row = 10, main = "Clustering based on gene expression")
resCut = cutree(res$tree_row, 2)
resc1 = names(which(resCut == 1))
resc2 = names(which(resCut == 2))

pheatmap(FFPE_STAR_Normalisation[rownames(variable), Adenocarcinoma], scale = "column", fontsize_col = 10, show_rownames = F, show_colnames = T, cluster_rows = F, annotation_col = ann_col, annotation_colors = ann_colors, fontsize_row = 5, main = "Clustering based on gene expression (Adenocarcinoma)")

















#Julien

annotZ <- rowAnnotation(Mutation = ihc1_l1_network$mutation, 
                        col = list(Mutation = c("other" = "black", "none"="white",
                                                "EGFR" = "red", "EGFR et al" = "orange",
                                                "KRAS" ="blue",  "KRAS et al"="lightblue", 
                                                "STK11" ="purple", "STK11 et al"="pink")))

Heatmap(scale(ihc1_l1_network[,2:5]), column_title = "Cell type estimation across IHC (rescaled)",
        show_column_dend = F, right_annotation = annotZ,
        name = "% of cell", show_row_names = T, show_heatmap_legend = T, row_names_gp = gpar(fontsize = 8))


## Jacobo
col_b <- colorRamp2(c(min(type_moffit["Basal", ]), max(type_moffit["Basal", ])), c("white", "chocolate1")) ## quantitative 
col_collison <- c("NA" = "white", "QM" = "orange", "Classical" = "deepskyblue1") ## category
column_ha <- HeatmapAnnotation(
  Serine = type_serine,
  PHGDH = ccle_sym.enz["PHGDH", ],
  CBS = ccle_sym.enz["CBS", ],
  PSAT1 = ccle_sym.enz["PSAT1", ],
  Classical_enr = type_moffit["Classical", ],
  Basal_enr = type_moffit["Basal", ],
  Collison = type_collison,
  col = list(
    Classical_enr = col_c,
    Basal_enr = col_b,
    Collison = col_collison,
    Serine = col_serine,
    PHGDH = col_PHGDH,
    CBS = col_CBS,
    PSAT1 = col_PSAT1
  )
)
Heatmap(ccle_sym.transclass,
        top_annotation = column_ha, # left_annotation = row_ha, #This one is only for when ccle_sym.transclass is Moffit
        cluster_columns = F
)



# ?viper
# ?msviper
# 
# run_msviper <- function(exprs_m, dorothea, use_aracne, ref, treat, ref_name, treat_name, minsize, ges.filter) {
#   # First we need to generate the phenotype table (AnnotatedDataFrame)
#   conditions <- rep("NA", ncol(exprs_m))
#   conditions[ref] <- ref_name
#   conditions[treat] <- treat_name
#   names(conditions) <- colnames(exprs_m)
#   conditions <- conditions[which(conditions != "NA")]
#   
#   phenotype <- data.frame(condition = factor(conditions))
#   rownames(phenotype) <- names(conditions)
#   
#   phenoData <- new("AnnotatedDataFrame", data = phenotype)
#   
#   exprs_m <- exprs_m[,which(colnames(exprs_m) %in% rownames(phenotype))] %>% as.matrix()
#   
#   # Create Expression set from phenotyble table and expression matrix
#   dset_viper <- ExpressionSet(assayData = exprs_m, phenoData = phenoData)
#   dset_viper$sampleID <- factor(colnames(exprs_m))
#   
#   # Aracne can be used to estimate the mor instead using the -1, 1 from dorothea
#   regulons <- NULL
#   if (use_aracne) {
#     regulons <- dorothea2aracne2viper_regulons(dorothea, dset_viper)
#   } else {
#     regulons <- dorothea2viper_regulons(dorothea)
#   }
#   
#   # We need to create the statistics signature from the conditions
#   signature <- rowTtest(dset_viper, "condition", treat_name, ref_name)
#   statistics_signature <- (qnorm(signature$p.value / 2, lower.tail = FALSE) * sign(signature$statistic))[, 1]
#   # Generate the null model with bootstrapping (1000 iterations)
#   nullmodel <- ttestNull(dset_viper, "condition", treat_name, ref_name, per = 1000, repos = T, verbose = F)
#   # Run msviper using the statistics signature, the regulons converted from dorothea table, the null model the minSize of regulon and the ges.filter
#   mrs <- msviper(ges = statistics_signature, regulon = regulons, nullmodel = nullmodel, minsize = minsize, ges.filter = ges.filter, verbose = F)
#   # Convert the msviper regulons to dorothea
#   dorothea_mrs_regulons <- viper_regulons2dorothea(mrs$regulon) %>%
#     mutate(state = ifelse(mor > 0, "activation", "inhibition"))
#   # Generate a table with the TFs, the regulon size, the NES score, the pval and the pval.fdr
#   mrs_table <- tibble(TF = names(mrs$es$p.value), size = mrs$es$size, nes = mrs$es$nes, pval = mrs$es$p.value, pval.fdr = p.adjust(mrs$es$p.value, method = "fdr")) %>% arrange(pval)
#   
#   list(mrs_table = mrs_table, mrs = mrs, regulons = dorothea_mrs_regulons)
# }
