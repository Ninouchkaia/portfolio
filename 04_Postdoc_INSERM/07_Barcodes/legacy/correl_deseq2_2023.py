#!/usr/bin/python
# -*- coding: utf-8 -*-
# jupyter notebook --browser='C:/Program Files (x86)/Google/Chrome/Application/chrome.exe %s'

import timeit
start = timeit.default_timer()

from scipy.cluster.hierarchy import linkage, dendrogram, fcluster
from scipy.spatial.distance import squareform


import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
 

# df = pd.read_csv("result6_fillna_control_renamed_filtered6.csv", sep=';', header=0, index_col=0)

df1 = pd.read_csv("merged_logfc_pval_filtered_deseq2_2023.csv", sep=';', header=0, index_col=0)
print (df1) #
df1 = df1.astype(float)
# df1 = df1.dropna(inplace=True)
df1 = df1.sort_index(axis=1, ascending=True, key=lambda x: x.str.lower())

df1 = df1.fillna('NA')
df1.drop(columns=df1.columns[0], axis=1,  inplace=True)
print(df1)
df1.to_csv('merged_logfc_pval_filtered_deseq2_2023_fillna.csv', sep=';', index=True)


matrix1 = df1.corr(
    method = 'pearson',  # The method of correlation
    min_periods = 1).round(2)     # Min number of observations required)

print(matrix1)

# matrix1 = matrix1.fillna(0)
# print(matrix1)

# matrix1 = matrix1.sort_index(axis=1, ascending=True, key=lambda x: x.str.lower())


# matrix1.to_csv('correl_matrix1_merged_pval_filtered_deseq2_fillna.csv', sep=';', index=True)

# m1 = sns.heatmap(matrix1,annot=False, yticklabels=True, xticklabels=True)
# m1.set_yticklabels(m1.get_ymajorticklabels(), size = 2)
# m1.set_xticklabels(m1.get_xmajorticklabels(), size = 2)
# m1.tick_params(width=0.5 )#left=False, bottom=False)
# # m1.set_title('deseq2_logfc_correlations')

# plt.savefig("deseq2_logfc_correlations_pval_filtered_ordered_2023.pdf", format='pdf')


# ##################  DENDROGRAM ########################################
# plt.figure(figsize=(12,5))
# dissimilarity = 1 - abs(matrix1)
# Z = linkage(squareform(dissimilarity), 'complete')

# with plt.rc_context({'lines.linewidth': 0.5}):
# 	dendrogram(Z, labels=df.columns, orientation='top', leaf_rotation=90, leaf_font_size = 1)

# plt.savefig("dendro.pdf", format='pdf')
# plt.show()
# ########################################################################

########### CORRPLOT + DENDROGRAM ####################################
# g = sns.clustermap(matrix1, method="complete", cmap='RdBu', annot=False, yticklabels=True, xticklabels=True,
#                annot_kws={"size": 7}, vmin=-1, vmax=1, figsize=(15,12));

# g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xmajorticklabels(), fontsize = 6)
# g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_ymajorticklabels(), fontsize = 6)
# g.ax_row_dendrogram.set_visible(False)



# plt.savefig("deseq2_pvalfiltered_logfc_correlations_clustered.pdf", format='pdf')

# plt.show()
#######################################################################
# import scipy.spatial as sp, scipy.cluster.hierarchy as hc

# row_dism = 1 - df1.T.corr()
# row_linkage = hc.linkage(sp.distance.squareform(row_dism), method='complete')
# col_dism = 1 - df1.corr()
# col_linkage = hc.linkage(sp.distance.squareform(col_dism), method='complete')

# sns.clustermap(df1,figsize=(5, 5),row_linkage=row_linkage, col_linkage=col_linkage)

stop = timeit.default_timer()
print(stop - start) 