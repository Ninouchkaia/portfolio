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
 

df = pd.read_csv("fold_changes_agreggated.csv", sep=';', header=0, index_col=0)
print (df) #

matrix = df.corr(
    method = 'pearson',  # The method of correlation
    min_periods = 1).round(2)     # Min number of observations required)

matrix.to_csv('fc_correlations.csv', sep='\t', index=True)

# ##################  DENDROGRAM ########################################
# plt.figure(figsize=(12,5))
# dissimilarity = 1 - abs(matrix)
# Z = linkage(squareform(dissimilarity), 'complete')

# with plt.rc_context({'lines.linewidth': 0.5}):
# 	dendrogram(Z, labels=df.columns, orientation='top', leaf_rotation=90, leaf_font_size = 1)

# plt.savefig("dendro.pdf", format='pdf')
# # plt.show()
# ########################################################################

# ############# CORRPLOT + DENDROGRAM ####################################
# # sns.clustermap(matrix, method="complete", cmap='RdBu', annot=True, 
# #                annot_kws={"size": 7}, vmin=-1, vmax=1, figsize=(15,12));
# # plt.show()
# ########################################################################

stop = timeit.default_timer()
print(stop - start) 