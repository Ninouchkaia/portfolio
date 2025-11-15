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

df1 = pd.read_csv("result6_fillna_control_renamed_filtered4_exp130921_T0_renamed_replicates_averages.csv", sep=';', header=0, index_col=0)
print (df1) #
df1 = df1.astype(float)
df1 = df1.sort_index(axis=1, ascending=True, key=lambda x: x.str.lower())

matrix1 = df1.corr(
    method = 'pearson',  # The method of correlation
    min_periods = 1).round(2)     # Min number of observations required)

# m1 = sns.heatmap(matrix1,annot=False, yticklabels=True, xticklabels=True)
# m1.set_yticklabels(m1.get_ymajorticklabels(), size = 2)
# m1.set_xticklabels(m1.get_xmajorticklabels(), size = 2)
# m1.tick_params(width=0.5 )#left=False, bottom=False)
# m1.set_title('avg_reads_correlations')

# plt.savefig("avg_reads_correlations.pdf", format='pdf')

df2 = pd.read_csv("result6_fillna_control_renamed_filtered4_exp130921_T0_renamed_replicates_averages_fc_zeros.csv", sep=';', header=0, index_col=0)
print (df2) #
df2 = df2.astype(float)
df2 = df2.sort_index(axis=1, ascending=True, key=lambda x: x.str.lower())

matrix2 = df2.corr(
    method = 'pearson',  # The method of correlation
    min_periods = 1).round(2)     # Min number of observations required)

# m2 = sns.heatmap(matrix2,annot=False, yticklabels=True, xticklabels=True)
# m2.set_yticklabels(m2.get_ymajorticklabels(), size = 2)
# m2.set_xticklabels(m2.get_xmajorticklabels(), size = 2)
# m2.tick_params(width=0.5 )#left=False, bottom=False)
# m2.set_title('fold_changes_correlations')
# plt.savefig("fold_changes_correlations.pdf", format='pdf')


diff = matrix1 - matrix2
h = sns.heatmap(diff,annot=False, yticklabels=True, xticklabels=True)
# h.set_yticklabels(h.get_ymajorticklabels(), size = 2)
# h.set_xticklabels(h.get_xmajorticklabels(), size = 2)
# h.tick_params(width=0.5 )#left=False, bottom=False)
# h.set_title('difference between correlations based on avg_reads vs. fold changes')
# plt.savefig("diff_correl.pdf", format='pdf')


clustered_diff = sns.clustermap(diff, method="complete", cmap='RdBu', annot=False, yticklabels=True, xticklabels=True,
               annot_kws={"size": 7}, vmin=-1, vmax=1, figsize=(15,12));

clustered_diff.ax_heatmap.set_xticklabels(clustered_diff.ax_heatmap.get_xmajorticklabels(), fontsize = 6)
clustered_diff.ax_heatmap.set_yticklabels(clustered_diff.ax_heatmap.get_ymajorticklabels(), fontsize = 6)
clustered_diff.fig.suptitle('"clustered_diff between correlations based on avg_reads vs. fold changes"') 


plt.savefig("clustered_diff.pdf", format='pdf')

plt.show()





# plt.show()

# matrix.to_csv('correl_result6_fillna_control_renamed_filtered4_exp130921_T0_renamed_replicates_averages.csv', sep=';', index=True)


##################  DENDROGRAM ########################################
# plt.figure(figsize=(12,5))
# dissimilarity = 1 - abs(matrix)
# Z = linkage(squareform(dissimilarity), 'complete')

# with plt.rc_context({'lines.linewidth': 0.5}):
# 	dendrogram(Z, labels=df.columns, orientation='top', leaf_rotation=90, leaf_font_size = 1)

# plt.savefig("dendro.pdf", format='pdf')
# # plt.show()
########################################################################

############# CORRPLOT + DENDROGRAM ####################################
# g = sns.clustermap(matrix, method="complete", cmap='RdBu', annot=False, yticklabels=True, xticklabels=True,
#                annot_kws={"size": 7}, vmin=-1, vmax=1, figsize=(15,12));

# g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xmajorticklabels(), fontsize = 6)
# g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_ymajorticklabels(), fontsize = 6)

# plt.savefig("corr_dendro_result6_fillna_control_renamed_filtered4_exp130921_T0_renamed_replicates_averages_fc_zeros.pdf", format='pdf')

# plt.show()
########################################################################

stop = timeit.default_timer()
print(stop - start) 