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
import networkx as nx
from igraph import *



df1 = pd.read_csv("result6_fillna_control_renamed_filtered4_exp130921_T0_renamed_replicates_averages.csv", sep=';', header=0, index_col=0)

print (df1) #
df1 = df1.astype(float)
df1 = df1.sort_index(axis=1, ascending=True, key=lambda x: x.str.lower())

matrix1 = df1.corr(
    method = 'pearson',  # The method of correlation
    min_periods = 1).round(2)     # Min number of observations required)


 
# Transform it in a links data frame (3 columns only):
links = matrix1.stack().reset_index()
links.columns = ['source', 'target', 'correlation']
 



# Keep only correlation over a threshold and remove self correlation (cor(A,A)=1)
links_filtered=links.loc[ (links['correlation'] > 0.75) & (links['source'] != links['target']) ]
# Build your graph
G = Graph.DataFrame(edges, directed=False, vertices=links['correlation'])




# plt.show()








stop = timeit.default_timer()
print(stop - start) 