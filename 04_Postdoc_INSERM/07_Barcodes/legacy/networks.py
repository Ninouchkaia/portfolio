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



matrix1 = pd.read_csv("correl_matrix1_merged_pval_filtered_deseq2.csv", sep=';', header=0, index_col=0)


 
# Transform it in a links data frame (3 columns only):
links = matrix1.stack().reset_index()
links.columns = ['source', 'target', 'correlation']
 



# Keep only correlation over a threshold and remove self correlation (cor(A,A)=1)
# links_filtered=links.loc[ (links['correlation'] > 0) & (links['source'] != links['target']) ]
# links_filtered=links.loc[ (links['source'] != links['target']) ]
links_filtered=links.loc[ (links['correlation'] > 0) & (links['source'] != links['target']) ]
# Build your graph
G=nx.from_pandas_edgelist(links_filtered, 'source', 'target', 'correlation')
# G=nx.from_pandas_edgelist(links, 'source', 'target', 'correlation')
my_nodes = list(G.nodes)
print(len(my_nodes))
# print(my_nodes)



mapping = {}
for node in my_nodes :
	print(node)
	mapping[node] = '_'.join(node.split("_")[0:2])
	print(mapping[node])

print(len(mapping))

G = nx.relabel_nodes(G, mapping)

# pos = nx.circular_layout(G, scale = 2, dim = 2)#, seed=7)  # positions for all nodes - seed for reproducibility
pos = nx.spring_layout(G, k=0.1)#, scale=2)

# pos = nx.random_layout(G)

np.random.shuffle(my_nodes)

# nodes

color_dict = {}
with open("color_mapping.tsv", 'r') as file_read :
	data = file_read.readlines()
	for line in data :
		line = line.replace("\n", "").split("\t")
		color = int(line[0])
		if line[1] == 'Fluor_006u' :
			drug = '5Fluor_006u'
		elif line[1] == 'Azacyt_1,5u' :
			drug = 'Azacyt_1.5u'
		elif line[1] == 'Bafilo_1,2n' :
			drug = 'Bafilo_1.2n'
		else :
			drug = line[1]
		color_dict[drug] = color
print(color_dict)
colors = []
for node in G.nodes :
	colors.append(((color_dict[node]+1)/100))
print(colors)

# nx.draw_networkx_nodes(G, pos, node_size=100, node_color='yellow', alpha = 0.5)
nx.draw_networkx_nodes(G, pos, node_size=300, node_color=colors)#, alpha = 0.5)

# edges
elarge = [(u, v) for (u, v, d) in G.edges(data=True) if d['correlation'] > 0.8]
# esmall = [(u, v) for (u, v, d) in G.edges(data=True) if 0.5 < d['correlation'] <= 0.9]
nx.draw_networkx_edges(G, pos, edgelist=elarge, width=0.1, edge_color='black')#, alpha=0.5)
# nx.draw_networkx_edges(G, pos, edgelist=esmall, width=0.5, alpha=0.5, edge_color="b", style="dashed")

# labels

label_dict = {}
for node in G.nodes :
	label = node.split("_")[0]
	label_dict[node] = label
print(label_dict)
nx.draw_networkx_labels(G, pos, labels = label_dict, font_size=4, font_family="sans-serif")





 
# # Plot the network:
# # nx.draw(G, pos=nx.circular_layout(G), with_labels=True, node_color='orange', node_size=25, edge_color='pink', linewidths=1, font_size=6)
# # nx.draw(G, pos=nx.fruchterman_reingold_layout(G), with_labels=True, node_color='orange', node_size=250, edge_color='pink', linewidths=1, font_size=6)
# nx.draw(G, pos=nx.spring_layout(G), with_labels=True, node_color='orange', node_size=250, edge_color='pink', linewidths=1, font_size=6)

# plt.show()

plt.savefig("network_correlation_colored.pdf", format='pdf')
plt.savefig("network_correlation_colored.png", format='png')








# df2 = pd.read_csv("result6_fillna_control_renamed_filtered4_exp130921_T0_renamed_replicates_averages_fc_zeros.csv", sep=';', header=0, index_col=0)
# print (df2) #
# df2 = df2.astype(float)
# df2 = df2.sort_index(axis=1, ascending=True, key=lambda x: x.str.lower())

# matrix2 = df2.corr(
#     method = 'pearson',  # The method of correlation
#     min_periods = 1).round(2)     # Min number of observations required)



stop = timeit.default_timer()
print(stop - start) 