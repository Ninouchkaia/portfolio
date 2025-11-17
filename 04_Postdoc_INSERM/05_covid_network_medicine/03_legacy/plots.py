#!/usr/bin/python
# -*- coding: utf-8 -*-
# jupyter notebook --browser='C:/Program Files (x86)/Google/Chrome/Application/chrome.exe %s'


import timeit
start = timeit.default_timer()

import sys
import os
import csv
from operator import itemgetter
import networkx as nx
from networkx.algorithms import community
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
from collections import Counter
from collections import defaultdict
from collections import OrderedDict
import statistics 
import re
import gzip
import math




# source_list = ['GO', 'Symptom', 'Disease','Drug', 'Human PPI (target)']
source_list = ['Drug', 'Symptom', 'Human PPI (target)','Disease','GO']



# zscore_dict = {}
# for source in source_list :
# 	node_type = source
# 	with open('COVID19_GDDS_Results_%s_to_CovidProteins_Relative-withZeros.tsv' % source, 'r') as file_read :
# 		data = file_read.readlines()
# 		for line in data[1:] :
# 			line = line.split('\t')
# 			node_name = line[0]
# 			zscore = float(line[-1].replace('\n', ''))
# 			if node_type not in zscore_dict :
# 				zscore_dict[node_type] = []
# 			zscore_dict[node_type].append(zscore)

# df = pd.DataFrame.from_dict(zscore_dict,  orient='index')
# print(len(df))
# # print(df)



deg_dict = {}
for source in source_list :
	node_type = source
	with open('COVID19_GDDS_%s_degrees_to_proteins.tsv' % source, 'r') as file_read :
		data = file_read.readlines()
		for line in data[1:] :
			line = line.split('\t')
			node_name = line[0]
			deg = float(line[1])
			if node_type not in deg_dict :
				deg_dict[node_type] = []
			deg_dict[node_type].append(deg)
df = pd.DataFrame.from_dict(deg_dict,  orient='index')



color_list = ["green", "blue", "orangered", "hotpink", "deepskyblue"]

fig, axes = plt.subplots(nrows=3, ncols=2, sharex=True, sharey=True, figsize=(15, 7))


source_n = 0 

for row in range(0,3) :
	for column in range(0,2) :
		if (row==0 and column==0) :
			for source in range(0,len(source_list)) :
				sns.kdeplot(df.transpose()[source_list[source]], color=color_list[source], shade=False, ax=axes[row][column], legend=False)
		else :
			sns.kdeplot(df.transpose()[source_list[source_n]], color=color_list[source_n], shade=True, legend=True, ax=axes[row][column])
			print(source_list[source_n], row, column)
			source_n = source_n + 1
# plt.savefig('all_zscores_density_final.png')
# plt.close()
fig.suptitle('Absolute degrees to proteins in covid networks - distributions',  fontsize=16)
plt.show()






### plot each plot separately
# for source in source_list :
# 	df[source].plot.kde(title='%s Zscores densities (%s-prot relative degrees)' % (source, source) )
# 	plt.savefig('%s_zscores_density.png' % source)
# 	# plt.show()
# 	plt.close()

# plot all densities in one plot
# df = df.transpose() # this is a trick to not be bothered about the columns not all having the same length.
# df.plot.kde(title='Zscores densities',colormap='Paired')
# plt.savefig('all_zscores_density_paired.png')
# plt.close()














# for j in range(0,len(color_list)) :
# 	sns.kdeplot(df.transpose()[source_list[j]], color=color_list[j], shade=True)
# plt.title('Zscores densities')
# plt.show()


# ncol = 2 # pick one dimension
# nrow = (len(df)+ ncol-1) / ncol # make sure enough subplots
# fig, ax = plt.subplots(nrows=nrow, ncols=ncol) # create the axes

# i = 0
# for source in source_list :
#   ix = np.unravel_index(i, ax.shape) # compute an appropriate index (1d or 2d)
#   df_transposed = df.transpose()
#   df_transposed[source].plot.kde(title='%s Zscores densities (%s-prot relative degrees)' % (source, source), ax=ax[ix] )
#   # ax[ix].plot(...)   # or direct axis object method plot (/scatter/bar/...)
#   plt.savefig('stacked_zscores_density.png')
#   i = i+1

# print(df.index.tolist())
# print(df)
# prop_cycle = plt.rcParams['axes.prop_cycle']
# colors = prop_cycle.by_key()['color']



# df1 = df.loc['GO']
# df1_trans = df1.transpose()
# df1_trans.plot.kde(ax=axes[0], legend=True, color='Paired', shade=True)
# df2 = df.loc['Symptom']
# df2_trans = df2.transpose()
# df2_trans.plot.kde(ax=axes[1], legend=True, color='Paired', shade=True)
# df3 = df.loc['Disease']
# df3_trans = df3.transpose()
# df3_trans.plot.kde(ax=axes[2], legend=True, color='Paired')
# df4 = df.loc['Drug']
# df4_trans = df4.transpose()
# df4_trans.plot.kde(ax=axes[3], legend=True, color='Paired')
# df5 = df.loc['Human PPI (target)']
# df5_trans = df5.transpose()
# df5_trans.plot.kde(ax=axes[4], legend=True, color='Paired')
# plt.show()



# ax = plt.gca()
# df = pd.DataFrame.from_dict(zscore_dict,  orient='index')
# df = df.transpose()
# print(df)
# for source in source_list:
# 	df_plot = df[source]
# 	# print(df_plot)
# 	df_plot.plot.kde(title='%s Zscores densities (%s-prot relative degrees)' % (source, source), ax=ax)
# 	# df_plot.plot.line(ax=ax)
# 	plt.savefig('zscores_density_test.png')




# df.plot.kde()
# df.groupby(['GO','Drug']).mean().unstack().plot()
# plt.show()

# COVID19_GDDS_proteins_degrees_to_Human PPI (target)

# target_list = ['GO', 'Symptom', 'Disease','Drug', 'Human PPI (target)']

# # degrees in covid-net
# deg_dict = {}
# for target in target_list :

# 	with open('full_analysis_for_proteins_undirected\\COVID19_GDDS_proteins_degrees_to_%s.tsv' % target, 'r') as file_read :
# 		data = file_read.readlines()
# 		for line in data :
# 			line = line.split('\t')
# 			protein = line[0]
# 			degree = float(line[1])
# 			degree_rel = float(line[2].replace('\n',''))
# 			if protein not in deg_dict :
# 				deg_dict[protein] = {}
# 			deg_dict[protein][target] = degree
 
# for protein in deg_dict :
# 	print(protein, deg_dict[protein])		

# with open('pulled_results_COVID19_GDDS_proteins_degrees_to_targets_undirected.tsv', 'a+') as file_write :
# 	#heading
# 	file_write.write('protein\t')
# 	for target in target_list[:-1] :
# 		file_write.write('%s\t' % target)
# 	file_write.write('%s\n' % target_list[-1])

# 	for protein in deg_dict :
# 		file_write.write('%s\t' % protein)
# 		for target in target_list[:-1] :
# 			file_write.write('%s\t' % deg_dict[protein][target])
# 		file_write.write('%s\n' % deg_dict[protein][target_list[-1]])









# target_list = ['GO', 'Symptom', 'Disease','Drug', 'Human PPI (target)']

# # degrees in covid-net
# deg_dict = {}
# for target in target_list :

# 	with open('full_analysis_for_proteins_undirected\\COVID19_GDDS_proteins_degrees_to_%s.tsv' % target, 'r') as file_read :
# 		data = file_read.readlines()
# 		for line in data :
# 			line = line.split('\t')
# 			protein = line[0]
# 			degree = float(line[1])
# 			degree_rel = float(line[2].replace('\n',''))
# 			if protein not in deg_dict :
# 				deg_dict[protein] = {}
# 			deg_dict[protein][target] = degree_rel
 
# for protein in deg_dict :
# 	print(protein, deg_dict[protein])		

# with open('pulled_results_COVID19_GDDS_proteins_relative_degrees_to_targets_undirected.tsv', 'a+') as file_write :
# 	#heading
# 	file_write.write('protein\t')
# 	for target in target_list[:-1] :
# 		file_write.write('%s\t' % target)
# 	file_write.write('%s\n' % target_list[-1])

# 	for protein in deg_dict :
# 		file_write.write('%s\t' % protein)
# 		for target in target_list[:-1] :
# 			file_write.write('%s\t' % deg_dict[protein][target])
# 		file_write.write('%s\n' % deg_dict[protein][target_list[-1]])









# target_list = ['GO', 'Symptom', 'Disease','Drug', 'Human PPI (target)']

# # degrees in covid-net
# score_dict = {}
# for target in target_list :

# 	with open('full_analysis_for_proteins_undirected\\COVID19_GDDS_Results_ForCovidProteins_to_%s_Relative-withZeros.tsv' % target, 'r') as file_read :
# 		data = file_read.readlines()
# 		for line in data[1:] :
# 			line = line.split('\t')
# 			protein = line[0]
# 			zscore = float(line[-1].replace('\n',''))
# 			if protein not in score_dict :
# 				score_dict[protein] = {}
# 			score_dict[protein][target] = zscore
 
# for protein in score_dict :
# 	print(protein, score_dict[protein])		

# with open('pulled_results_COVID19_GDDS_proteins_relative_zscores_to_targets_undirected.tsv', 'a+') as file_write :
# 	#heading
# 	file_write.write('protein\t')
# 	for target in target_list[:-1] :
# 		file_write.write('%s\t' % target)
# 	file_write.write('%s\n' % target_list[-1])

# 	for protein in score_dict :
# 		file_write.write('%s\t' % protein)
# 		for target in target_list[:-1] :
# 			file_write.write('%s\t' % score_dict[protein][target])
# 		file_write.write('%s\n' % score_dict[protein][target_list[-1]])









# from scipy.stats import spearmanr
# def corrfunc(x,y, ax=None, **kws):
#     """Plot the correlation coefficient in the top left hand corner of a plot."""
#     r, _ = spearmanr(x, y)
#     ax = ax or plt.gca()
#     # Unicode for lowercase rho (ρ)
#     rho = '\u03C1'
#     ax.annotate(f'{rho}_spearman = {r:.2f}', xy=(.1, .9), xycoords=ax.transAxes)

# # data = pd.read_csv('pulled_results_COVID19_GDDS_proteins_relative_degrees_to_targets_undirected.tsv', usecols=['Drug', 'Human PPI (target)'], sep='\t')
# data = pd.read_csv('pulled_results_COVID19_GDDS_proteins_relative_degrees_to_targets_undirected.tsv', sep='\t')

# # Correlation Matrix Heatmap
# # f, ax = plt.subplots(figsize=(10, 6))
# # corr = data.corr(method='spearman')
# # hm = sns.heatmap(round(corr,2), annot=True, ax=ax, cmap="coolwarm",fmt='.2f',
# #                  linewidths=.05)
# # f.subplots_adjust(top=0.93)
# # t= f.suptitle('Protein Abs Degrees Undirected Correlation Heatmap', fontsize=14)
# # plt.show()


# g = sns.pairplot(data, corner=True, diag_kind="kde", kind="reg", dropna=True)
# # g.map_lower(corrfunc)
# g.fig.suptitle('Protein Relative Degrees Undirected Pairplot', x=0.2, y=1, ha='left')
# plt.show()




# # Pair-wise Scatter Plots
# cols = ['GO', 'Symptom', 'Disease','Drug', 'Human PPI (target)']
# pp = sns.pairplot(data[cols], height=1.8, aspect=1.8,
#                   plot_kws=dict(edgecolor="k", linewidth=0.5),
#                   diag_kind="kde", diag_kws=dict(shade=True))

# fig = pp.fig 
# fig.subplots_adjust(top=0.93, wspace=0.3)
# t = fig.suptitle('Protein Relative Degree Zscores Pairwise Plots', fontsize=14)
# fig.show()


# sns.set(style="ticks", color_codes=True)
# iris = sns.load_dataset("iris")
# print(iris)
# g = sns.pairplot(iris)
# g = sns.pairplot(iris, hue="species")
# g = sns.pairplot(data, corner=True, diag_kind="kde", kind="reg")
# g.fig.suptitle('Protein Abs Degrees Pairplot', x=0.2, y=1, ha='left')
# plt.show()

# sns.set(style="ticks", color_codes=True)
# # g = sns.pairplot(iris)
# # g = sns.pairplot(iris, hue="species")
# # g = sns.pairplot(iris, corner=True, diag_kind="kde", kind="reg")
# g = sns.PairGrid(data)
# g = g.map(plt.scatter)
# plt.show()














# targets = [1,2,3,4,5]
# b = {}
# data = ['prot1','prot2','prot3']
# for target in targets :
# 	for prot in data :
# 		deg = 'deg_%s_%s' % (prot,target)
# 		degrel = 'degrel_%s_%s' % (prot,target)
# 		b[prot] = {target : (deg, degrel)}

# for i in b :
# 	print(i, b[i])







# with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\COVID19_GDDS_nodes.csv', 'r') as nodecsv: # Open the file                       
# 	nodereader = csv.reader(nodecsv) # Read the csv  
# 	nodes = [n for n in nodereader][1:]                     
# 	node_names = [n[0] for n in nodes] # Get a list of only the node names 

# with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\COVID19_GDDS_edges.csv', 'r') as edgecsv: # Open the file
# 	edgereader = csv.reader(edgecsv) # Read the csv     
# 	edges = [tuple(e) for e in edgereader][1:] # Retrieve the data

# G = nx.Graph()
# G.add_edges_from(edges)


# descr_dict = {}
# for node in nodes: 
# 	descr_dict[node[0]] = node[2]

# prots_linked_to_targets = []
# for u,v,c in G.edges(data=True) :	
# 	if (descr_dict[u] == 'Human PPI (target)' and descr_dict[v] == 'Drug') :
# 		prots_linked_to_targets.append(u)
# if 'ACPP' in prots_linked_to_targets :
# 	print('found ACPP')






stop = timeit.default_timer()
print(stop - start)  

