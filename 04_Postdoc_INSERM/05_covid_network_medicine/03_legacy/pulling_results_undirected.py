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









from scipy.stats import spearmanr
def corrfunc(x,y, ax=None, **kws):
    """Plot the correlation coefficient in the top left hand corner of a plot."""
    r, _ = spearmanr(x, y)
    ax = ax or plt.gca()
    # Unicode for lowercase rho (ρ)
    rho = '\u03C1'
    ax.annotate(f'{rho}_spearman = {r:.2f}', xy=(.1, .9), xycoords=ax.transAxes)

# data = pd.read_csv('pulled_results_COVID19_GDDS_proteins_relative_degrees_to_targets_undirected.tsv', usecols=['Drug', 'Human PPI (target)'], sep='\t')
data = pd.read_csv('pulled_results_COVID19_GDDS_proteins_relative_degrees_to_targets_undirected.tsv', sep='\t')

# Correlation Matrix Heatmap
# f, ax = plt.subplots(figsize=(10, 6))
# corr = data.corr(method='spearman')
# hm = sns.heatmap(round(corr,2), annot=True, ax=ax, cmap="coolwarm",fmt='.2f',
#                  linewidths=.05)
# f.subplots_adjust(top=0.93)
# t= f.suptitle('Protein Abs Degrees Undirected Correlation Heatmap', fontsize=14)
# plt.show()


g = sns.pairplot(data, corner=True, diag_kind="kde", kind="reg", dropna=True)
# g.map_lower(corrfunc)
g.fig.suptitle('Protein Relative Degrees Undirected Pairplot', x=0.2, y=1, ha='left')
plt.show()




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