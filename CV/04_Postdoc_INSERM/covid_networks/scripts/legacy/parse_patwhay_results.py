#!/usr/bin/python
# -*- coding: utf-8 -*-
# jupyter notebook --browser='C:/Program Files (x86)/Google/Chrome/Application/chrome.exe %s'


import timeit
start = timeit.default_timer()

import sys
import os
import csv
from fractions import Fraction
from collections import OrderedDict

pathway_dict, pathway_dict_simple = {}, {}
with open("enrichPathway-2ndorder-pvalueCutoff0.005.tsv", 'r') as file_read :
	data = file_read.readlines()
	for line in data[1:] :
		line = line.split("\t")
		virus = line[0]
		pathway_id = line[1]
		pathway_description = line[2]
		# GeneRatio = float(Fraction((line[3])))	
		GeneRatio = line[3]	
		BgRatio	= float(Fraction((line[4])))	
		if line[5] != "NA" :
			pvalue = float(line[5])
		else : 
			pvalue = line[5]
		if line[6] != "NA" :
			padjust = float(line[6])
		else :
			padjust = line[6]	
		if line[7] != "NA" :	
			qvalue = float(line[7])
		else :
			qvalue = line[7]
		geneID = line[8]
		count = int(line[9].replace("\n",""))
		if virus not in pathway_dict :
			pathway_dict[virus] = [(pathway_description, GeneRatio, padjust)]
			pathway_dict_simple[virus] = [pathway_description]
		else :
			pathway_dict[virus].append((pathway_description, GeneRatio, padjust))
			pathway_dict_simple[virus].append(pathway_description)

# print(pathway_dict)

# get overlapping pathways between lacrosse and covid
overlap = list(set(pathway_dict_simple["Human_immunodeficiency_virus_type_1_group_M_subtype_B"]) & set(pathway_dict_simple["Human_SARS_coronavirus_2"]))
for i in overlap :
	print(i)

shared_pathways = {}
for virus in pathway_dict_simple :
	n_overlap = len(list(set(pathway_dict_simple[virus]) & set(pathway_dict_simple["Human_SARS_coronavirus_2"])))
	shared_pathways[virus] = n_overlap

for i in shared_pathways :
	print(i, shared_pathways[i])

# print(pathway_dict["Human_SARS_coronavirus_2"])
# my_list = []
# for i in pathway_dict :
# 	my_list.append(len(pathway_dict[i]))
# my_list=set(my_list)
# print(pathway_dict["Human_SARS_coronavirus_2"])

# res = sorted(pathway_dict["Human_SARS_coronavirus_2"], key=lambda x: x[2] )

# with open("sarscov2_second_order_enriched_pathway0.005.output.tsv" , 'a+') as file_write :
# 	for i in res :
# 		file_write.write("%s\t%s\t%s\n" % (i[0], i[2], i[1]))


stop = timeit.default_timer()
print(stop - start)  


