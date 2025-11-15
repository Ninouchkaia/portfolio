#!/usr/bin/python
# -*- coding: utf-8 -*-
# jupyter notebook --browser='C:/Program Files (x86)/Google/Chrome/Application/chrome.exe %s'

import timeit
start = timeit.default_timer()

from scipy.cluster.hierarchy import linkage, dendrogram, fcluster
from scipy.spatial.distance import squareform
import glob
import statistics

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
 

with open("result6_fillna_control_renamed_filtered4.csv") as file_read:
	data = file_read.readlines()
	### parse headers
	drugs_list, exp_list = [], []
	drugs_index_dict, exp_index_dict = {}, {}
	barcode_dict = {}
	index_content_dict = {}
	drug_dict = {}
	exp_dict = {}
	for line in data[:1] :
		line = line.replace("\n","").split(";")
		for condition in line[1:] :
			conditions = condition.split("_")
			drug = conditions[0][:-1]
			replicate = conditions[0][-1]
			print(drug,replicate)
			drugs_list.append(drug)
			exp = conditions[2]
			exp_list.append(exp)

			if drug not in drugs_index_dict :
				drugs_index_dict[drug] = [line.index(condition)]
			else :
				drugs_index_dict[drug].append(line.index(condition))
			if exp not in exp_index_dict :
				exp_index_dict[exp] = [line.index(condition)]
			else :
				exp_index_dict[exp].append(line.index(condition))

			index_content_dict[line.index(condition)] = (drug,replicate,exp)

			################# fill the drug_dict ############################
			if drug not in drug_dict :
				drug_dict[drug] = {}
				
			if exp not in drug_dict[drug] :
				drug_dict[drug][exp] = {}
			
			if replicate not in drug_dict[drug][exp] :
				drug_dict[drug][exp][replicate] = [line.index(condition)]
			else :
				drug_dict[drug][exp][replicate].append(line.index(condition))
			#################################################################

			################# fill the exp_dict #############################
			if exp not in exp_dict :
				exp_dict[exp] = {}
				
			if drug not in exp_dict[exp] :
				exp_dict[exp][drug] = {}
			
			if replicate not in exp_dict[exp][drug] :
				exp_dict[exp][drug][replicate] = [line.index(condition)]
			else :
				exp_dict[exp][drug][replicate].append(line.index(condition))			

				


	drugs_list = list(set(drugs_list))
	exp_list = list(set(exp_list))

print("drugs_list", drugs_list)
print("exp_list", exp_list)	

for drug in drugs_index_dict :
	print(drug, drugs_index_dict[drug])

for exp in exp_index_dict :
	print(exp, exp_index_dict[exp])

# for index in index_content_dict :
# 	# print(index, index_content_dict[index])
# 	print(index, index_content_dict[index][0], index_content_dict[index][1], index_content_dict[index][2])

for drug in drug_dict :
	print(drug, drug_dict[drug])

for exp in exp_dict :
	print(exp, exp_dict[exp])

#### build file displaying the reads numbers as fold changes compared to 
with open("result6_fillna_control_renamed_filtered4.csv") as file_read:
	data = file_read.readlines()
	for line in data[1:] :
		line = line.replace("\n","").split(";")
		barcode = line[0]
		reads = line[1:]


stop = timeit.default_timer()
print(stop - start) 