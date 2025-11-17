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


# covid_disease_list = []
# with open('COVID19_GDDS_Results_Diseases_to_Proteins_rel.tsv', 'r') as file_read :
# 	data = file_read.readlines()
# 	for line in data[1:] :
# 		line = line.split('\t')
# 		covid_disease = line[0]
# 		covid_disease_list.append(covid_disease)
# print(len(set(covid_disease_list)))	#4176	
# covid_disease_set = set(covid_disease_list)


# top_disease_list = []
# with open('ranked_diseases.tsv', 'r') as file_read :
# 	data = file_read.readlines()
# 	for line in data[1:359] :
# 		line = line.split('\t')
# 		covid_disease = line[0]
# 		top_disease_list.append(covid_disease)
# print(len(set(top_disease_list)))	#358	
# top_disease_set = set(top_disease_list)


# zhou_disease_list = []
# with open('ZhouSymptoms_41467_2014_BFncomms5212_MOESM1045_ESM.txt', 'r') as file_read :
# 	data = file_read.readlines()
# 	for line in data[1:] :
# 		line = line.split('\t')
# 		zhou_disease = line[1]
# 		zhou_disease_list.append(zhou_disease)
# print(len(set(zhou_disease_list)))	#4219
# zhou_disease_set = set(zhou_disease_list)

# overlap_diseases = covid_disease_set.intersection(zhou_disease_set)
# print(len(overlap_diseases)) #940

# overlap_top_diseases = zhou_disease_set.intersection(top_disease_set)
# print(len(overlap_top_diseases)) #84
# print(overlap_top_diseases) #84






# covid_symptom_list = []
# with open('COVID19_GDDS_Results_Symptoms_to_Proteins_rel.tsv', 'r') as file_read :
# 	data = file_read.readlines()
# 	for line in data[1:] :
# 		line = line.split('\t')
# 		covid_symptom = line[0]
# 		covid_symptom_list.append(covid_symptom)
# print(len(set(covid_symptom_list)))	#2157	
# covid_symptom_set = set(covid_symptom_list)

# zhou_symptom_list = []
# with open('ZhouSymptoms_41467_2014_BFncomms5212_MOESM1045_ESM.txt', 'r') as file_read :
# 	data = file_read.readlines()
# 	for line in data[1:] :
# 		line = line.split('\t')
# 		zhou_symptom = line[0]
# 		zhou_symptom_list.append(zhou_symptom)
# print(len(set(zhou_symptom_list)))	#322
# zhou_symptom_set = set(zhou_symptom_list)

# overlap_symptoms = covid_symptom_set.intersection(zhou_symptom_set)
# print(len(overlap_symptoms)) #60

# overlap_symptoms = zhou_symptom_set.intersection(covid_symptom_set)
# print(len(overlap_symptoms)) #60





covid_disease_list = []
covid_disease_scores = {}
with open('COVID19_GDDS_Results_Diseases_to_Proteins_rel.tsv', 'r') as file_read :
	data = file_read.readlines()
	for line in data[1:] :
		line = line.split('\t')
		covid_disease = line[0]
		score = float(line[-1].replace('\n',''))
		covid_disease_list.append(covid_disease)
		covid_disease_scores[covid_disease] = score
print(len(set(covid_disease_list)))	#4176	
covid_disease_set = set(covid_disease_list)


# covid_disease_scores_sorted = {k: v for k, v in sorted(covid_disease_scores.items(), key=lambda item: item[1])}
# print(covid_disease_scores_sorted)

covid_disease_scores_filtered = {k: v for k, v in sorted(covid_disease_scores.items(), key=lambda item: item[1]) if v > 2}
# print(covid_disease_scores_filtered)
print(len(covid_disease_scores_filtered)) #358 > 1 and 86 > 2, lets go with this one.
top_covid_diseases = [k for k in covid_disease_scores_filtered]
print(top_covid_diseases)
top_covid_diseases_set = set(top_covid_diseases)

with open("top_disease.tsv", 'a+') as file_write :
	for disease in covid_disease_scores_filtered :
		file_write.write("%s\t%s\n" % (disease,covid_disease_scores_filtered[disease])) 










# zhou_disease_list = []
# zhou_disease_symptom_dict = {}
# with open('ZhouSymptoms_41467_2014_BFncomms5212_MOESM1045_ESM.txt', 'r') as file_read :
# 	data = file_read.readlines()
# 	for line in data[1:] :
# 		line = line.split('\t')
# 		zhou_symptom = line[0]
# 		zhou_disease = line[1]
# 		zhou_disease_list.append(zhou_disease)
# 		if zhou_disease not in zhou_disease_symptom_dict :
# 			zhou_disease_symptom_dict[zhou_disease] = []
# 		zhou_disease_symptom_dict[zhou_disease].append(zhou_symptom)

# print(len(set(zhou_disease_list)))	#4219
# zhou_disease_set = set(zhou_disease_list)

# overlap_diseases = top_covid_diseases_set.intersection(zhou_disease_set)
# print(len(overlap_diseases)) #21
# print(overlap_diseases)
# #{'Brain Edema', 'Muscle Spasticity', 'Cerebral Hemorrhage', 'Nijmegen Breakage Syndrome', 'Pseudoxanthoma Elasticum', 'Contracture', 'Blindness, Cortical', 'Airway Obstruction', 'Learning Disorders', 'Hyperbilirubinemia, Neonatal', 'Mitochondrial Encephalomyopathies', 'Hepatitis C', 'Lymphohistiocytosis, Hemophagocytic', 'Neoplasm, Residual', 'Ehlers-Danlos Syndrome', 'Infection', 'Tricuspid Valve Prolapse', 'Severe Acute Respiratory Syndrome', 'Liver Failure', 'Retroviridae Infections', 'Adrenal Gland Neoplasms'}
# disease_overlap_dict = {}
# for disease in overlap_diseases :
# 	disease_overlap_dict[disease] = covid_disease_scores[disease]
# disease_overlap_dict_sorted = {k: v for k, v in sorted(disease_overlap_dict.items(), key=lambda item: item[1])}

# symptoms_occurrences = []
# for disease in disease_overlap_dict_sorted :
# 	#print(disease, disease_overlap_dict_sorted[disease], zhou_disease_symptom_dict[disease])
# 	symptoms_occurrences.append(zhou_disease_symptom_dict[disease])

# # print(symptoms_occurrences)

# symptoms_counts = []
# for symptoms_list in symptoms_occurrences :
# 	for symptom in symptoms_list :
# 		symptoms_counts.append(symptom)
# # print(symptoms_counts)

# symptoms_counts = sorted(symptoms_counts)
# # print(symptoms_counts)

# symptoms_counts_dict = defaultdict( int )
# for symptom in  symptoms_counts:
#     symptoms_counts_dict[symptom] += 1
# symptoms_counts_dict_sorted = {k: v for k, v in sorted(symptoms_counts_dict.items(), key=lambda item: item[1])}
# for symptom in symptoms_counts_dict_sorted :
# 	print(symptom, symptoms_counts_dict_sorted[symptom])


stop = timeit.default_timer()
print(stop - start)  

