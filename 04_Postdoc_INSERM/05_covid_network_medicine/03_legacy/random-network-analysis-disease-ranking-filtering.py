#!/usr/bin/python
# -*- coding: utf-8 -*-
# jupyter notebook --browser='C:/Program Files (x86)/Google/Chrome/Application/chrome.exe %s'


import timeit
start = timeit.default_timer()

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





#####################################################################################################
#########################  Compute network Diseases degree distribution #################################
##########################################  based on all edges  #####################################
#####################################################################################################

with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\COVID19_GDDS_nodes.csv', 'r') as nodecsv: # Open the file                       
	nodereader = csv.reader(nodecsv) # Read the csv  
	nodes = [n for n in nodereader][1:]                     
	node_names = [n[0] for n in nodes] # Get a list of only the node names 

with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\COVID19_GDDS_edges.csv', 'r') as edgecsv: # Open the file
	edgereader = csv.reader(edgecsv) # Read the csv     
	edges = [tuple(e) for e in edgereader][1:] # Retrieve the data

G = nx.Graph()
G.add_nodes_from(node_names)
G.add_edges_from(edges)

descr_dict = {}
for node in nodes: 
	descr_dict[node[0]] = node[2]

print(nx.info(G))
# Number of nodes: 16007
# Number of edges: 53917
# Average degree:   6.7367

# print(nx.degree(G)) # this is the list containing as many tuples as nodes, indicating each node's degree.

average_degree = statistics.mean([ tpl[1] for tpl in nx.degree(G) ])
# print(average_degree) # 6.736677703504717

diseases_degrees = nx.degree(G, nbunch=[n for n in G.nodes if descr_dict[n] == 'Disease'])
# print("degree of disease nodes only", diseases_degrees)

average_degree_of_diseases = statistics.mean([ tpl[1] for tpl in diseases_degrees ])
print("average degree of diseases nodes", average_degree_of_diseases) # 4.329980842911877

# #####################################################################################################
# ################################# CALCULATION ON EDGES WITH FILTERING  ###########################
# #####################################################################################################
# # let's count edges only in order to calculate the nodes' degrees based on the number of edges they make in the graph that connect them to proteins that are themselves connected to viral genes

#first we have to build the set of proteins of interest
proteins_of_interest = set()
for u,v,c in G.edges(data=True) :	
	if (descr_dict[u] == 'Viral Gene' and descr_dict[v] == 'Human PPI (target)') :
		proteins_of_interest.add(v)
print(len(proteins_of_interest)) # 332 OK

# now we keep only the disease that interact with those proteins of interest
diseases_linked_to_some_proteins_of_interest = []
for u,v,c in G.edges(data=True) :	
	if (descr_dict[u] == 'Disease' and v in proteins_of_interest) :
		diseases_linked_to_some_proteins_of_interest.append(u) # this never happens
	if (descr_dict[v] == 'Disease' and u in proteins_of_interest) :
		diseases_linked_to_some_proteins_of_interest.append(v)
print(len(diseases_linked_to_some_proteins_of_interest), len(set(diseases_linked_to_some_proteins_of_interest))) # 9430 3164
# ({'Drug': 5703, 'Disease': 4176, 'GO': 3487, 'Symptom': 2157, 'Human PPI (target)': 457, 'Viral Gene': 27})

# # print(sorted(diseases_linked_to_some_proteins_of_interest))

disease_counts = sorted(diseases_linked_to_some_proteins_of_interest)

# # # for the diseases
disease_counts_dict = defaultdict( int )
for disease in disease_counts:
    disease_counts_dict[disease] += 1

my_diseases_nodes = [n for n in G.nodes if descr_dict[n] == 'Disease']
for disease in my_diseases_nodes :
	if disease_counts_dict.get(disease) == None :
		# print("%s is missing in the linked diseases : Adding it with a frequency of 0.\n" % disease)
		disease_counts_dict[disease]= 0

average_degree_of_my_diseases_nodes = statistics.mean([ disease_counts_dict[disease] for disease in disease_counts_dict ])
print("average_degree_of_my_diseases_nodes", average_degree_of_my_diseases_nodes)
# 2.2581417624521074

max_degree_of_my_diseases_nodes = max([ disease_counts_dict[disease] for disease in disease_counts_dict ])
print("max_degree_of_my_diseases_nodes", max_degree_of_my_diseases_nodes) #93


ordered_disease_counts_dict = OrderedDict(sorted(disease_counts_dict.items()))

# with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\Further_Analysis\\Diseases\\First_order_proteins\\COVID19_GDDS_diseases_degrees_to_first_order_proteins.tsv', 'a+') as file_write:
# 	for disease in ordered_disease_counts_dict :
# 		file_write.write("%s\t%s\t%s\n" % (disease, ordered_disease_counts_dict[disease], ordered_disease_counts_dict[disease] / sum(ordered_disease_counts_dict.values())))


# # ##################################################################
# # ############ SAME ANALYSIS FOR RANDOM NETWORKS ###################
# # # ##################################################################

# with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\Further_Analysis\\Diseases\\First_order_proteins\\Log_degree_analysis_on_random_networks_for_diseases_proteins_of_interest.txt', 'a+') as log_write:

# 	nodefiles_numbers = []
# 	edgefiles_numbers = []
# 	for filename in os.listdir("A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\random_networks\\Mock_networks_21Apr"):
# 		if 'nodes_' in filename :
# 			filename = filename.split(".csv")[0]
# 			nodefiles_numbers.append(int(filename[19:]))
# 		if 'edges_' in filename :
# 			filename = filename.split(".csv")[0]
# 			edgefiles_numbers.append(int(filename[19:]))

# 	log_write.write("len(nodefiles_numbers) = %s\n" % len(nodefiles_numbers)) #2051
# 	log_write.write("len(edgefiles_numbers) = %s\n" % len(edgefiles_numbers)) #2051
# 	log_write.write("difference between edgefiles_numbers and nodefiles_numbers = %s\n" % len(set(nodefiles_numbers).difference(set(edgefiles_numbers)))) #0

# 	nodefiles_numbers = sorted(nodefiles_numbers)
# 	edgefiles_numbers = sorted(edgefiles_numbers)

# 	for iteration in nodefiles_numbers[1:] :
# 		print("\n########################## ITERATION %s ###########################\n" % iteration)
# 		log_write.write("\n########################## ITERATION %s ###########################\n" % iteration)
	
# 		with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\random_networks\\Mock_networks_21Apr\\COVID19_GDDS_nodes_%s.csv' % str(iteration), 'r') as nodecsv: # Open the file                   
# 			nodereader = csv.reader(nodecsv) # Read the csv  
# 			nodes = [n for n in nodereader][1:]                     
# 			node_names = [n[0] for n in nodes] # Get a list of only the node names 

# 		with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\random_networks\\Mock_networks_21Apr\\COVID19_GDDS_edges_%s.csv' % str(iteration), 'r') as edgecsv: # Open the file
# 			edgereader = csv.reader(edgecsv) # Read the csv     
# 			edges = [tuple(e) for e in edgereader][1:] # Retrieve the data

# 			log_write.write("len(node_names) : %s\n" % len(node_names)) #16007 --> 13448
# 			log_write.write("len(edges) : %s\n" % len(edges)) #53917 --> 37835

# 		G = nx.Graph()
# 		G.add_nodes_from(node_names)
# 		G.add_edges_from(edges)

# 		descr_dict = {}
# 		description_set=set()
# 		for node in nodes: 
# 			descr_dict[node[0]] = node[2]
# 			description_set.add(node[2])

# 		# print(description_set) #{'Symptom', 'Human PPI (target)', 'GO', 'Disease', 'Drug'}
# 		# recounted = Counter(list(descr_dict.values()))
# 		# print("%s\n" % recounted)

# 		# print(nx.info(G))
# 		log_write.write("\nRANDOM GRAPH %s\n" % iteration)
# 		log_write.write("%s\n" % nx.info(G))

# 		#####################################################################################################
# 		################################# GET EDGES CONNECTING DISEASES TO FIRST ORDER PROTEINS ################
# 		#####################################################################################################
# 		#first we have to build the set of proteins of interest
# 		proteins_of_interest = set([n for n in G.nodes if descr_dict[n] == 'Human PPI (target)'])
# 		# print(proteins_of_interest)

# 		# now we keep only the diseases that interact with those proteins of interest
# 		diseases_linked_to_some_proteins_of_interest = []
# 		for u,v,c in G.edges(data=True) :	
# 			if (descr_dict[u] == 'Disease' and v in proteins_of_interest) :
# 				diseases_linked_to_some_proteins_of_interest.append(u) # this never happens
# 			if (descr_dict[v] == 'Disease' and u in proteins_of_interest) :
# 				diseases_linked_to_some_proteins_of_interest.append(v)
# 		# print(len(diseases_linked_to_some_proteins_of_interest), len(set(diseases_linked_to_some_proteins_of_interest))) # 8200 3740
		
# 		disease_counts = sorted(diseases_linked_to_some_proteins_of_interest)

# 		# for the diseases
# 		disease_counts_dict = defaultdict( int )
# 		for disease in disease_counts:
# 		    disease_counts_dict[disease] += 1

# 		my_diseases_nodes = [n for n in G.nodes if descr_dict[n] == 'Disease']
# 		for disease in my_diseases_nodes :
# 			if disease_counts_dict.get(disease) == None :
# 				# print("%s is missing in the linked disease : Adding it with a frequency of 0.\n" % disease)
# 				disease_counts_dict[disease]= 0

# 		average_degree_of_my_diseases_nodes = statistics.mean([ disease_counts_dict[disease] for disease in disease_counts_dict ])
# 		# print("average_degree_of_my_diseases_nodes", average_degree_of_my_diseases_nodes)
# 		log_write.write("Average degree of diseases to a protein : %s\n" % average_degree_of_my_diseases_nodes)

# 		ordered_disease_counts_dict = OrderedDict(sorted(disease_counts_dict.items()))

# 		with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\Further_Analysis\\Random_Networks\\Diseases\\First_order_proteins\\COVID19_GDDS_diseases_degrees_to_first_order_proteins_%s.tsv' % str(iteration), 'a+') as file_write:
# 			for disease in ordered_disease_counts_dict :
# 				if (sum(ordered_disease_counts_dict.values()) != 0) :
# 					file_write.write("%s\t%s\t%s\n" % (disease, ordered_disease_counts_dict[disease], ordered_disease_counts_dict[disease] / sum(ordered_disease_counts_dict.values())))
# 				else :
# 					file_write.write("%s\t%s\t%s\n" % (disease, ordered_disease_counts_dict[disease], ordered_disease_counts_dict[disease]))



# # ##################################################################
# # ############ COMPUTE DEGREES MEANS AND Z-SCORES ##################
# # ######## INCLUDING ABSENCE OF DISEASES AS ZERO DEGREES ###########
# # ##################################################################


# # For covid network, compute the mean and the standard deviation of diseases degrees and normalized degrees
covid_diseases = []
covid_disease_degree_dict = {}
covid_disease_norm_degree_dict = {}

total_disease_degree_dict = defaultdict(list)
total_disease_normalized_degree_dict = defaultdict(list)

with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\Further_Analysis\\Diseases\\First_order_proteins\\COVID19_GDDS_diseases_degrees_to_first_order_proteins.tsv', 'r') as data :
	for line in data :
		line = line.split('\t')
		disease = line[0]
		covid_diseases.append(disease)
		disease_degree = int(line[1])
		disease_norm_degree = float(line[2].replace('\n', ''))
		if covid_disease_degree_dict.get(disease) == None :
			covid_disease_degree_dict[disease] = disease_degree
		else :
			print("warning : %s duplicate ?!" % disease)
		covid_disease_norm_degree_dict[disease] = disease_norm_degree

		if total_disease_degree_dict.get(disease) == None :
			total_disease_degree_dict[disease] = [disease_degree]
		else :
			total_disease_degree_dict[disease].append(disease_degree)

		if total_disease_normalized_degree_dict.get(disease) == None :
			total_disease_normalized_degree_dict[disease] = [disease_norm_degree]
		else :
			total_disease_normalized_degree_dict[disease].append(disease_norm_degree)

mean_disease_deg = statistics.mean(covid_disease_degree_dict.values())
mean_norm_disease_deg = statistics.mean(covid_disease_norm_degree_dict.values())
sd_disease_deg = statistics.pstdev(covid_disease_degree_dict.values())
sd_norm_disease_deg = statistics.pstdev(covid_disease_norm_degree_dict.values())
print(mean_disease_deg, sd_disease_deg, mean_norm_disease_deg, sd_norm_disease_deg) #1.4378397334736104 2.251958221466831 0.0001753463089601964 0.00027462905139839403

print(len(covid_diseases),len(covid_disease_degree_dict),len(covid_disease_norm_degree_dict)) #5703 5703 5703


nodefiles_numbers = []
edgefiles_numbers = []
for filename in os.listdir("A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\random_networks\\Mock_networks_21Apr"):
	if 'nodes_' in filename :
		filename = filename.split(".csv")[0]
		nodefiles_numbers.append(int(filename[19:]))
	if 'edges_' in filename :
		filename = filename.split(".csv")[0]
		edgefiles_numbers.append(int(filename[19:]))

# log_write.write("len(nodefiles_numbers) = %s\n" % len(nodefiles_numbers)) #2051
# log_write.write("len(edgefiles_numbers) = %s\n" % len(edgefiles_numbers)) #2051
# log_write.write("difference between edgefiles_numbers and nodefiles_numbers = %s\n" % len(set(nodefiles_numbers).difference(set(edgefiles_numbers)))) #0
print("len(nodefiles_numbers) = %s\n" % len(nodefiles_numbers)) #2051
print("len(edgefiles_numbers) = %s\n" % len(edgefiles_numbers)) #2051
print("difference between edgefiles_numbers and nodefiles_numbers = %s\n" % len(set(nodefiles_numbers).difference(set(edgefiles_numbers)))) #0

nodefiles_numbers = sorted(nodefiles_numbers)
edgefiles_numbers = sorted(edgefiles_numbers)


for iteration in nodefiles_numbers[1:] :
	# log_write.write("\n########################## ITERATION %s ###########################\n" % iteration)
	print("\n########################## ITERATION %s ###########################\n" % iteration)

	with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\Further_Analysis\\Random_Networks\\Diseases\\First_order_proteins\\COVID19_GDDS_diseases_degrees_to_first_order_proteins_%s.tsv' % str(iteration), 'r') as data: # Open the file  
		for line in data :
			line = line.split('\t')
			disease = line[0]
			degree = float(line[1])
			normalized_degree = float(line[2].replace('\n', ''))

			if total_disease_degree_dict.get(disease) == None :
				total_disease_degree_dict[disease] = [degree]
			else :
				total_disease_degree_dict[disease].append(degree)

			if total_disease_normalized_degree_dict.get(disease) == None :
				total_disease_normalized_degree_dict[disease] = [normalized_degree]
			else :
				total_disease_normalized_degree_dict[disease].append(normalized_degree)

print(len(total_disease_degree_dict)) # 72069
print(len(total_disease_normalized_degree_dict)) # 72069

### make the files with all calculations
with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\Further_Analysis\\Diseases\\First_order_proteins\\COVID19_GDDS_Results_Diseases_first_order_proteins-withZeros.tsv', 'a+') as file_write:
	file_write.write("disease\t# of nets it is found in\tmin-degree\tmax-degree\tmean-degree\tmedian-degree\tsd-degree\n")
	for disease in total_disease_degree_dict :
		nets_number = len(total_disease_degree_dict[disease])
		min_deg = min(total_disease_degree_dict[disease])
		max_deg = max(total_disease_degree_dict[disease])
		mean_deg = statistics.mean(total_disease_degree_dict[disease])
		median_deg = statistics.median(total_disease_degree_dict[disease])
		sd_deg = statistics.pstdev(total_disease_degree_dict[disease])
		file_write.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (disease, nets_number, min_deg, max_deg, mean_deg, median_deg, sd_deg))


with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\Further_Analysis\\Diseases\\\First_order_proteins\\COVID19_GDDS_Results_Diseases_first_order_proteins_Relative-withZeros.tsv', 'a+') as file_write:
	file_write.write("disease\t# of nets it is found in\tmin-degree-rel\tmax-degree-rel\tmean-degree-rel\tmedian-degree-rel\tsd-degree-rel\n")
	for disease in total_disease_normalized_degree_dict :
		nets_number = len(total_disease_normalized_degree_dict[disease])
		min_deg = min(total_disease_normalized_degree_dict[disease])
		max_deg = max(total_disease_normalized_degree_dict[disease])
		mean_deg = statistics.mean(total_disease_normalized_degree_dict[disease])
		median_deg = statistics.median(total_disease_normalized_degree_dict[disease])
		sd_deg = statistics.pstdev(total_disease_normalized_degree_dict[disease])
		file_write.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (disease, nets_number, min_deg, max_deg, mean_deg, median_deg, sd_deg))

### compute Z-score for the covid diseases
# For each disease found in covid network, compute 2 Z-scores.

with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\Further_Analysis\\Diseases\\First_order_proteins\\COVID19_GDDS_Results_ForCovidDiseases-first_order_proteins-withZeros.tsv', 'a+') as file_write:
	file_write.write("covid disease\t# of nets it is found in\tmin-degree\tmax-degree\tmean-degree\tmedian-degree\tsd-degree\tcovid-degree\tZ-score\n")
	for disease in covid_diseases :
		nets_number = len(total_disease_degree_dict[disease])
		min_deg = min(total_disease_degree_dict[disease])
		max_deg = max(total_disease_degree_dict[disease])
		mean_deg = statistics.mean(total_disease_degree_dict[disease])
		median_deg = statistics.median(total_disease_degree_dict[disease])
		sd_deg = statistics.pstdev(total_disease_degree_dict[disease])
		covid_deg = covid_disease_degree_dict[disease]
		if sd_deg != 0 :
			z_score = ( covid_deg - mean_deg) / sd_deg
		else : 
			z_score = 'NaN'
		file_write.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (disease, nets_number, min_deg, max_deg, mean_deg, median_deg, sd_deg, covid_deg, z_score))


with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\Further_Analysis\\Diseases\\First_order_proteins\\COVID19_GDDS_Results_ForCovidDiseases_first_order_proteins_Relative-withZeros.tsv', 'a+') as file_write:
	file_write.write("covid disease\t# of nets it is found in\tmin-degree-rel\tmax-degree-rel\tmean-degree-rel\tmedian-degree-rel\tsd-degree-rel\tcovid-degree-rel\tZ-score-rel\n")
	for disease in covid_diseases :
		nets_number = len(total_disease_normalized_degree_dict[disease])
		min_deg = min(total_disease_normalized_degree_dict[disease])
		max_deg = max(total_disease_normalized_degree_dict[disease])
		mean_deg = statistics.mean(total_disease_normalized_degree_dict[disease])
		median_deg = statistics.median(total_disease_normalized_degree_dict[disease])
		sd_deg = statistics.pstdev(total_disease_normalized_degree_dict[disease])
		covid_deg = covid_disease_norm_degree_dict[disease]
		if sd_deg != 0 :
			z_score = ( covid_deg - mean_deg) / sd_deg
		else : 
			z_score = 'NaN'
		file_write.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (disease, nets_number, min_deg, max_deg, mean_deg, median_deg, sd_deg, covid_deg, z_score))










stop = timeit.default_timer()
print(stop - start) 