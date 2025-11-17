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
#########################  Compute network Drug degree distribution #################################
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
print(average_degree) # 6.736677703504717

drugs_degrees = nx.degree(G, nbunch=[n for n in G.nodes if descr_dict[n] == 'Drug'])
# print("degree of drugs nodes only", drugs_degrees)

average_degree_of_drugs = statistics.mean([ tpl[1] for tpl in drugs_degrees ])
print("average degree of drugs nodes", average_degree_of_drugs) # 2.2826582500438364

# #####################################################################################################
# ################################# CALCULATION ON EDGES WITHOUT FILTERING  ###########################
# #####################################################################################################
# # let's count edges only in order to calculate the nodes' degrees based on the number of edges they make in the graph

drugs_linked_to_some_nodes = []
for u,v,c in G.edges(data=True) :	
	if (descr_dict[u] == 'Drug') :
		drugs_linked_to_some_nodes.append(u) # this never happens
	if (descr_dict[v] == 'Drug') :
		drugs_linked_to_some_nodes.append(v)
print(len(drugs_linked_to_some_nodes), len(set(drugs_linked_to_some_nodes))) # 13018 5703
# ({'Drug': 5703, 'Disease': 4176, 'GO': 3487, 'Symptom': 2157, 'Human PPI (target)': 457, 'Viral Gene': 27})

# print(sorted(drugs_linked_to_some_nodes))

drug_counts = sorted(drugs_linked_to_some_nodes)

# # for the drugs
drug_counts_dict = defaultdict( int )
for drug in drug_counts:
    drug_counts_dict[drug] += 1

my_drugs_nodes = [n for n in G.nodes if descr_dict[n] == 'Drug']
for drug in my_drugs_nodes :
	if drug_counts_dict.get(drug) == None :
		print("%s is missing in the linked drugs : Adding it with a frequency of 0.\n" % drug)
		drug_counts_dict[drug]= 0

average_degree_of_my_drugs_nodes = statistics.mean([ drug_counts_dict[drug] for drug in drug_counts_dict ])
# print("average_degree_of_my_drugs_nodes", average_degree_of_my_drugs_nodes)
# 2.2826582500438364

ordered_drug_counts_dict = OrderedDict(sorted(drug_counts_dict.items()))

# with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\Further_Analysis\\Drugs\\COVID19_GDDS_drugs_degrees_to_all_nodes.tsv', 'a+') as file_write:
# 	for drug in ordered_drug_counts_dict :
# 		file_write.write("%s\t%s\t%s\n" % (drug, ordered_drug_counts_dict[drug], ordered_drug_counts_dict[drug] / sum(ordered_drug_counts_dict.values())))


# ##################################################################
# ############ SAME ANALYSIS FOR RANDOM NETWORKS ###################
# ##################################################################

# with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\Further_Analysis\\Drugs\\Log_degree_analysis_on_random_networks_for_drugs_all_nodes.txt', 'a+') as log_write:

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
# 		for node in nodes: 
# 			descr_dict[node[0]] = node[2]

# 		print(nx.info(G))
# 		log_write.write("\nRANDOM GRAPH %s\n" % iteration)
# 		log_write.write("%s\n" % nx.info(G))

# 		#####################################################################################################
# 		################################# GET EDGES CONNECTING DRUGS TO ANY NODES  ################
# 		#####################################################################################################
# 		# We fetch all edges based on drug-nodes or nodes-drug links, count these edges, 
# 		# in order to calculate the nodes' degrees based on their number.

# 		drugs_linked_to_some_nodes = []
# 		for u,v,c in G.edges(data=True) :	
# 			if (descr_dict[u] == 'Drug') :
# 				drugs_linked_to_some_nodes.append(u)
# 			if (descr_dict[v] == 'Drug') :
# 				drugs_linked_to_some_nodes.append(v)
# 		print(len(drugs_linked_to_some_nodes), len(set(drugs_linked_to_some_nodes))) #
# 		log_write.write("Number of edges linking drugs to any nodes : %s \n" % len(drugs_linked_to_some_nodes))
# 		log_write.write("Number of drugs with at least one edge to a node: %s \n" % len(set(drugs_linked_to_some_nodes)))
		
# 		drug_counts = sorted(drugs_linked_to_some_nodes)

# 		# for the drugs
# 		drug_counts_dict = defaultdict( int )
# 		for drug in drug_counts:
# 		    drug_counts_dict[drug] += 1

# 		my_drugs_nodes = [n for n in G.nodes if descr_dict[n] == 'Drug']
# 		for drug in my_drugs_nodes :
# 			if drug_counts_dict.get(drug) == None :
# 				# print("%s is missing in the linked drugs : Adding it with a frequency of 0.\n" % drug)
# 				drug_counts_dict[drug]= 0

# 		average_degree_of_my_drugs_nodes = statistics.mean([ drug_counts_dict[drug] for drug in drug_counts_dict ])
# 		print("average_degree_of_my_drugs_nodes", average_degree_of_my_drugs_nodes)
# 		log_write.write("Average degree of drugs to any nodes : %s\n" % average_degree_of_my_drugs_nodes)

# 		ordered_drug_counts_dict = OrderedDict(sorted(drug_counts_dict.items()))

# 		with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\Further_Analysis\\Random_Networks\\Drugs\\COVID19_GDDS_drugs_degrees_to_all_nodes_%s.tsv' % str(iteration), 'a+') as file_write:
# 			for drug in ordered_drug_counts_dict :
# 				file_write.write("%s\t%s\t%s\n" % (drug, ordered_drug_counts_dict[drug], ordered_drug_counts_dict[drug] / sum(ordered_drug_counts_dict.values())))




# ##################################################################
# ############ COMPUTE DEGREES MEANS AND Z-SCORES ##################
# ######## INCLUDING ABSENCE OF DISEASES AS ZERO DEGREES ###########
# ##################################################################


# For covid network, compute the mean and the standard deviation of drugs degrees and normalized degrees
covid_drugs = []
covid_drug_degree_dict = {}
covid_drug_norm_degree_dict = {}

total_drug_degree_dict = defaultdict(list)
total_drug_normalized_degree_dict = defaultdict(list)

with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\Further_Analysis\\Drugs\\COVID19_GDDS_drugs_degrees_to_all_nodes.tsv', 'r') as data :
	for line in data :
		line = line.split('\t')
		drug = line[0]
		covid_drugs.append(drug)
		drug_degree = int(line[1])
		drug_norm_degree = float(line[2].replace('\n', ''))
		if covid_drug_degree_dict.get(drug) == None :
			covid_drug_degree_dict[drug] = drug_degree
		else :
			print("warning : %s duplicate ?!" % drug)
		covid_drug_norm_degree_dict[drug] = drug_norm_degree

		if total_drug_degree_dict.get(drug) == None :
			total_drug_degree_dict[drug] = [drug_degree]
		else :
			total_drug_degree_dict[drug].append(drug_degree)

		if total_drug_normalized_degree_dict.get(drug) == None :
			total_drug_normalized_degree_dict[drug] = [drug_norm_degree]
		else :
			total_drug_normalized_degree_dict[drug].append(drug_norm_degree)

mean_drug_deg = statistics.mean(covid_drug_degree_dict.values())
mean_norm_drug_deg = statistics.mean(covid_drug_norm_degree_dict.values())
sd_drug_deg = statistics.pstdev(covid_drug_degree_dict.values())
sd_norm_drug_deg = statistics.pstdev(covid_drug_norm_degree_dict.values())
print(mean_drug_deg, sd_drug_deg, mean_norm_drug_deg, sd_norm_drug_deg) #2.2826582500438364 2.9197086133682784 0.0001753463089601964 0.0002242824253624426

print(len(covid_drugs),len(covid_drug_degree_dict),len(covid_drug_norm_degree_dict)) #5703 5703 5703


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

	with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\Further_Analysis\\Random_Networks\\Drugs\\COVID19_GDDS_drugs_degrees_to_all_nodes_%s.tsv' % str(iteration), 'r') as data: # Open the file  
		for line in data :
			line = line.split('\t')
			drug = line[0]
			degree = float(line[1])
			normalized_degree = float(line[2].replace('\n', ''))

			if total_drug_degree_dict.get(drug) == None :
				total_drug_degree_dict[drug] = [degree]
			else :
				total_drug_degree_dict[drug].append(degree)

			if total_drug_normalized_degree_dict.get(drug) == None :
				total_drug_normalized_degree_dict[drug] = [normalized_degree]
			else :
				total_drug_normalized_degree_dict[drug].append(normalized_degree)

print(len(total_drug_degree_dict)) # 72069
print(len(total_drug_normalized_degree_dict)) # 72069

### make the files with all calculations
with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\Further_Analysis\\Drugs\\COVID19_GDDS_Results_Drugs_All_Nodes-withZeros.tsv', 'a+') as file_write:
	file_write.write("drug\t# of nets it is found in\tmin-degree\tmax-degree\tmean-degree\tmedian-degree\tsd-degree\n")
	for drug in total_drug_degree_dict :
		nets_number = len(total_drug_degree_dict[drug])
		min_deg = min(total_drug_degree_dict[drug])
		max_deg = max(total_drug_degree_dict[drug])
		mean_deg = statistics.mean(total_drug_degree_dict[drug])
		median_deg = statistics.median(total_drug_degree_dict[drug])
		sd_deg = statistics.pstdev(total_drug_degree_dict[drug])
		file_write.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (drug, nets_number, min_deg, max_deg, mean_deg, median_deg, sd_deg))


with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\Further_Analysis\\Drugs\\COVID19_GDDS_Results_Drugs_All_Nodes_Relative-withZeros.tsv', 'a+') as file_write:
	file_write.write("drug\t# of nets it is found in\tmin-degree-rel\tmax-degree-rel\tmean-degree-rel\tmedian-degree-rel\tsd-degree-rel\n")
	for drug in total_drug_normalized_degree_dict :
		nets_number = len(total_drug_normalized_degree_dict[drug])
		min_deg = min(total_drug_normalized_degree_dict[drug])
		max_deg = max(total_drug_normalized_degree_dict[drug])
		mean_deg = statistics.mean(total_drug_normalized_degree_dict[drug])
		median_deg = statistics.median(total_drug_normalized_degree_dict[drug])
		sd_deg = statistics.pstdev(total_drug_normalized_degree_dict[drug])
		file_write.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (drug, nets_number, min_deg, max_deg, mean_deg, median_deg, sd_deg))

### compute Z-score for the covid diseases
# For each drug found in covid network, compute 2 Z-scores.

with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\Further_Analysis\\Drugs\\COVID19_GDDS_Results_ForCovidDrugs-withZeros.tsv', 'a+') as file_write:
	file_write.write("covid drug\t# of nets it is found in\tmin-degree\tmax-degree\tmean-degree\tmedian-degree\tsd-degree\tcovid-degree\tZ-score\n")
	for drug in covid_drugs :
		nets_number = len(total_drug_degree_dict[drug])
		min_deg = min(total_drug_degree_dict[drug])
		max_deg = max(total_drug_degree_dict[drug])
		mean_deg = statistics.mean(total_drug_degree_dict[drug])
		median_deg = statistics.median(total_drug_degree_dict[drug])
		sd_deg = statistics.pstdev(total_drug_degree_dict[drug])
		covid_deg = covid_drug_degree_dict[drug]
		if sd_deg != 0 :
			z_score = ( covid_deg - mean_deg) / sd_deg
		else : 
			z_score = 'NaN'
		file_write.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (drug, nets_number, min_deg, max_deg, mean_deg, median_deg, sd_deg, covid_deg, z_score))


with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\Further_Analysis\\Drugs\\COVID19_GDDS_Results_ForCovidDrugs_Relative-withZeros.tsv', 'a+') as file_write:
	file_write.write("covid drug\t# of nets it is found in\tmin-degree-rel\tmax-degree-rel\tmean-degree-rel\tmedian-degree-rel\tsd-degree-rel\tcovid-degree-rel\tZ-score-rel\n")
	for drug in covid_drugs :
		nets_number = len(total_drug_normalized_degree_dict[drug])
		min_deg = min(total_drug_normalized_degree_dict[drug])
		max_deg = max(total_drug_normalized_degree_dict[drug])
		mean_deg = statistics.mean(total_drug_normalized_degree_dict[drug])
		median_deg = statistics.median(total_drug_normalized_degree_dict[drug])
		sd_deg = statistics.pstdev(total_drug_normalized_degree_dict[drug])
		covid_deg = covid_drug_norm_degree_dict[drug]
		if sd_deg != 0 :
			z_score = ( covid_deg - mean_deg) / sd_deg
		else : 
			z_score = 'NaN'
		file_write.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (drug, nets_number, min_deg, max_deg, mean_deg, median_deg, sd_deg, covid_deg, z_score))










stop = timeit.default_timer()
print(stop - start) 