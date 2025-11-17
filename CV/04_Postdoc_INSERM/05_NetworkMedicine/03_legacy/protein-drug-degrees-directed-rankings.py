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
#########################  Compute network Protein degree distribution #################################
##########################################  based on edges to proteins #####################################
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

proteins_degrees = nx.degree(G, nbunch=[n for n in G.nodes if descr_dict[n] == 'Human PPI (target)'])
# print("degree of protein nodes only", proteins_degrees)

average_degree_of_proteins = statistics.mean([ tpl[1] for tpl in proteins_degrees ])
print("average degree of proteins nodes", average_degree_of_proteins) # 

# Number of nodes: 16007
# Number of edges: 53917
# Average degree:   6.7367
# average degree of proteins nodes 115.66739606126914


# # # # #####################################################################################################
# # # # ################################# CALCULATION ON EDGES WITH FILTERING PROT-DRUG #####################
# # # # #####################################################################################################

# now we keep only the edges of proteins that interact with other proteins
prots_linked_to_drugs = []
for u,v,c in G.edges(data=True) :	
	if (descr_dict[u] == 'Human PPI (target)' and descr_dict[v] == 'Drug') :
		prots_linked_to_drugs.append(u)
print(len(prots_linked_to_drugs), len(set(prots_linked_to_drugs))) # 13018 169

# ({'Drug': 5703, 'Disease': 4176, 'GO': 3487, 'Symptom': 2157, 'Human PPI (target)': 457, 'Viral Gene': 27})

# print(sorted(prots_linked_to_drugs))

protein_counts = sorted(prots_linked_to_drugs)

protein_counts_dict = defaultdict( int )
for protein in protein_counts:
    protein_counts_dict[protein] += 1

my_proteins_nodes = [n for n in G.nodes if descr_dict[n] == 'Human PPI (target)']
for protein in my_proteins_nodes :
	if protein_counts_dict.get(protein) == None :
		# print("%s is missing in the linked proteins : Adding it with a frequency of 0.\n" % protein)
		protein_counts_dict[protein]= 0

average_degree_of_my_proteins_nodes = statistics.mean([ protein_counts_dict[protein] for protein in protein_counts_dict ])
print("average_degree_of_my_proteins_nodes", average_degree_of_my_proteins_nodes)
# 4.37417943107221

max_degree_of_my_proteins_nodes = max([ protein_counts_dict[protein] for protein in protein_counts_dict ])
print("max_degree_of_my_proteins_nodes", max_degree_of_my_proteins_nodes)
# average_degree_of_my_proteins_nodes 28.485776805251643
# max_degree_of_my_proteins_nodes 1105



ordered_protein_counts_dict = OrderedDict(sorted(protein_counts_dict.items()))


# with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\Further_Analysis\\Proteins\\COVID19_GDDS_proteins_degrees_to_drugs.tsv', 'a+') as file_write:
# 	for protein in ordered_protein_counts_dict :
# 		file_write.write("%s\t%s\t%s\n" % (protein, ordered_protein_counts_dict[protein], ordered_protein_counts_dict[protein] / sum(ordered_protein_counts_dict.values())))


# # # ##################################################################
# # # ############ SAME ANALYSIS FOR RANDOM NETWORKS ###################
# # # # ##################################################################

# with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\Further_Analysis\\Proteins\\Log_degree_analysis_on_random_networks_for_proteins_drugs.txt', 'a+') as log_write:

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
# 		################################# GET EDGES CONNECTING proteinS TO FIRST ORDER PROTEINS ################
# 		#####################################################################################################

# 		prots_linked_to_drugs = []
# 		for u,v,c in G.edges(data=True) :	
# 			if (descr_dict[u] == 'Human PPI (target)' and descr_dict[v] == 'Drug') :
# 				prots_linked_to_drugs.append(u)
		
# 		protein_counts = sorted(prots_linked_to_drugs)

# 		protein_counts_dict = defaultdict( int )
# 		for protein in protein_counts:
# 		    protein_counts_dict[protein] += 1

# 		my_proteins_nodes = [n for n in G.nodes if descr_dict[n] == 'Human PPI (target)']
# 		for protein in my_proteins_nodes :
# 			if protein_counts_dict.get(protein) == None :
# 				# print("%s is missing in the linked proteins : Adding it with a frequency of 0.\n" % protein)
# 				protein_counts_dict[protein]= 0

# 		average_degree_of_my_proteins_nodes = statistics.mean([ protein_counts_dict[protein] for protein in protein_counts_dict ])
# 		log_write.write("Average degree of proteins linked to drug(s) : %s\n" % average_degree_of_my_proteins_nodes)

# 		ordered_protein_counts_dict = OrderedDict(sorted(protein_counts_dict.items()))

# 		with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\Further_Analysis\\Proteins\\prot_drug\\COVID19_GDDS_proteins_degrees_to_drugs_%s.tsv' % str(iteration), 'a+') as file_write:
# 			for protein in ordered_protein_counts_dict :
# 				if (sum(ordered_protein_counts_dict.values()) != 0) :
# 					file_write.write("%s\t%s\t%s\n" % (protein, ordered_protein_counts_dict[protein], ordered_protein_counts_dict[protein] / sum(ordered_protein_counts_dict.values())))
# 				else :
# 					file_write.write("%s\t%s\t%s\n" % (protein, ordered_protein_counts_dict[protein], ordered_protein_counts_dict[protein]))



# # # # ##################################################################
# # # # ############ COMPUTE DEGREES MEANS AND Z-SCORES ##################
# # # # ######## INCLUDING ABSENCE OF drugs AS ZERO DEGREES ###########
# # # # ##################################################################


# # For covid network, compute the mean and the standard deviation of proteins degrees and normalized degrees
covid_proteins = []
covid_protein_degree_dict = {}
covid_protein_norm_degree_dict = {}

total_protein_degree_dict = defaultdict(list)
total_protein_normalized_degree_dict = defaultdict(list)

with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\Further_Analysis\\Proteins\\COVID19_GDDS_proteins_degrees_to_drugs.tsv', 'r') as data :
	for line in data :
		line = line.split('\t')
		protein = line[0]
		covid_proteins.append(protein)
		protein_degree = int(line[1])
		protein_norm_degree = float(line[2].replace('\n', ''))
		if covid_protein_degree_dict.get(protein) == None :
			covid_protein_degree_dict[protein] = protein_degree
		else :
			print("warning : %s duplicate ?!" % protein)
		covid_protein_norm_degree_dict[protein] = protein_norm_degree

		if total_protein_degree_dict.get(protein) == None :
			total_protein_degree_dict[protein] = [protein_degree]
		else :
			total_protein_degree_dict[protein].append(protein_degree)

		if total_protein_normalized_degree_dict.get(protein) == None :
			total_protein_normalized_degree_dict[protein] = [protein_norm_degree]
		else :
			total_protein_normalized_degree_dict[protein].append(protein_norm_degree)

mean_protein_deg = statistics.mean(covid_protein_degree_dict.values())
mean_norm_protein_deg = statistics.mean(covid_protein_norm_degree_dict.values())
sd_protein_deg = statistics.pstdev(covid_protein_degree_dict.values())
sd_norm_protein_deg = statistics.pstdev(covid_protein_norm_degree_dict.values())
print(mean_protein_deg, sd_protein_deg, mean_norm_protein_deg, sd_norm_protein_deg) #
print(len(covid_proteins),len(covid_protein_degree_dict),len(covid_protein_norm_degree_dict)) #
# 28.485776805251643 119.87329340458344 0.002188183807439825 0.00920827265360143
# 457 457 457



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

# nodefiles_numbers = sorted(nodefiles_numbers)
# edgefiles_numbers = sorted(edgefiles_numbers)


for iteration in nodefiles_numbers[1:] :
	# log_write.write("\n########################## ITERATION %s ###########################\n" % iteration)
	print("\n########################## ITERATION %s ###########################\n" % iteration)

	with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\Further_Analysis\\Proteins\\prot_drug\\COVID19_GDDS_proteins_degrees_to_drugs_%s.tsv' % str(iteration), 'r') as data: # Open the file  
		for line in data :
			line = line.split('\t')
			protein = line[0]
			degree = float(line[1])
			normalized_degree = float(line[2].replace('\n', ''))

			if total_protein_degree_dict.get(protein) == None :
				total_protein_degree_dict[protein] = [degree]
			else :
				total_protein_degree_dict[protein].append(degree)

			if total_protein_normalized_degree_dict.get(protein) == None :
				total_protein_normalized_degree_dict[protein] = [normalized_degree]
			else :
				total_protein_normalized_degree_dict[protein].append(normalized_degree)

print(len(total_protein_degree_dict)) # 19928

print(len(total_protein_normalized_degree_dict)) # 19928

### make the files with all calculations
with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\Further_Analysis\\Proteins\\COVID19_GDDS_Results_Proteins_to_drugs-withZeros.tsv', 'a+') as file_write:
	file_write.write("protein\t# of nets it is found in\tmin-degree\tmax-degree\tmean-degree\tmedian-degree\tsd-degree\n")
	for protein in total_protein_degree_dict :
		nets_number = len(total_protein_degree_dict[protein])
		min_deg = min(total_protein_degree_dict[protein])
		max_deg = max(total_protein_degree_dict[protein])
		mean_deg = statistics.mean(total_protein_degree_dict[protein])
		median_deg = statistics.median(total_protein_degree_dict[protein])
		sd_deg = statistics.pstdev(total_protein_degree_dict[protein])
		file_write.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (protein, nets_number, min_deg, max_deg, mean_deg, median_deg, sd_deg))


with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\Further_Analysis\\Proteins\\COVID19_GDDS_Results_Proteins_to_drugs_Relative-withZeros.tsv', 'a+') as file_write:
	file_write.write("protein\t# of nets it is found in\tmin-degree-rel\tmax-degree-rel\tmean-degree-rel\tmedian-degree-rel\tsd-degree-rel\n")
	for protein in total_protein_normalized_degree_dict :
		nets_number = len(total_protein_normalized_degree_dict[protein])
		min_deg = min(total_protein_normalized_degree_dict[protein])
		max_deg = max(total_protein_normalized_degree_dict[protein])
		mean_deg = statistics.mean(total_protein_normalized_degree_dict[protein])
		median_deg = statistics.median(total_protein_normalized_degree_dict[protein])
		sd_deg = statistics.pstdev(total_protein_normalized_degree_dict[protein])
		file_write.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (protein, nets_number, min_deg, max_deg, mean_deg, median_deg, sd_deg))

### compute Z-score for the covid proteins
# For each protein found in covid network, compute 2 Z-scores.

with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\Further_Analysis\\Proteins\\COVID19_GDDS_Results_ForCovidProteins_to_drugs-withZeros.tsv', 'a+') as file_write:
	file_write.write("covid protein\t# of nets it is found in\tmin-degree\tmax-degree\tmean-degree\tmedian-degree\tsd-degree\tcovid-degree\tZ-score\n")
	for protein in covid_proteins :
		nets_number = len(total_protein_degree_dict[protein])
		min_deg = min(total_protein_degree_dict[protein])
		max_deg = max(total_protein_degree_dict[protein])
		mean_deg = statistics.mean(total_protein_degree_dict[protein])
		median_deg = statistics.median(total_protein_degree_dict[protein])
		sd_deg = statistics.pstdev(total_protein_degree_dict[protein])
		covid_deg = covid_protein_degree_dict[protein]
		if sd_deg != 0 :
			z_score = ( covid_deg - mean_deg) / sd_deg
		else : 
			z_score = 'NaN'
		file_write.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (protein, nets_number, min_deg, max_deg, mean_deg, median_deg, sd_deg, covid_deg, z_score))


with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\Further_Analysis\\Proteins\\COVID19_GDDS_Results_ForCovidProteins_to_drugs_Relative-withZeros.tsv', 'a+') as file_write:
	file_write.write("covid protein\t# of nets it is found in\tmin-degree-rel\tmax-degree-rel\tmean-degree-rel\tmedian-degree-rel\tsd-degree-rel\tcovid-degree-rel\tZ-score-rel\n")
	for protein in covid_proteins :
		nets_number = len(total_protein_normalized_degree_dict[protein])
		min_deg = min(total_protein_normalized_degree_dict[protein])
		max_deg = max(total_protein_normalized_degree_dict[protein])
		mean_deg = statistics.mean(total_protein_normalized_degree_dict[protein])
		median_deg = statistics.median(total_protein_normalized_degree_dict[protein])
		sd_deg = statistics.pstdev(total_protein_normalized_degree_dict[protein])
		covid_deg = covid_protein_norm_degree_dict[protein]
		if sd_deg != 0 :
			z_score = ( covid_deg - mean_deg) / sd_deg
		else : 
			z_score = 'NaN'
		file_write.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (protein, nets_number, min_deg, max_deg, mean_deg, median_deg, sd_deg, covid_deg, z_score))










stop = timeit.default_timer()
print(stop - start) 