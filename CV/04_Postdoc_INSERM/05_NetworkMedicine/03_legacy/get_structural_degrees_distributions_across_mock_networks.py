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



with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\COVID19_GDDS_nodes.csv', 'r') as nodecsv: # Open the file                       
		nodereader = csv.reader(nodecsv) # Read the csv  
		nodes = [n for n in nodereader][1:]                     
		node_names = [n[0] for n in nodes] # Get a list of only the node names 

descr_dict_covid = {}
for node in nodes: 
	descr_dict_covid[node[0]] = node[2]

# covid_proteins =  [node for node in descr_dict_covid if descr_dict_covid[node] == 'Human PPI (target)']
# covid_GO = [node for node in descr_dict_covid if descr_dict_covid[node] == 'GO']
# covid_symptom = [node for node in descr_dict_covid if descr_dict_covid[node] == 'Symptom']
# covid_disease = [node for node in descr_dict_covid if descr_dict_covid[node] == 'Disease']
covid_drug = [node for node in descr_dict_covid if descr_dict_covid[node] == 'Drug']


nodefiles_numbers = []
edgefiles_numbers = []
for filename in os.listdir("A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\random_networks\\Mock_networks_21Apr"):
	if 'nodes_' in filename :
		filename = filename.split(".csv")[0]
		nodefiles_numbers.append(int(filename[19:]))
	if 'edges_' in filename :
		filename = filename.split(".csv")[0]
		edgefiles_numbers.append(int(filename[19:]))

nodefiles_numbers = sorted(nodefiles_numbers)
edgefiles_numbers = sorted(edgefiles_numbers)

networks_to_avoid = [3, 11, 22, 23, 33, 34, 38, 57, 72, 77, 123, 127, 148, 156, 157, 183, 205, 209, 210, 217, 232, 252, 268, 281, 291, 312, 341, 353, 362, 364, 413, 421, 433, 443, 451, 494, 525, 532, 569, 592, 629, 658, 662, 674, 676, 701, 705, 712, 736, 757, 788, 796, 825, 827, 838, 873, 876, 908, 910, 916, 917, 918, 938, 958, 965, 1088, 1125, 1128, 1145, 1146, 1165, 1171, 1183, 1185, 1192, 1195, 1211, 1244, 1254, 1257, 1300, 1353, 1354, 1364, 1372, 1394, 1395, 1401, 1410, 1427, 1454, 1477, 1497, 1517, 1532, 1549, 1551, 1575, 1578, 1611, 1624, 1629, 1637, 1650, 1655, 1697, 1700, 1713, 1736, 1745, 1762, 1790, 1800, 1815, 1824, 1827, 1866, 1872, 1886, 1887, 2053, 2072, 2106, 2115, 2128, 2160, 2217, 2343, 2383, 2388, 2413, 2458, 2465]
nodefiles_numbers = [item for item in nodefiles_numbers if item not in networks_to_avoid]
print(nodefiles_numbers)


### for covid proteins
# covid_proteins_structural_degree_distribution,covid_proteins_structural_strength_distribution = {},{} 
# for iteration in nodefiles_numbers[1:] :
# 	print("\n########################## ITERATION %s ###########################\n" % iteration)
# 	with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\Further_Analysis\\full_analysis_for_source_to_proteins_undirected\\Human PPI (target)-prot\\COVID19_GDDS_Human PPI (target)_degrees_to_proteins_%s.tsv' % str(iteration), 'r') as file_read :                 
# 		data = file_read.readlines()
# 		for line in data :
# 			line = line.split('\t')
# 			protein = line[0]
# 			if protein in covid_proteins :
# 				structural_degree = int(line[1])
# 				structural_strength = float(line[2].replace('\n',''))
# 				if protein not in covid_proteins_structural_degree_distribution :
# 					covid_proteins_structural_degree_distribution[protein] = []
# 				if protein not in covid_proteins_structural_strength_distribution :
# 					covid_proteins_structural_strength_distribution[protein] = []
# 				covid_proteins_structural_degree_distribution[protein].append(structural_degree)
# 				covid_proteins_structural_strength_distribution[protein].append(structural_strength)
# with open('pvalues_from_zscores\\covid_proteins_structural_degree_distribution_in_mock_networks.tsv', 'a+') as file_write :
# 	file_write.write('covid_protein\tstructural_degree_in_mock_networks\n')
# 	for protein in covid_proteins_structural_degree_distribution :
# 		file_write.write('%s\t' % protein)
# 		for degree in covid_proteins_structural_degree_distribution[protein][:-1] :
# 			file_write.write('%s\t' % degree)
# 		file_write.write('%s\n' % covid_proteins_structural_degree_distribution[protein][-1])
# with open('pvalues_from_zscores\\covid_proteins_structural_strength_distribution_in_mock_networks.tsv', 'a+') as file_write :
# 	file_write.write('covid_protein\tstructural_strength_in_mock_networks\n')
# 	for protein in covid_proteins_structural_strength_distribution :
# 		file_write.write('%s\t' % protein)
# 		for degree in covid_proteins_structural_strength_distribution[protein][:-1] :
# 			file_write.write('%s\t' % degree)
# 		file_write.write('%s\n' % covid_proteins_structural_strength_distribution[protein][-1])


### for covid GO
# covid_GO_degree_to_proteins_distribution,covid_GO_structural_strength_to_proteins_distribution = {},{} 
# for iteration in nodefiles_numbers[1:] :
# 	print("\n########################## ITERATION %s ###########################\n" % iteration)
# 	with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\Further_Analysis\\full_analysis_for_source_to_proteins_undirected\\GO-prot\\COVID19_GDDS_GO_degrees_to_proteins_%s.tsv' % str(iteration), 'r') as file_read :                 
# 		data = file_read.readlines()
# 		for line in data :
# 			line = line.split('\t')
# 			GO = line[0]
# 			if GO in covid_GO :
# 				degree_to_proteins = int(line[1])
# 				structural_strength_to_proteins = float(line[2].replace('\n',''))
# 				if GO  not in covid_GO_degree_to_proteins_distribution :
# 					covid_GO_degree_to_proteins_distribution[GO] = []
# 				if GO not in covid_GO_structural_strength_to_proteins_distribution :
# 					covid_GO_structural_strength_to_proteins_distribution[GO] = []
# 				covid_GO_degree_to_proteins_distribution[GO].append(degree_to_proteins)
# 				covid_GO_structural_strength_to_proteins_distribution[GO].append(structural_strength_to_proteins)
# with open('pvalues_from_zscores\\covid_GO_degree_to_proteins_distribution_in_mock_networks.tsv', 'a+') as file_write :
# 	file_write.write('covid_GO\tdegree_to_proteins_in_mock_networks\n')
# 	for GO in covid_GO_degree_to_proteins_distribution :
# 		file_write.write('%s\t' % GO)
# 		for degree in covid_GO_degree_to_proteins_distribution[GO][:-1] :
# 			file_write.write('%s\t' % degree)
# 		file_write.write('%s\n' % covid_GO_degree_to_proteins_distribution[GO][-1])
# with open('pvalues_from_zscores\\covid_GO_structural_strength_to_proteins_distribution_in_mock_networks.tsv', 'a+') as file_write :
# 	file_write.write('covid_GO\tstructural_strength_to_proteins_in_mock_networks\n')
# 	for GO in covid_GO_structural_strength_to_proteins_distribution :
# 		file_write.write('%s\t' % GO)
# 		for degree in covid_GO_structural_strength_to_proteins_distribution[GO][:-1] :
# 			file_write.write('%s\t' % degree)
# 		file_write.write('%s\n' % covid_GO_structural_strength_to_proteins_distribution[GO][-1])


# ####for covid symptom
# covid_symptom_degree_to_proteins_distribution,covid_symptom_structural_strength_to_proteins_distribution = {},{} 
# for iteration in nodefiles_numbers[1:] :
# 	print("\n########################## ITERATION %s ###########################\n" % iteration)
# 	with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\Further_Analysis\\full_analysis_for_source_to_proteins_undirected\\Symptom-prot\\COVID19_GDDS_symptom_degrees_to_proteins_%s.tsv' % str(iteration), 'r') as file_read :                 
# 		data = file_read.readlines()
# 		for line in data :
# 			line = line.split('\t')
# 			symptom = line[0]
# 			if symptom in covid_symptom :
# 				degree_to_proteins = int(line[1])
# 				structural_strength_to_proteins = float(line[2].replace('\n',''))
# 				if symptom  not in covid_symptom_degree_to_proteins_distribution :
# 					covid_symptom_degree_to_proteins_distribution[symptom] = []
# 				if symptom not in covid_symptom_structural_strength_to_proteins_distribution :
# 					covid_symptom_structural_strength_to_proteins_distribution[symptom] = []
# 				covid_symptom_degree_to_proteins_distribution[symptom].append(degree_to_proteins)
# 				covid_symptom_structural_strength_to_proteins_distribution[symptom].append(structural_strength_to_proteins)
# with open('pvalues_from_zscores\\covid_symptom_degree_to_proteins_distribution_in_mock_networks.tsv', 'a+') as file_write :
# 	file_write.write('covid_symptom\tdegree_to_proteins_in_mock_networks\n')
# 	for symptom in covid_symptom_degree_to_proteins_distribution :
# 		file_write.write('%s\t' % symptom)
# 		for degree in covid_symptom_degree_to_proteins_distribution[symptom][:-1] :
# 			file_write.write('%s\t' % degree)
# 		file_write.write('%s\n' % covid_symptom_degree_to_proteins_distribution[symptom][-1])
# with open('pvalues_from_zscores\\covid_symptom_structural_strength_to_proteins_distribution_in_mock_networks.tsv', 'a+') as file_write :
# 	file_write.write('covid_symptom\tstructural_strength_to_proteins_in_mock_networks\n')
# 	for symptom in covid_symptom_structural_strength_to_proteins_distribution :
# 		file_write.write('%s\t' % symptom)
# 		for degree in covid_symptom_structural_strength_to_proteins_distribution[symptom][:-1] :
# 			file_write.write('%s\t' % degree)
# 		file_write.write('%s\n' % covid_symptom_structural_strength_to_proteins_distribution[symptom][-1])


# ## for covid disease		
# covid_disease_degree_to_proteins_distribution,covid_disease_structural_strength_to_proteins_distribution = {},{} 
# for iteration in nodefiles_numbers[1:] :
# 	print("\n########################## ITERATION %s ###########################\n" % iteration)
# 	with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\Further_Analysis\\full_analysis_for_source_to_proteins_undirected\\Disease-prot\\COVID19_GDDS_Disease_degrees_to_proteins_%s.tsv' % str(iteration), 'r') as file_read :                 
# 		data = file_read.readlines()
# 		for line in data :
# 			line = line.split('\t')
# 			disease = line[0]
# 			if disease in covid_disease :
# 				degree_to_proteins = int(line[1])
# 				structural_strength_to_proteins = float(line[2].replace('\n',''))
# 				if disease  not in covid_disease_degree_to_proteins_distribution :
# 					covid_disease_degree_to_proteins_distribution[disease] = []
# 				if disease not in covid_disease_structural_strength_to_proteins_distribution :
# 					covid_disease_structural_strength_to_proteins_distribution[disease] = []
# 				covid_disease_degree_to_proteins_distribution[disease].append(degree_to_proteins)
# 				covid_disease_structural_strength_to_proteins_distribution[disease].append(structural_strength_to_proteins)
# with open('pvalues_from_zscores\\covid_disease_degree_to_proteins_distribution_in_mock_networks.tsv', 'a+') as file_write :
# 	file_write.write('covid_disease\tdegree_to_proteins_in_mock_networks\n')
# 	for disease in covid_disease_degree_to_proteins_distribution :
# 		file_write.write('%s\t' % disease)
# 		for degree in covid_disease_degree_to_proteins_distribution[disease][:-1] :
# 			file_write.write('%s\t' % degree)
# 		file_write.write('%s\n' % covid_disease_degree_to_proteins_distribution[disease][-1])
# with open('pvalues_from_zscores\\covid_disease_structural_strength_to_proteins_distribution_in_mock_networks.tsv', 'a+') as file_write :
# 	file_write.write('covid_disease\tstructural_strength_to_proteins_in_mock_networks\n')
# 	for disease in covid_disease_structural_strength_to_proteins_distribution :
# 		file_write.write('%s\t' % disease)
# 		for degree in covid_disease_structural_strength_to_proteins_distribution[disease][:-1] :
# 			file_write.write('%s\t' % degree)
# 		file_write.write('%s\n' % covid_disease_structural_strength_to_proteins_distribution[disease][-1])




## for covid drug		
covid_drug_degree_to_proteins_distribution,covid_drug_structural_strength_to_proteins_distribution = {},{} 
for iteration in nodefiles_numbers[1:] :
	print("\n########################## ITERATION %s ###########################\n" % iteration)
	with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\Further_Analysis\\full_analysis_for_source_to_proteins_undirected\\Drug-prot\\COVID19_GDDS_Drug_degrees_to_proteins_%s.tsv' % str(iteration), 'r') as file_read :                 
		data = file_read.readlines()
		for line in data :
			line = line.split('\t')
			drug = line[0]
			if drug in covid_drug :
				degree_to_proteins = int(line[1])
				structural_strength_to_proteins = float(line[2].replace('\n',''))
				if drug  not in covid_drug_degree_to_proteins_distribution :
					covid_drug_degree_to_proteins_distribution[drug] = []
				if drug not in covid_drug_structural_strength_to_proteins_distribution :
					covid_drug_structural_strength_to_proteins_distribution[drug] = []
				covid_drug_degree_to_proteins_distribution[drug].append(degree_to_proteins)
				covid_drug_structural_strength_to_proteins_distribution[drug].append(structural_strength_to_proteins)
with open('covid_drug_degree_to_proteins_distribution_in_mock_networks.tsv', 'a+') as file_write :
	file_write.write('covid_drug\tdegree_to_proteins_in_mock_networks\n')
	for drug in covid_drug_degree_to_proteins_distribution :
		file_write.write('%s\t' % drug)
		for degree in covid_drug_degree_to_proteins_distribution[drug][:-1] :
			file_write.write('%s\t' % degree)
		file_write.write('%s\n' % covid_drug_degree_to_proteins_distribution[drug][-1])
with open('covid_drug_structural_strength_to_proteins_distribution_in_mock_networks.tsv', 'a+') as file_write :
	file_write.write('covid_drug\tstructural_strength_to_proteins_in_mock_networks\n')
	for drug in covid_drug_structural_strength_to_proteins_distribution :
		file_write.write('%s\t' % drug)
		for degree in covid_drug_structural_strength_to_proteins_distribution[drug][:-1] :
			file_write.write('%s\t' % degree)
		file_write.write('%s\n' % covid_drug_structural_strength_to_proteins_distribution[drug][-1])







# # dirname = "full_analysis_for_source_to_proteins_undirected"
# # os.mkdir(dirname)
# source_list = ['Symptom', 'Disease','Drug', 'Human PPI (target)']
# # source_list = ['GO']

# with open('Log_full_analysis_for_source_to_proteins_undirected.txt', 'a+') as full_analysis_log_write :

# 	#####################################################################################################
# 	#########################  Compute network source to protein degree distribution ####################
# 	##########################################  based on edges to protein ###############################
# 	#####################################################################################################

# 	with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\COVID19_GDDS_nodes.csv', 'r') as nodecsv: # Open the file                       
# 		nodereader = csv.reader(nodecsv) # Read the csv  
# 		nodes = [n for n in nodereader][1:]                     
# 		node_names = [n[0] for n in nodes] # Get a list of only the node names 

# 	with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\COVID19_GDDS_edges.csv', 'r') as edgecsv: # Open the file
# 		edgereader = csv.reader(edgecsv) # Read the csv     
# 		edges = [tuple(e) for e in edgereader][1:] # Retrieve the data

# 	G_covid = nx.Graph()
# 	G_covid.add_nodes_from(node_names)
# 	G_covid.add_edges_from(edges)

# 	descr_dict_covid = {}
# 	for node in nodes: 
# 		descr_dict_covid[node[0]] = node[2]

# 	print(nx.info(G_covid))
# 	# full_analysis_log_write.write(nx.info(G_covid))
# 	# Number of nodes: 16007
# 	# Number of edges: 53917
# 	# Average degree:   6.7367

# 	# print(nx.degree(G_covid)) # this is the list containing as many tuples as nodes, indicating each node's degree.

# 	average_degree = statistics.mean([ tpl[1] for tpl in nx.degree(G_covid) ])
# 	# print(average_degree) # 6.736677703504717

# 	for source_node in source_list :

# 		source_degrees = nx.degree(G_covid, nbunch=[n for n in G_covid.nodes if descr_dict_covid[n] == source_node])
# 		# print("degree of source nodes only", source_degrees)

# 		average_degree_of_sources = statistics.mean([ tpl[1] for tpl in source_degrees ])
# 		print("average degree of %s nodes" % source_node, average_degree_of_sources) # 
# 		# full_analysis_log_write.write("average degree of %s nodes : %s\n" % (source_node, average_degree_of_sources))
# 	# average degree of GO nodes 2.64123888729567
# 	# average degree of Symptom nodes 6.644413537320353
# 	# average degree of Disease nodes 4.329980842911877
# 	# average degree of Drug nodes 2.2826582500438364
# 	# average degree of Human PPI (target) nodes 115.66739606126914



# 	# # # # # #####################################################################################################
# 	# # # # # ################################# CALCULATION ON EDGES WITH FILTERING source - PROT #####################
# 	# # # # # #####################################################################################################

# 	for source_node in source_list :
# 		subfoldername = 'full_analysis_for_source_to_proteins_undirected/%s-prot' % (source_node)
# 		os.mkdir(subfoldername)
# 		print("Searching %s - Human PPI (target) links\n" % source_node)
# 		full_analysis_log_write.write("Searching %s - Human PPI (target) links\n" % source_node)
		
# 		# now we keep only the edges of source_node that interact with proteins
# 		sources_linked_to_proteins = []
# 		for u,v,c in G_covid.edges(data=True) :	
# 			if (descr_dict_covid[u] == source_node and descr_dict_covid[v] == 'Human PPI (target)') :
# 				sources_linked_to_proteins.append(u)
# 			if (descr_dict_covid[v] == source_node and descr_dict_covid[u] == 'Human PPI (target)') :
# 				sources_linked_to_proteins.append(v)	
# 		print("len(%s_linked_to_proteins), len(set(%s_linked_to_proteins))" % (source_node, source_node), len(sources_linked_to_proteins), len(set(sources_linked_to_proteins))) #
# 		full_analysis_log_write.write("len(%s_linked_to_proteins), len(set(%s_linked_to_proteins)) : %s, %s\n" % (source_node, source_node, len(sources_linked_to_proteins), len(set(sources_linked_to_proteins)))) #


# 		source_counts = sorted(sources_linked_to_proteins)

# 		source_counts_dict = defaultdict( int )
# 		for source in  source_counts:
# 		    source_counts_dict[source] += 1

# 		my_sources_nodes = [n for n in G_covid.nodes if descr_dict_covid[n] == source_node]
# 		for source in  my_sources_nodes :
# 			if source_counts_dict.get(source) == None :
# 				source_counts_dict[source]= 0

# 		average_degree_of_my_sources_nodes = statistics.mean([ source_counts_dict[source] for source in  source_counts_dict ])
# 		print("average_degree_of_my_sources_nodes", average_degree_of_my_sources_nodes)
# 		full_analysis_log_write.write("average_degree_of_my_%s_nodes : %s\n" % (source_node, average_degree_of_my_sources_nodes))

# 		max_degree_of_my_sources_nodes = max([ source_counts_dict[source] for source in  source_counts_dict ])
# 		print("max_degree_of_my_sources_nodes", max_degree_of_my_sources_nodes)

# 		ordered_source_counts_dict = OrderedDict(sorted(source_counts_dict.items()))


# 		with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\Further_Analysis\\full_analysis_for_source_to_proteins_undirected\\COVID19_GDDS_%s_degrees_to_proteins.tsv' % source_node, 'a+') as file_write:
# 			for source in ordered_source_counts_dict :
# 				file_write.write("%s\t%s\t%s\n" % (source, ordered_source_counts_dict[source], ordered_source_counts_dict[source] / sum(ordered_source_counts_dict.values())))


# 		# # # # # ##################################################################
# 		# # # # # ############ SAME ANALYSIS FOR RANDOM NETWORKS ###################
# 		# # # # # # ##################################################################
		
# 		with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\Further_Analysis\\full_analysis_for_source_to_proteins_undirected\\Log_degree_analysis_on_random_networks_for_%s_proteins.txt' % source_node, 'a+') as log_write:
# 			log_write.write("In covid net : len(sources_linked_to_proteins) and len(set(sources_linked_to_proteins)) are : %s and %s \n" % (len(sources_linked_to_proteins), len(set(sources_linked_to_proteins))) ) # # 9208 435
# 			log_write.write("In covid net : average_degree_of_my_sources_nodes = %s\n" % average_degree_of_my_sources_nodes)
# 			log_write.write("In covid net : max_degree_of_my_sources_nodes = %s\n" % max_degree_of_my_sources_nodes)

# 			nodefiles_numbers = []
# 			edgefiles_numbers = []
# 			for filename in os.listdir("A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\random_networks\\Mock_networks_21Apr"):
# 				if 'nodes_' in filename :
# 					filename = filename.split(".csv")[0]
# 					nodefiles_numbers.append(int(filename[19:]))
# 				if 'edges_' in filename :
# 					filename = filename.split(".csv")[0]
# 					edgefiles_numbers.append(int(filename[19:]))

# 			# log_write.write("len(nodefiles_numbers) = %s\n" % len(nodefiles_numbers)) #2051
# 			# log_write.write("len(edgefiles_numbers) = %s\n" % len(edgefiles_numbers)) #2051
# 			# log_write.write("difference between edgefiles_numbers and nodefiles_numbers = %s\n" % len(set(nodefiles_numbers).difference(set(edgefiles_numbers)))) #0

# 			nodefiles_numbers = sorted(nodefiles_numbers)
# 			edgefiles_numbers = sorted(edgefiles_numbers)

# 			networks_to_avoid = [3, 11, 22, 23, 33, 34, 38, 57, 72, 77, 123, 127, 148, 156, 157, 183, 205, 209, 210, 217, 232, 252, 268, 281, 291, 312, 341, 353, 362, 364, 413, 421, 433, 443, 451, 494, 525, 532, 569, 592, 629, 658, 662, 674, 676, 701, 705, 712, 736, 757, 788, 796, 825, 827, 838, 873, 876, 908, 910, 916, 917, 918, 938, 958, 965, 1088, 1125, 1128, 1145, 1146, 1165, 1171, 1183, 1185, 1192, 1195, 1211, 1244, 1254, 1257, 1300, 1353, 1354, 1364, 1372, 1394, 1395, 1401, 1410, 1427, 1454, 1477, 1497, 1517, 1532, 1549, 1551, 1575, 1578, 1611, 1624, 1629, 1637, 1650, 1655, 1697, 1700, 1713, 1736, 1745, 1762, 1790, 1800, 1815, 1824, 1827, 1866, 1872, 1886, 1887, 2053, 2072, 2106, 2115, 2128, 2160, 2217, 2343, 2383, 2388, 2413, 2458, 2465]
# 			nodefiles_numbers = [item for item in nodefiles_numbers if item not in networks_to_avoid]

# 			for iteration in nodefiles_numbers[1:] :
# 				print("\n########################## ITERATION %s -- %s ###########################\n" % (iteration, source_node))
# 				log_write.write("\n########################## ITERATION %s ###########################\n" % iteration)
			
# 				with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\random_networks\\Mock_networks_21Apr\\COVID19_GDDS_nodes_%s.csv' % str(iteration), 'r') as nodecsv: # Open the file                   
# 					nodereader = csv.reader(nodecsv) # Read the csv  
# 					nodes = [n for n in nodereader][1:]                     
# 					node_names = [n[0] for n in nodes] # Get a list of only the node names 

# 				with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\random_networks\\Mock_networks_21Apr\\COVID19_GDDS_edges_%s.csv' % str(iteration), 'r') as edgecsv: # Open the file
# 					edgereader = csv.reader(edgecsv) # Read the csv     
# 					edges = [tuple(e) for e in edgereader][1:] # Retrieve the data

# 					log_write.write("len(node_names) : %s\n" % len(node_names)) #16007 --> 13448
# 					log_write.write("len(edges) : %s\n" % len(edges)) #53917 --> 37835

# 				G_random = nx.Graph()
# 				G_random.add_nodes_from(node_names)
# 				G_random.add_edges_from(edges)

# 				descr_dict_random = {}
# 				description_random_set=set()
# 				for node in nodes: 
# 					descr_dict_random[node[0]] = node[2]
# 					description_random_set.add(node[2])

# 				# print(nx.info(G_random))
# 				log_write.write("\nRANDOM GRAPH %s\n" % iteration)
# 				log_write.write("%s\n" % nx.info(G_random))

# 				#####################################################################################################
# 				################################# GET EDGES CONNECTING sources TO proteins ################
# 				#####################################################################################################

# 				sources_linked_to_proteins = []
# 				for u,v,c in G_random.edges(data=True) :	
# 					if (descr_dict_random[u] == source_node and descr_dict_random[v] == 'Human PPI (target)') :
# 						sources_linked_to_proteins.append(u)
# 					if (descr_dict_random[v] == source_node and descr_dict_random[u] == 'Human PPI (target)') :
# 						sources_linked_to_proteins.append(v)		
				
# 				source_counts = sorted(sources_linked_to_proteins)

# 				source_counts_dict = defaultdict( int )
# 				for source in  source_counts:
# 				    source_counts_dict[source] += 1

# 				my_sources_nodes = [n for n in G_random.nodes if descr_dict_random[n] == source]
# 				for source in  my_sources_nodes :
# 					if source_counts_dict.get(source) == None :
# 						# print("%s is missing in the linked sources : Adding it with a frequency of 0.\n" % protein)
# 						source_counts_dict[source]= 0

# 				average_degree_of_my_sources_nodes = statistics.mean([ source_counts_dict[source] for source in  source_counts_dict ])
# 				log_write.write("Average degree of %s linked to Human PPI (target) : %s\n" % (source_node, average_degree_of_my_sources_nodes))

# 				ordered_source_counts_dict = OrderedDict(sorted(source_counts_dict.items()))

# 				with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\Further_Analysis\\full_analysis_for_source_to_proteins_undirected\\%s-prot\\COVID19_GDDS_%s_degrees_to_proteins_%s.tsv' % (source_node, source_node, str(iteration)), 'a+') as file_write:
# 					for source in ordered_source_counts_dict :
# 						if (sum(ordered_source_counts_dict.values()) != 0) :
# 							file_write.write("%s\t%s\t%s\n" % (source, ordered_source_counts_dict[source], ordered_source_counts_dict[source] / sum(ordered_source_counts_dict.values())))
# 						else :
# 							file_write.write("%s\t%s\t%s\n" % (source, ordered_source_counts_dict[source], ordered_source_counts_dict[source]))



# 		# # # # ##################################################################
# 		# # # # ############ COMPUTE DEGREES MEANS AND Z-SCORES ##################
# 		# # # # ######## INCLUDING ABSENCE OF source AS ZERO DEGREES ###########
# 		# # # # ##################################################################


# 		# # For covid network, compute the mean and the standard deviation of proteins degrees and normalized degrees
# 		covid_sources = []
# 		covid_source_degree_dict = {}
# 		covid_source_norm_degree_dict = {}

# 		total_source_degree_dict = defaultdict(list)
# 		total_source_normalized_degree_dict = defaultdict(list)

# 		with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\Further_Analysis\\full_analysis_for_source_to_proteins_undirected\\COVID19_GDDS_%s_degrees_to_proteins.tsv' % source_node, 'r') as data :
# 			for line in data :
# 				line = line.split('\t')
# 				source = line[0]
# 				covid_sources.append(source)
# 				source_degree = int(line[1])
# 				source_norm_degree = float(line[2].replace('\n', ''))
# 				if covid_source_degree_dict.get(source) == None :
# 					covid_source_degree_dict[source] = source_degree
# 				else :
# 					print("warning : %s duplicate ?!" % protein)
# 				covid_source_norm_degree_dict[source] = source_norm_degree

# 				if total_source_degree_dict.get(source) == None :
# 					total_source_degree_dict[source] = [source_degree]
# 				else :
# 					total_source_degree_dict[source].append(source_degree)

# 				if total_source_normalized_degree_dict.get(source) == None :
# 					total_source_normalized_degree_dict[source] = [source_norm_degree]
# 				else :
# 					total_source_normalized_degree_dict[source].append(source_norm_degree)

# 		mean_source_deg = statistics.mean(covid_source_degree_dict.values())
# 		mean_norm_source_deg = statistics.mean(covid_source_norm_degree_dict.values())
# 		sd_source_deg = statistics.pstdev(covid_source_degree_dict.values())
# 		sd_norm_source_deg = statistics.pstdev(covid_source_norm_degree_dict.values())
# 		print(mean_source_deg, sd_source_deg, mean_norm_source_deg, sd_norm_source_deg) #
# 		print(len(covid_sources),len(covid_source_degree_dict),len(covid_source_norm_degree_dict)) #
# 		# 20.14879649890591 18.9696222570052 0.002188183807439825 0.00206012405050013
# 		# 457 457 457






# 		nodefiles_numbers = []
# 		edgefiles_numbers = []
# 		for filename in os.listdir("A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\random_networks\\Mock_networks_21Apr"):
# 			if 'nodes_' in filename :
# 				filename = filename.split(".csv")[0]
# 				nodefiles_numbers.append(int(filename[19:]))
# 			if 'edges_' in filename :
# 				filename = filename.split(".csv")[0]
# 				edgefiles_numbers.append(int(filename[19:]))

# 		# log_write.write("len(nodefiles_numbers) = %s\n" % len(nodefiles_numbers)) #2051
# 		# log_write.write("len(edgefiles_numbers) = %s\n" % len(edgefiles_numbers)) #2051
# 		# log_write.write("difference between edgefiles_numbers and nodefiles_numbers = %s\n" % len(set(nodefiles_numbers).difference(set(edgefiles_numbers)))) #0
# 		print("len(nodefiles_numbers) = %s\n" % len(nodefiles_numbers)) #2051
# 		print("len(edgefiles_numbers) = %s\n" % len(edgefiles_numbers)) #2051
# 		print("difference between edgefiles_numbers and nodefiles_numbers = %s\n" % len(set(nodefiles_numbers).difference(set(edgefiles_numbers)))) #0

# 		# nodefiles_numbers = sorted(nodefiles_numbers)
# 		# edgefiles_numbers = sorted(edgefiles_numbers)

# 		networks_to_avoid = [3, 11, 22, 23, 33, 34, 38, 57, 72, 77, 123, 127, 148, 156, 157, 183, 205, 209, 210, 217, 232, 252, 268, 281, 291, 312, 341, 353, 362, 364, 413, 421, 433, 443, 451, 494, 525, 532, 569, 592, 629, 658, 662, 674, 676, 701, 705, 712, 736, 757, 788, 796, 825, 827, 838, 873, 876, 908, 910, 916, 917, 918, 938, 958, 965, 1088, 1125, 1128, 1145, 1146, 1165, 1171, 1183, 1185, 1192, 1195, 1211, 1244, 1254, 1257, 1300, 1353, 1354, 1364, 1372, 1394, 1395, 1401, 1410, 1427, 1454, 1477, 1497, 1517, 1532, 1549, 1551, 1575, 1578, 1611, 1624, 1629, 1637, 1650, 1655, 1697, 1700, 1713, 1736, 1745, 1762, 1790, 1800, 1815, 1824, 1827, 1866, 1872, 1886, 1887, 2053, 2072, 2106, 2115, 2128, 2160, 2217, 2343, 2383, 2388, 2413, 2458, 2465]
# 		nodefiles_numbers = [item for item in nodefiles_numbers if item not in networks_to_avoid]

# 		for iteration in nodefiles_numbers[1:] :
# 			# log_write.write("\n########################## ITERATION %s ###########################\n" % iteration)
# 			print("\n########################## ITERATION %s ###########################\n" % iteration)

# 			with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\Further_Analysis\\full_analysis_for_source_to_proteins_undirected\\%s-prot\\COVID19_GDDS_%s_degrees_to_proteins_%s.tsv' % (source_node, source_node, str(iteration)), 'r') as data: # Open the file  
# 				for line in data :
# 					line = line.split('\t')
# 					source = line[0]
# 					degree = float(line[1])
# 					normalized_degree = float(line[2].replace('\n', ''))

# 					if total_source_degree_dict.get(source) == None :
# 						total_source_degree_dict[source] = [degree]
# 					else :
# 						total_source_degree_dict[source].append(degree)

# 					if total_source_normalized_degree_dict.get(source) == None :
# 						total_source_normalized_degree_dict[source] = [normalized_degree]
# 					else :
# 						total_source_normalized_degree_dict[source].append(normalized_degree)

# 		print(len(total_source_degree_dict)) # 19928

# 		print(len(total_source_normalized_degree_dict)) # 19928

# 		### make the files with all calculations
# 		with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\Further_Analysis\\full_analysis_for_source_to_proteins_undirected\\COVID19_GDDS_Results_%s_to_Proteins-withZeros.tsv' % source_node, 'a+') as file_write:
# 			file_write.write("%s\t# of nets it is found in\tmin-degree\tmax-degree\tmean-degree\tmedian-degree\tsd-degree\n" % source_node)
# 			for source in  total_source_degree_dict :
# 				nets_number = len(total_source_degree_dict[source])
# 				min_deg = min(total_source_degree_dict[source])
# 				max_deg = max(total_source_degree_dict[source])
# 				mean_deg = statistics.mean(total_source_degree_dict[source])
# 				median_deg = statistics.median(total_source_degree_dict[source])
# 				sd_deg = statistics.pstdev(total_source_degree_dict[source])
# 				file_write.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (source, nets_number, min_deg, max_deg, mean_deg, median_deg, sd_deg))


# 		with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\Further_Analysis\\full_analysis_for_source_to_proteins_undirected\\COVID19_GDDS_Results_%s_to_Proteins_Relative-withZeros.tsv' % source_node, 'a+') as file_write:
# 			file_write.write("%s\t# of nets it is found in\tmin-degree-rel\tmax-degree-rel\tmean-degree-rel\tmedian-degree-rel\tsd-degree-rel\n" % source_node)
# 			for source in  total_source_normalized_degree_dict :
# 				nets_number = len(total_source_normalized_degree_dict[source])
# 				min_deg = min(total_source_normalized_degree_dict[source])
# 				max_deg = max(total_source_normalized_degree_dict[source])
# 				mean_deg = statistics.mean(total_source_normalized_degree_dict[source])
# 				median_deg = statistics.median(total_source_normalized_degree_dict[source])
# 				sd_deg = statistics.pstdev(total_source_normalized_degree_dict[source])
# 				file_write.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (source, nets_number, min_deg, max_deg, mean_deg, median_deg, sd_deg))

# 		### compute Z-score for the covid sources
# 		# For each source found in covid network, compute 2 Z-scores.

# 		with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\Further_Analysis\\full_analysis_for_source_to_proteins_undirected\\COVID19_GDDS_Results_%s_to_CovidProteins-withZeros.tsv' % source_node, 'a+') as file_write:
# 			file_write.write("covid %s\t# of nets it is found in\tmin-degree\tmax-degree\tmean-degree\tmedian-degree\tsd-degree\tcovid-degree\tZ-score\n" % source_node)
# 			for source in  covid_sources :
# 				nets_number = len(total_source_degree_dict[source])
# 				min_deg = min(total_source_degree_dict[source])
# 				max_deg = max(total_source_degree_dict[source])
# 				mean_deg = statistics.mean(total_source_degree_dict[source])
# 				median_deg = statistics.median(total_source_degree_dict[source])
# 				sd_deg = statistics.pstdev(total_source_degree_dict[source])
# 				covid_deg = covid_source_degree_dict[source]
# 				if sd_deg != 0 :
# 					z_score = ( covid_deg - mean_deg) / sd_deg
# 				else : 
# 					z_score = 'NaN'
# 				file_write.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (source, nets_number, min_deg, max_deg, mean_deg, median_deg, sd_deg, covid_deg, z_score))


# 		with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\Further_Analysis\\full_analysis_for_source_to_proteins_undirected\\COVID19_GDDS_Results_%s_to_CovidProteins_Relative-withZeros.tsv' % source_node, 'a+') as file_write:
# 			file_write.write("covid %s\t# of nets it is found in\tmin-degree-rel\tmax-degree-rel\tmean-degree-rel\tmedian-degree-rel\tsd-degree-rel\tcovid-degree-rel\tZ-score-rel\n" % source_node)
# 			for source in  covid_sources :
# 				nets_number = len(total_source_normalized_degree_dict[source])
# 				min_deg = min(total_source_normalized_degree_dict[source])
# 				max_deg = max(total_source_normalized_degree_dict[source])
# 				mean_deg = statistics.mean(total_source_normalized_degree_dict[source])
# 				median_deg = statistics.median(total_source_normalized_degree_dict[source])
# 				sd_deg = statistics.pstdev(total_source_normalized_degree_dict[source])
# 				covid_deg = covid_source_norm_degree_dict[source]
# 				if sd_deg != 0 :
# 					z_score = ( covid_deg - mean_deg) / sd_deg
# 				else : 
# 					z_score = 'NaN'
# 				file_write.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (source, nets_number, min_deg, max_deg, mean_deg, median_deg, sd_deg, covid_deg, z_score))










stop = timeit.default_timer()
print(stop - start)  