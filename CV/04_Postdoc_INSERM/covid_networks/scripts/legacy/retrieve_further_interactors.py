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
from collections import Counter
from collections import defaultdict
from collections import OrderedDict

rootdir = 'A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid-2\\SecondVersion\\networks2'
viruses_folders_list = []
for subdir, dirs, files in os.walk(rootdir):
	for name in dirs:
	    viruses_folders_list.append(os.path.join(rootdir, name))
# print(viruses_folders_list)


for path in viruses_folders_list :
	print(path.split("\\")[-1])
	### first get the list of 
	# with open("%s\\direct_interactors.txt" % path, 'a+') as file_write :
	# with open("%s\\direct_and_second_range_interactors.txt" % path, 'a+') as file_write :
	with open("%s\\direct_and_second_and_third_range_interactors.txt" % path, 'a+') as file_write :

		file_write.write("%s\n" % path.split("\\")[-1])
		nodes_path = '%s\\nodes_full.csv' % path
		edges_path = '%s\\edges_full.csv' % path
		# make viral_proteins dict
		descr_dict_nodes, sapiens_dict, virus_dict = {},{},{}
		with open(nodes_path, 'r') as nodecsv: # Open the file    
			nodereader = csv.reader(nodecsv, delimiter=" ") # Read the csv  
			next(nodereader)                   			
			for line in nodereader :
				# print(line)
				node_id = line[0]
				node_name = line[1]
				node_type = line[2].replace('\n', '')
				descr_dict_nodes[node_id] = node_name
				if (node_type == '1') :
					sapiens_dict[node_id] = node_name
				elif (node_type == '0') :
					virus_dict[node_id] = node_name
				else :
					print(line, "Warning : unidentified nodes !\n")
					file_write.write("%s\tWarning : unidentified nodes !\n" % line)
		# print(virus_dict)
		# print(descr_dict_nodes.keys())


		with open(edges_path, 'r') as edgecsv: # Open the file
			edgereader = csv.reader(edgecsv, delimiter=" ") # Read the csv  
			next(edgereader)   
			edges = [tuple(e)[0:2] for e in edgereader] # Retrieve the data
			# print(edges)

		G_interactome = nx.Graph()
		G_interactome.add_nodes_from(descr_dict_nodes.keys())
		G_interactome.add_edges_from(edges)

		print(nx.info(G_interactome))
		file_write.write("%s\n" % nx.info(G_interactome))

		#### fetch viral direct interactors (1st order)
		direct_interactors = []
		for u,v,c in G_interactome.edges(data=True) :
			if (u in virus_dict) :
				direct_interactors.append(v)
			elif (v in virus_dict) :
				direct_interactors.append(u)

		### fetch viral secondary interactors (2nd order)
		secondary_interactors = []
		for u,v,c in G_interactome.edges(data=True) :
			if (u in direct_interactors) :
				secondary_interactors.append(v)
			elif (v in direct_interactors) :
				secondary_interactors.append(u)
		# print(secondary_interactors)

		#### fetch viral tertiary interactors (3rd order)
		tertiary_interactors = []
		for u,v,c in G_interactome.edges(data=True) :
			if (u in secondary_interactors) :
				tertiary_interactors.append(v)
			elif (v in secondary_interactors) :
				tertiary_interactors.append(u)
		tertiary_interactors = set(tertiary_interactors)

		# merge first and second order interactors
		# first_and_second_order_interactors = direct_interactors + secondary_interactors

		# merge first and second and third order interactors
		first_and_second_and_third_order_interactors = list(direct_interactors) + list(secondary_interactors) + list(tertiary_interactors)


		## only first range :
		# for protein in list(set(list(direct_interactors))) :
		# 	print(protein, descr_dict_nodes[protein])
		# 	file_write.write("%s\t%s\n" % (protein, descr_dict_nodes[protein]))

		## for 1st and 2nd range :
		# for protein in list(set(list(first_and_second_order_interactors))) :
		# 	print(protein, descr_dict_nodes[protein])
		# 	file_write.write("%s\t%s\n" % (protein, descr_dict_nodes[protein]))		

		## for 1st and 2nd and 3rd range :
		for protein in first_and_second_and_third_order_interactors :
			print(protein, descr_dict_nodes[protein])
			file_write.write("%s\t%s\n" % (protein, descr_dict_nodes[protein]))	

	



stop = timeit.default_timer()
print(stop - start)  