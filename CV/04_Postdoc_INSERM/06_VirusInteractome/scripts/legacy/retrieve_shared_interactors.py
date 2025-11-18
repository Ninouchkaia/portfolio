#!/usr/bin/python
# -*- coding: utf-8 -*-
# jupyter notebook --browser='C:/Program Files (x86)/Google/Chrome/Application/chrome.exe %s'


import timeit
start = timeit.default_timer()

import sys
import os
import csv

rootdir = 'A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid-2\\SecondVersion\\networks2'
viruses_folders_list = []
for subdir, dirs, files in os.walk(rootdir):
	for name in dirs:
	    viruses_folders_list.append(os.path.join(rootdir, name))
# print(viruses_folders_list)


interactors_dict = {}
full_interactors_list = []
for path in viruses_folders_list :
	virus_name = path.split("\\")[-1]
	print(virus_name)

	my_protein_list = []
	myfile = "%s\\nodes.csv" % path
	with open(myfile, 'r') as file_read :
		data = file_read.readlines()
		for line in data[1:] :
			line = line.split(',')
			protein = line[1]
			organism = line[2].replace('\n','')
			if organism == 'Homo sapiens' :
				my_protein_list.append(protein)
				full_interactors_list.append(protein)
				interactors_dict[virus_name] = my_protein_list
full_interactors_list = list(set(full_interactors_list))
print(len(full_interactors_list))

covid_prots_direct = interactors_dict["Human_SARS_coronavirus_2"]

# covid_prots = interactors_dict["Human_SARS_coronavirus_2"]
# shared_prot = {}
# virus_sharing_prots = {}

# for virus in interactors_dict :
# 	# print(virus, interactors_dict[virus])
# 	if virus != "Human_SARS_coronavirus_2" :
# 		for prot in interactors_dict[virus] :
# 			if prot in covid_prots :
# 				if virus not in virus_sharing_prots :
# 					virus_sharing_prots[virus] = [prot]
# 				else : 
# 					virus_sharing_prots[virus].append(prot)
# 				# print(virus, prot)
# 				if prot not in shared_prot:
# 					shared_prot[prot] = [virus]
# 				else :
# 					shared_prot[prot].append(virus)


# # for prot in shared_prot :
# # 	print(prot, shared_prot[prot])

# for virus in virus_sharing_prots :
# 	print(virus)

# # with open("virus_gene_table.txt", "a+") as file_write :
# # 	for virus in interactors_dict :
# # 		for gene in interactors_dict[virus]:
# # 			file_write.write("%s\t%s\n" % (virus,gene))














interactors_dict2 = {}

with open("virus_gene_table_second_order.txt", 'r') as file_read :
	data = file_read.readlines()
	for line in data :
		line = line.split('\t')
		virus = line[0]
		prot = line[1].replace('\n','')
		if virus not in interactors_dict2 :
			interactors_dict2[virus] = [prot]
		else : 
			interactors_dict2[virus].append(prot)


covid_prots = interactors_dict2["Human_SARS_coronavirus_2"]
shared_prot = {}
virus_sharing_prots = {}
for virus in interactors_dict2 :
	# print(virus, interactors_dict2[virus])
	if virus != "Human_SARS_coronavirus_2" :
		for prot in interactors_dict2[virus] :
			if prot in covid_prots_direct :
				if virus not in virus_sharing_prots :
					virus_sharing_prots[virus] = [prot]
				else : 
					virus_sharing_prots[virus].append(prot)
				# print(virus, prot)
				if prot not in shared_prot:
					shared_prot[prot] = [virus]
				else :
					shared_prot[prot].append(virus)

for virus in virus_sharing_prots :
	print(virus, len(virus_sharing_prots[virus]))


# for prot in shared_prot :
# 	print(prot, shared_prot[prot])


stop = timeit.default_timer()
print(stop - start)  





