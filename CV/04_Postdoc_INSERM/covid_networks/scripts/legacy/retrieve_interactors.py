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
print(viruses_folders_list)


# interactors_dict = {}
# full_interactors_list = []
# for path in viruses_folders_list :
# 	virus_name = path.split("\\")[-1]
# 	print(virus_name)

# 	my_protein_list = []
# 	myfile = "%s\\nodes.csv" % path
# 	with open(myfile, 'r') as file_read :
# 		data = file_read.readlines()
# 		for line in data[1:] :
# 			line = line.split(',')
# 			protein = line[1]
# 			organism = line[2].replace('\n','')
# 			if organism == 'Homo sapiens' :
# 				my_protein_list.append(protein)
# 				full_interactors_list.append(protein)
# 				interactors_dict[virus_name] = my_protein_list
# full_interactors_list = list(set(full_interactors_list))
# print(len(full_interactors_list))
# for virus in interactors_dict :
# 	print(virus, interactors_dict[virus])


# with open("genes_list.txt", "a+") as file_write :
# 	for virus in interactors_dict :
# 		file_write.write("%s\t" % virus)
# 		for gene in interactors_dict[virus][0:-1] :
# 			file_write.write("%s\t" % gene)
# 		file_write.write("%s\n" % interactors_dict[virus][-1])

# with open("virus_gene_table.txt", "a+") as file_write :
# 	for virus in interactors_dict :
# 		for gene in interactors_dict[virus]:
# 			file_write.write("%s\t%s\n" % (virus,gene))



# full_gene_list_dict = {}
# for path in viruses_folders_list :
# 	virus_name = path.split("\\")[-1]
# 	print(virus_name)
# 	my_gene_list = []
# 	myfile = "%s\\direct_and_second_range_interactors.txt" % path
# 	with open(myfile, 'r') as file_read :
# 		data = file_read.readlines()
# 		for line in data[6:] :
# 			line = line.split('\t')
# 			gene = line[1].replace('\n','')
# 			my_gene_list.append(gene)
# 			full_gene_list_dict[virus_name] = my_gene_list

# with open("first_and_second_order_genes_list.txt", "a+") as file_write :
# 	for virus in full_gene_list_dict :
# 		file_write.write("%s\t" % virus)
# 		for gene in full_gene_list_dict[virus][0:-1] :
# 			file_write.write("%s\t" % gene)
# 		file_write.write("%s\n" % full_gene_list_dict[virus][-1])

# with open("virus_gene_table_second_order.txt", "a+") as file_write :
# 	for virus in full_gene_list_dict :
# 		for gene in full_gene_list_dict[virus]:
# 			file_write.write("%s\t%s\n" % (virus,gene))



full_gene_list_dict = {}
for path in viruses_folders_list :
	virus_name = path.split("\\")[-1]
	print(virus_name)
	my_gene_list = []
	myfile = "%s\\direct_and_second_and_third_range_interactors.txt" % path
	with open(myfile, 'r') as file_read :
		data = file_read.readlines()
		for line in data[6:] :
			line = line.split('\t')
			gene = line[1].replace('\n','')
			my_gene_list.append(gene)
			full_gene_list_dict[virus_name] = my_gene_list

with open("first_and_second_and_third_order_genes_list.txt", "a+") as file_write :
	for virus in full_gene_list_dict :
		file_write.write("%s\t" % virus)
		for gene in full_gene_list_dict[virus][0:-1] :
			file_write.write("%s\t" % gene)
		file_write.write("%s\n" % full_gene_list_dict[virus][-1])

with open("virus_gene_table_third_order.txt", "a+") as file_write :
	for virus in full_gene_list_dict :
		for gene in full_gene_list_dict[virus]:
			file_write.write("%s\t%s\n" % (virus,gene))


stop = timeit.default_timer()
print(stop - start)  