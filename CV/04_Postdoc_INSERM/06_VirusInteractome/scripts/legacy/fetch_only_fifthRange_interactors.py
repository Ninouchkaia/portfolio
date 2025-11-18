#!/usr/bin/python
# -*- coding: utf-8 -*-
# jupyter notebook --browser='C:/Program Files (x86)/Google/Chrome/Application/chrome.exe %s'


import timeit
start = timeit.default_timer()


import gseapy
import os
import sys

rootdir = 'A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid-2\\Virus_host_interactomes_thresh25\\thresh=0.25'


viruses_folders_list = []
for subdir, dirs, files in os.walk(rootdir):
	for name in dirs:
	    viruses_folders_list.append(os.path.join(rootdir, name))
print(viruses_folders_list)

full_gene_list_dict = {}
for path in viruses_folders_list :
	virus_name = path.split("\\")[-1]
	print(virus_name)
	my_gene_list = []
	myfile = "%s\\only_Fifth_range_interactors.txt" % path
	print(myfile)
	with open(myfile, 'r') as file_read :
		data = file_read.readlines()
		for line in data[6:] :
			line = line.split('\t')
			gene = line[1].replace('\n','')
			my_gene_list.append(gene)
			full_gene_list_dict[virus_name] = my_gene_list

with open("only_Fifth_order_genes_list.txt", "a+") as file_write :
	for virus in full_gene_list_dict :
		file_write.write("%s\t" % virus)
		for gene in full_gene_list_dict[virus][0:-1] :
			file_write.write("%s\t" % gene)
		file_write.write("%s\n" % full_gene_list_dict[virus][-1])

stop = timeit.default_timer()
print(stop - start)  
