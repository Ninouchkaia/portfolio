#!/usr/bin/python
# -*- coding: utf-8 -*-
# jupyter notebook --browser='C:/Program Files (x86)/Google/Chrome/Application/chrome.exe %s'


import timeit
start = timeit.default_timer()


import gseapy
import os
import sys
import pandas as pd

# assign a list object to enrichr
# my_gene_list = ['SCARA3', 'LOC100044683', 'CMBL', 'CLIC6', 'IL13RA1', 'TACSTD2', 'DKKL1', 'CSF1',
#      'SYNPO2L', 'TINAGL1', 'PTX3', 'BGN', 'HERC1', 'EFNA1', 'CIB2', 'PMP22', 'TMEM173']

# gseapy.enrichr(gene_list=l, description='pathway', gene_sets='KEGG_2016', outfile='test')


rootdir = 'A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid-2\\Virus_host_interactomes_thresh25\\thresh=0.25'

# for subdir, dirs, files in os.walk(rootdir):
#     for file in files:
#         print(os.path.join(subdir, file))

viruses_folders_list = []
for subdir, dirs, files in os.walk(rootdir):
	for name in dirs:
	    viruses_folders_list.append(os.path.join(rootdir, name))
# print(viruses_folders_list)

full_gene_list_dict = {}
full_interactors_list = []
for path in viruses_folders_list :
	virus_name = path.split("\\")[-1]
	my_gene_list = []
	myfile = "%s\\only_Fourth_range_interactors.txt" % path
	with open(myfile, 'r') as file_read :
		data = file_read.readlines()
		for line in data[6:] :
			line = line.split('\t')
			gene = line[1].replace('\n','')
			my_gene_list.append(gene)
			full_interactors_list.append(gene)
			full_gene_list_dict[virus_name] = my_gene_list
full_interactors_list = list(set(full_interactors_list))



df = pd.DataFrame(columns = full_gene_list_dict.keys(), index = full_interactors_list)

for virus in full_gene_list_dict :
	for gene in full_gene_list_dict[virus] :
		df.loc[gene,virus] = 1
# print(df.loc['TP53'])


with open("gene_virus_table_OnlyFourthRange.txt", 'a+') as file_write :
	for virus in full_gene_list_dict :
		for gene in full_gene_list_dict[virus] :
			file_write.write("%s\t%s\n" % (virus, gene))

# with open("gene_virus_table_without_covid.txt", 'a+') as file_write :
# 	for virus in full_gene_list_dict :
# 		if (virus != "covid19") :
# 			for gene in full_gene_list_dict[virus] :
# 				file_write.write("%s\t%s\n" % (virus, gene))


stop = timeit.default_timer()
print(stop - start)  

