#!/usr/bin/python
# -*- coding: utf-8 -*-
# jupyter notebook --browser='C:/Program Files (x86)/Google/Chrome/Application/chrome.exe %s'


import timeit
start = timeit.default_timer()


import gseapy
import os
import sys

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

for path in viruses_folders_list :
	my_gene_list = []
	if (path.split("\\")[-1] == 'covid19') :
		myfile = "%s\\direct_interactors.txt" % path
		with open(myfile, 'r') as file_read :
			data = file_read.readlines()
			for line in data[6:] :
				line = line.split('\t')
				gene = line[1].replace('\n','')
				my_gene_list.append(gene)
				
				enr = gseapy.enrichr(gene_list=my_gene_list,
	                 gene_sets=['GO_Biological_Process_2013','GO_Cellular_Component_2013','GO_Molecular_Function_2013',],
	                 # gene_sets=['Reactome_2016',],
	                 organism='Human', # don't forget to set organism to the one you desired! e.g. Yeast
	                 description='test_name',
	                 outdir='%s/enrichr_Reactome_2016' % path,
	                 # no_plot=True,
	                 cutoff=1 # test dataset, use lower value from range(0,1)
	                )


stop = timeit.default_timer()
print(stop - start)  