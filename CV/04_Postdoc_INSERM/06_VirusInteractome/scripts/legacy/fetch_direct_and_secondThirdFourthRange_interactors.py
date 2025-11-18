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
print(viruses_folders_list)

full_gene_list_dict = {}
for path in viruses_folders_list :
	virus_name = path.split("\\")[-1]
	print(virus_name)
	my_gene_list = []
	myfile = "%s\\direct_and_secondThirdFourth_range_interactors.txt" % path
	print(myfile)
	with open(myfile, 'r') as file_read :
		data = file_read.readlines()
		for line in data[6:] :
			line = line.split('\t')
			gene = line[1].replace('\n','')
			my_gene_list.append(gene)
			full_gene_list_dict[virus_name] = my_gene_list

with open("first_and_secondThirdFourth_order_genes_list.txt", "a+") as file_write :
	for virus in full_gene_list_dict :
		file_write.write("%s\t" % virus)
		for gene in full_gene_list_dict[virus][0:-1] :
			file_write.write("%s\t" % gene)
		file_write.write("%s\n" % full_gene_list_dict[virus][-1])

stop = timeit.default_timer()
print(stop - start)  


# viruses <- gcSample[,1]
# Freq <- count.fields("genes_list.txt", sep = '\t')
# myhist <-list(breaks=viruses, counts=Freq, density=Freq/diff(viruses),xname="viruses")
# class(myhist) <- "histogram"
# plot(myhist)

# # barchart with added parameters
# viruses <- c(gcSample[,1])
# viral_interactors <- c(count.fields("genes_list.txt", sep = '\t'))
# barplot(viral_interactors)

# barplot(viral_interactors,
# main = "Number of direct human protein interactors per virus",
# xlab = "Virus name",
# ylab = "Number of proteins",
# names.arg = viruses,
# col = "darkred",
# horiz = FALSE)


# count.fields(dat, sep = ',')
# # [1] 2 3 2 2 2 2 3 3 7
# max(count.fields(dat, sep = ','))
# # [1] 7


# column_number <- max(count.fields("genes_list.txt", sep = '\t'))

# read.table(dat, header = FALSE, sep = ",", col.names = paste0("V",seq_len(7)), fill = TRUE)

# gcSampleNina <- read.table("genes_list.txt", header = FALSE, sep = "\t", col.names = paste0("V",seq_len(max(count.fields("genes_list.txt", sep = '\t')))),, dec = ".", fill = TRUE)

# dfNina <- data.frame(t(gcSampleNina[-1]))
# colnames(dfNina) <- gcSampleNina[, 1]