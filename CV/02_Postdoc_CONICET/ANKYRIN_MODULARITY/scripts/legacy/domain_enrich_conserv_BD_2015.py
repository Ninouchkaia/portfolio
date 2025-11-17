#!/usr/bin/python
# -*- coding: utf-8 -*-

import timeit
start = timeit.default_timer()

import sys
import os
from collections import Counter
import re
from Bio import SeqIO


my_fasta = sys.argv[1]

num_homologs = int(sys.argv[2])

conserv_dict = {}
with open("domain_conservation_percentages_in_BD_%s.txt" % (num_homologs), 'rU') as file_open :
	data = file_open.readlines()
	for line in data :
		line = line.split()
		domain_name = line[0] 
		domain_conservation = float(line[1].replace("\n",""))
		conserv_dict[domain_name] = domain_conservation

average_domain_conservation = sum(conserv_dict.values()) / len(conserv_dict)
print average_domain_conservation


with open("Pfam_domains_in_BD_2038_Zscores.txt", 'rU') as file_open :
	data = file_open.readlines()
	with open("Pfam_domains_in_%s_MaxHomologs_%s_Zscores_color.txt" % (my_fasta,num_homologs), 'a+') as file_write :
		file_write.write("domain_name\tdomain_enrichment_log(obs/exp)\tZscores_enrichment\tdomain_conservation_over_%sHomologs\tcolor\n" % num_homologs)
		for line in data[1:] :
			line = line.split("\t")
			domain_name = line[0]
			domain_enrichment = float(line[4].replace(",", "."))
			domain_enrichment_Zscore = float((line[5].replace(",", ".")).replace("\n",""))

			
			if domain_name in conserv_dict :
				domain_conservation = conserv_dict[domain_name]
				if domain_enrichment_Zscore > 1 and domain_conservation > average_domain_conservation:
					color = "\"red\""
				else :
					color = "\"black\""
				file_write.write("%s\t%f\t%f\t%f\t%s\n" % (domain_name, domain_enrichment, domain_enrichment_Zscore, domain_conservation, color))
			else :
				print "domains %s" % domain_name, "not found in dict"




stop = timeit.default_timer()
print stop - start

