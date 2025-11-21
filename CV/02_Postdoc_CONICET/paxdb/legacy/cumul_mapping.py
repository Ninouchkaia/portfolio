#!/usr/bin/python
# -*- coding: utf-8 -*-

#script that cumulates mapping table from 1st try and 2nd try.
# first versions were generated with the id mapping tool from uniprot website
# second versions were generated from a mapping given by pax-db website

import timeit
start = timeit.default_timer()
import os
import sys
import re
from commands import getoutput # permet d'obtenir l'output d'une commande bash.

species_list = ["10090-M.musculus", "160490-S.pyogenes", "224308-Spectral_counting_B.subtili", "267671-L.interrogans", "3702-A.thaliana", "449447-M.aeruginosa", "4896-S.pombe", "4932-S.cerevisiae", "511145-E.coli", "593117-Spectral_counting_T.gammatolerans", "6239-C.elegans", "7227-D.melanogaster", "83332-M.tuberculosis", "9031-Spectral_counting_G.gallus", "9606-H.sapiens", "9913-B.taurus", "99287-Spectral_counting_S.typhimurium"]

for species in species_list:
	print species
	add_to_file_dict = {}
	add_to_file = open("/home/nina/scripts/paxdb/abundance_datasets/uniprot_id/%s_dataset.txt_uniprot_id" % species, 'a+')
	add_to_file_list = add_to_file.readlines()
	for line in add_to_file_list :
		line = line.split("\t")
		uniprot_id_of_add_to_file = line[1].replace('\n', '') # remove '\n' only
		add_to_file_dict[line[0]] = uniprot_id_of_add_to_file
	mapping_table = open("/home/nina/scripts/paxdb/abundance_datasets/try_1/mapping_tables/%s_uniprot.txt" % species)
	mapping_table = mapping_table.readlines()
	#print add_to_file_dict
	for i in range(1, len(mapping_table)) :
		mapping_table_splitted = mapping_table[i].split("\t")
		string_id_of_mapping_table = mapping_table_splitted[0]
		uniprot_id_of_mapping_table = mapping_table_splitted[1].replace('\n', '') # remove '\n' only
		if string_id_of_mapping_table in add_to_file_dict :
			if uniprot_id_of_mapping_table != add_to_file_dict[string_id_of_mapping_table] :
				 #print mapping_table[i], "OKOKOKO"
				 add_to_file.write("%s" % mapping_table[i])
			
		else :
			#print mapping_table[i]
			add_to_file.write("%s" % mapping_table[i])
			
	add_to_file.close()

stop = timeit.default_timer()
print stop - start 