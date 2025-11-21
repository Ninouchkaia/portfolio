#!/usr/bin/python
# -*- coding: utf-8 -*-

#script that maps all string ids from abundance datasets through reading a mapping table given by pax-db website.

import timeit
start = timeit.default_timer()
import os
import sys
import re
from commands import getoutput # permet d'obtenir l'output d'une commande bash.

list_of_datasets = ["10090-M.musculus_dataset.txt", "160490-S.pyogenes_dataset.txt", "224308-Spectral_counting_B.subtili_dataset.txt", "267671-L.interrogans_dataset.txt", "3702-A.thaliana_dataset.txt", "449447-M.aeruginosa_dataset.txt", "4896-S.pombe_dataset.txt", "4932-S.cerevisiae_dataset.txt", "511145-E.coli_dataset.txt", "593117-Spectral_counting_T.gammatolerans_dataset.txt", "6239-C.elegans_dataset.txt", "7227-D.melanogaster_dataset.txt", "83332-M.tuberculosis_dataset.txt", "9031-Spectral_counting_G.gallus_dataset.txt", "9606-H.sapiens_dataset.txt", "9913-B.taurus_dataset.txt", "99287-Spectral_counting_S.typhimurium_dataset.txt"]
#list_of_datasets = ["449447-M.aeruginosa_dataset.txt"]

for abundance_dataset in list_of_datasets :
	print abundance_dataset

	dataset = open("/home/nina/scripts/paxdb/abundance_datasets/datasets/%s" % abundance_dataset)
	dataset = dataset.readlines()

	mapping_table = open("/home/nina/scripts/paxdb/abundance_datasets/paxdb-uniprot-links-v3_0_0.txt")
	mapping_table = mapping_table.readlines()
	
	mapped_file_write = open("/home/nina/scripts/paxdb/abundance_datasets/uniprot_id/%s_uniprot_id" % abundance_dataset, 'a+')
	
	unmapped_file_write = open("/home/nina/scripts/paxdb/abundance_datasets/uniprot_id_not/%s_uniprot_id_not" % abundance_dataset, 'a+')
	
	fasta_retrieve_write = open("/home/nina/scripts/paxdb/abundance_datasets/fasta_retrieve/%s_fasta_retrieve" % abundance_dataset, 'a+')
	
	mapped_dict = {}
	for mapping in mapping_table :
		mapping = mapping.replace('\n', '') # remove '\n' only
		mapping = mapping.split("\t")
		string_id = mapping[0]
		uniprot_id = mapping[1]
		mapped_dict[mapping[0]] = mapping[1]
	
	for line in dataset :
		line = line.replace('\n', '') # remove '\n' only
		line = line.split("\t")
		string_id_to_map = line[1]

		if string_id_to_map in mapped_dict :
			if (len(mapped_dict[string_id_to_map]) > 7) :
				mapped_file_write.write("%s\t%s\n" % (string_id_to_map, mapped_dict[string_id_to_map][0:6]))
				mapped_file_write.write("%s\t%s\n" % (string_id_to_map, mapped_dict[string_id_to_map][7:13]))
				fasta_retrieve_write.write("%s\n" % (mapped_dict[string_id_to_map][0:6]))
				fasta_retrieve_write.write("%s\n" % (mapped_dict[string_id_to_map][7:13]))
			else:
				mapped_file_write.write("%s\t%s\n" % (string_id_to_map, mapped_dict[string_id_to_map]))
				fasta_retrieve_write.write("%s\n" % (mapped_dict[string_id_to_map]))
		else :
			unmapped_file_write.write("%s\n" % string_id_to_map)


	mapped_file_write.close()
	unmapped_file_write.close()

stop = timeit.default_timer()
print stop - start 