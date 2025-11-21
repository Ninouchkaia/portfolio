#!/usr/bin/python
# -*- coding: utf-8 -*-

import timeit
start = timeit.default_timer()
import os
import sys
import re
from Bio import SeqIO # to parse the fasta file
from Bio import SeqUtils # Miscellaneous functions for dealing with sequences.
from Bio.SeqUtils import  ProtParam


species_list = ["10090-M.musculus", "160490-S.pyogenes", "224308-Spectral_counting_B.subtili", "267671-L.interrogans", "3702-A.thaliana", "449447-M.aeruginosa", "4896-S.pombe", "4932-S.cerevisiae", "511145-E.coli", "593117-Spectral_counting_T.gammatolerans", "6239-C.elegans", "7227-D.melanogaster", "83332-M.tuberculosis", "9031-Spectral_counting_G.gallus", "9606-H.sapiens", "9913-B.taurus", "99287-Spectral_counting_S.typhimurium"]

handle_fasta = open("/home/nina/scripts/paxdb/abundance_datasets/fasta_files/processed_paxdb-protein-sequences.fa") #opens the fasta file containing the sequences and their information
handle_fasta_list = handle_fasta.readlines()

for species in species_list :
	print species
	file_write = open("/home/nina/scripts/paxdb/abundance_datasets/fasta_files/%s_processed_fasta.fa" % species, 'a+')
	file_open = open("/home/nina/scripts/paxdb/abundance_datasets/datasets/%s_dataset.txt" % species, "rU") #opens the abundance dataset
	list_of_id_in_dataset = file_open.readlines()
	list_of_identifiant_prot = []
	for line in list_of_id_in_dataset :
		line = line.split("\t")
		string_id = line[1]
		list_of_identifiant_prot.append(string_id)

	for identifiant in list_of_identifiant_prot :
		for i in range(0, len(handle_fasta_list)) : # for each uniprot in the fasta file : (for another obscur reason, the 2 for loops have to go in this order) 
			if identifiant in handle_fasta_list[i]:
				file_write.write("%s%s" % (handle_fasta_list[i], handle_fasta_list[i+1]))

			

handle_fasta.close()

stop = timeit.default_timer()
print stop - start 