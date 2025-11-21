#!/usr/bin/python
# -*- coding: utf-8 -*-

#purpose of the script : re-organize pre-existing tables by putting the last 4 summary lines on top.

#update from version5 (and version6)
#this version is based on version5 rather than 6. It no longer generates a table with frequencies * abundance values.
# this version adds up 4 lines on the top of the tables generated with script version5
# 1. Sum of each amino-acid present in the tables
# 2. Fraction of the total amino-acids represented by each sum of amino-acids
# 3. Sum of each amino-acid ponderated by the corresponding abundance
# 4. Fraction of the total amino-acids represented by each sum of amino-acids ponderated by the corresponding abundance
# it reads the tables generated with script version5 and puts the 4 last lines of each table at their beginning

#update from version4
# we cumulated first and second version of the mapping tables. first versions were generated with the id mapping tool from uniprot website
# second versions were generated from a mapping given by pax-db website

#update from version2
# we avoid to open twice the fasta file and find the matched uniprot id only once
# in order to reach the end of the script more rapidly (originally it took 1350 seconds for the 17 datasets > with this script : 1130 seconds)

#description of the original version
# this script reads a fasta file of proteins and counts the amino acid occurences in each protein.
# this script also reads a dataset of protein abundance data. The abundance datasets are mainly composed of 2 columns : string ID - abundance score
# this script matches the string IDs to uniprot IDs through the use of a mapping file generated with the id mapping tool from uniprot.org.



import timeit
start = timeit.default_timer()
import os
import sys
import re
from Bio import SeqIO # to parse the fasta file
from Bio import SeqUtils # Miscellaneous functions for dealing with sequences.
from Bio.SeqUtils import  ProtParam

# the list of the 17 datasets to be analyzed
species_list = ["10090-M.musculus", "160490-S.pyogenes", "267671-L.interrogans", "3702-A.thaliana", "449447-M.aeruginosa", "4896-S.pombe", "4932-S.cerevisiae", "511145-E.coli", "6239-C.elegans", "7227-D.melanogaster", "83332-M.tuberculosis", "9606-H.sapiens", "9913-B.taurus"]
#species_list = ["224308-Spectral_counting_B.subtili"]
#species_list = ["224308-Spectral_counting_B.subtili", "593117-Spectral_counting_T.gammatolerans", "9031-Spectral_counting_G.gallus", "99287-Spectral_counting_S.typhimurium"]


for species in species_list:
	print species

	## we start by creating a list of uniprot ids that were mapped to string ids, from text files located in the folder mapping_tables 
	table = open("/home/nina/scripts/paxdb/abundance_datasets/tables/%s_table.txt" % species, "rU")
	table = table.readlines()
	file_write = open("/home/nina/scripts/paxdb/abundance_datasets/tables/tables_2/%s_table2.txt" % species,'a+')
	
	first_line = table[0]
	file_write.write(first_line)
	
	four_last_lines = table[-4:]
	for line in four_last_lines :
		file_write.write(line)
	
	for line in table[1:-4] :
		file_write.write(line)
	
	 



	


stop = timeit.default_timer()
print stop - start 