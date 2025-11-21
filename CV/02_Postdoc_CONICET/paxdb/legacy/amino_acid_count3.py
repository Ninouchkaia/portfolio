#!/usr/bin/python
# -*- coding: utf-8 -*-

#update from version2
# we avoid to open twice the fasta file and find the matched uniprot id only once
# in order to reach the end of the script more rapidly (originally it took 1350 seconds for the 17 datasets > with this script : 1130 seconds)

#description of the original version
# this script reads a fasta file of proteins and counts the amino acid occurences in each protein.
# this script also reads a dataset of protein abundance data. The abundance datasets are mainly composed of 2 columns : string ID - abundance score
# this script matches the string IDs to uniprot IDs through the use of a mapping file generated with the id mapping tool from uniprot.org.
# it creates different files that serve to make a final file containing the uniprot id names of the proteins, their abundance (cf. pax-db databse), and the counts of amino acids for each protein.
# it has to be executed where the final files will be saved
# 


import timeit
start = timeit.default_timer()
import os
import sys
import re
from Bio import SeqIO # to parse the fasta file
from Bio import SeqUtils # Miscellaneous functions for dealing with sequences.
from Bio.SeqUtils import  ProtParam

# the list of the 17 datasets to be analyzed
#species_list = ["10090-M.musculus", "160490-S.pyogenes", "224308-Spectral_counting_B.subtili", "267671-L.interrogans", "3702-A.thaliana", "449447-M.aeruginosa", "4896-S.pombe", "4932-S.cerevisiae", "511145-E.coli", "593117-Spectral_counting_T.gammatolerans", "6239-C.elegans", "7227-D.melanogaster", "83332-M.tuberculosis", "9031-Spectral_counting_G.gallus", "9606-H.sapiens", "9913-B.taurus", "99287-Spectral_counting_S.typhimurium"]
species_list = ["9606-H.sapiens"]


for species in species_list:
	print species

	## we start by creating a list of uniprot ids that were mapped to string ids, from text files located in the folder mapping_tables 
	handle_mapping_table = open("/home/nina/scripts/paxdb/abundance_datasets/mapping_tables/%s_uniprot.txt" % species, "rU")
	read_mapping_table = handle_mapping_table.readlines()
	handle_fasta = open("/home/nina/scripts/paxdb/abundance_datasets/fasta/%s.fasta" % species, "rU") #opens the fasta file containing the sequences and their information
	handle_dataset = open("/home/nina/scripts/paxdb/abundance_datasets/datasets/%s_dataset.txt" % species, "rU") #opens the abundance dataset
	read_dataset = handle_dataset.readlines()
	abundance_list = []
	file_write5 = open("%s_5.txt" % species,'a+') 
	#print the keys of the dict (the amino acids) in the alphabetical order
	file_write5.write("\t") #leaves the first case empty for the column of protein names to be written on next lines (first column)
	file_write5.write("A\t")
	file_write5.write("C\t")
	file_write5.write("D\t")
	file_write5.write("E\t")
	file_write5.write("F\t")
	file_write5.write("G\t")
	file_write5.write("H\t")
	file_write5.write("I\t")
	file_write5.write("K\t")
	file_write5.write("L\t")
	file_write5.write("M\t")
	file_write5.write("N\t")
	file_write5.write("P\t")
	file_write5.write("Q\t")
	file_write5.write("R\t")
	file_write5.write("S\t")
	file_write5.write("T\t")
	file_write5.write("V\t")
	file_write5.write("W\t")
	file_write5.write("Y\t")
	file_write5.write("abundance")
	file_write5.write("\n")
	for record in SeqIO.parse(handle_fasta, "fasta") : # for each uniprot in the fasta file : (for another obscur reason, the 2 for loops have to go in this order)
		for mapping in read_mapping_table :
			mapping = mapping.replace('\n', '') # remove '\n' only
			mapping = mapping.split("\t")
			if record.id[3:9] == mapping[1][0:6] : # for a given uniprot id : if the uniprot id from the fasta file and the mapping file is identical
				for abundance_data in read_dataset :
					#abundance_data = abundance_data.replace('\n', '') # remove '\n' only
					abundance_data = abundance_data.split("\t")
					#print abundance_data
					if mapping[0] == abundance_data[1]: # for the same given uniprot id, corresponds a string id : if the string id from the mapping file and the abundance file is identical
						abundance_list.append("%s\t%s\t%s" % (mapping[1], abundance_data[1], abundance_data[2])) # we retrieve the uniprot id of the mapping file, the string id from the abundance data file and the abundance data for the given uniprot/string id
		
						my_seq = record.seq #get each sequence as a string
						X = ProtParam.ProteinAnalysis("%s" % my_seq)  # creates the variable containing the sequence, to be passed in the count_amino_acids function.
						dict_seq = X.count_amino_acids() # will count frequency of amino acid in each sequence and return as a dict
						

						
						#print the list of proteins in a column
						file_write5.write("%s\t" % record.id ) 
						
						#print the counts of each amino acid of each protein as a table
						file_write5.write("%s\t" % dict_seq["A"]) 
						file_write5.write("%s\t" % dict_seq["C"])
						file_write5.write("%s\t" % dict_seq["D"])
						file_write5.write("%s\t" % dict_seq["E"])
						file_write5.write("%s\t" % dict_seq["F"])
						file_write5.write("%s\t" % dict_seq["G"])
						file_write5.write("%s\t" % dict_seq["H"])
						file_write5.write("%s\t" % dict_seq["I"])
						file_write5.write("%s\t" % dict_seq["K"])
						file_write5.write("%s\t" % dict_seq["L"])
						file_write5.write("%s\t" % dict_seq["M"])
						file_write5.write("%s\t" % dict_seq["N"])
						file_write5.write("%s\t" % dict_seq["P"])
						file_write5.write("%s\t" % dict_seq["Q"])
						file_write5.write("%s\t" % dict_seq["R"])
						file_write5.write("%s\t" % dict_seq["S"])
						file_write5.write("%s\t" % dict_seq["T"])
						file_write5.write("%s\t" % dict_seq["V"])
						file_write5.write("%s\t" % dict_seq["W"])
						file_write5.write("%s\t" % dict_seq["Y"])
			
						#print the abundance values for each protein as a column 
						file_write5.write("%s" % (abundance_data[2]))
						#file_write5.write("%s\n" % (abundance_data[2])) # for some datasets we have to jump a line otherwise things get written one a single line

	handle_mapping_table.close()				  
	handle_fasta.close()
	handle_dataset.close()


stop = timeit.default_timer()
print stop - start 