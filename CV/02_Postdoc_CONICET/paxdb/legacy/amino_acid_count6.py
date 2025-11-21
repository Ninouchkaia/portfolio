#!/usr/bin/python
# -*- coding: utf-8 -*-


#update from version5
# this version of the script counts amino acids multiplied by abundance values of the corresponding protein

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
	handle_mapping_table = open("/home/nina/scripts/paxdb/abundance_datasets/uniprot_id/%s_dataset.txt_uniprot_id" % species, "rU")
	read_mapping_table = handle_mapping_table.readlines()
	handle_fasta = open("/home/nina/scripts/paxdb/abundance_datasets/fasta_files/%s.fasta" % species, "rU") #opens the fasta file containing the sequences and their information
	handle_dataset = open("/home/nina/scripts/paxdb/abundance_datasets/datasets/%s_dataset.txt" % species, "rU") #opens the abundance dataset
	read_dataset = handle_dataset.readlines()
	abundance_dict = {}
	file_write = open("/home/nina/scripts/paxdb/abundance_datasets/tables/processed_tables/%s_processed_table.txt" % species,'a+') 
	#print the keys of the dict (the amino acids) in the alphabetical order
	file_write.write("uniprot_id\t") #leaves the first case empty for the column of protein names to be written on next lines (first column)
	file_write.write("string_id\t") #leaves the second case empty for the column of protein names in string id to be written on next lines (second column)
	file_write.write("A\t")
	file_write.write("C\t")
	file_write.write("D\t")
	file_write.write("E\t")
	file_write.write("F\t")
	file_write.write("G\t")
	file_write.write("H\t")
	file_write.write("I\t")
	file_write.write("K\t")
	file_write.write("L\t")
	file_write.write("M\t")
	file_write.write("N\t")
	file_write.write("P\t")
	file_write.write("Q\t")
	file_write.write("R\t")
	file_write.write("S\t")
	file_write.write("T\t")
	file_write.write("V\t")
	file_write.write("W\t")
	file_write.write("Y\t")
	file_write.write("abundance")
	file_write.write("\n")
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
						#abundance_dict[abundance_data[1]] = abundance_data[2] #dictionnary of abundance values
		
						my_seq = record.seq #get each sequence as a string
						X = ProtParam.ProteinAnalysis("%s" % my_seq)  # creates the variable containing the sequence, to be passed in the count_amino_acids function.
						dict_seq = X.count_amino_acids() # will count frequency of amino acid in each sequence and return as a dict
						

						
						#print the list of proteins in the first column
						file_write.write("%s\t" % record.id ) 
						
						#print the list of proteins as string id in the second column
						file_write.write("%s\t" % abundance_data[1] ) 
						
						#print the counts of each amino acid of each protein as a table
						file_write.write("%s\t" % (float(dict_seq["A"])*(float(abundance_data[2])))) 
						file_write.write("%s\t" % (float(dict_seq["C"])*(float(abundance_data[2]))))
						file_write.write("%s\t" % (float(dict_seq["D"])*(float(abundance_data[2]))))
						file_write.write("%s\t" % (float(dict_seq["E"])*(float(abundance_data[2]))))
						file_write.write("%s\t" % (float(dict_seq["F"])*(float(abundance_data[2]))))
						file_write.write("%s\t" % (float(dict_seq["G"])*(float(abundance_data[2]))))
						file_write.write("%s\t" % (float(dict_seq["H"])*(float(abundance_data[2]))))
						file_write.write("%s\t" % (float(dict_seq["I"])*(float(abundance_data[2]))))
						file_write.write("%s\t" % (float(dict_seq["K"])*(float(abundance_data[2]))))
						file_write.write("%s\t" % (float(dict_seq["L"])*(float(abundance_data[2]))))
						file_write.write("%s\t" % (float(dict_seq["M"])*(float(abundance_data[2]))))
						file_write.write("%s\t" % (float(dict_seq["N"])*(float(abundance_data[2]))))
						file_write.write("%s\t" % (float(dict_seq["P"])*(float(abundance_data[2]))))
						file_write.write("%s\t" % (float(dict_seq["Q"])*(float(abundance_data[2]))))
						file_write.write("%s\t" % (float(dict_seq["R"])*(float(abundance_data[2]))))
						file_write.write("%s\t" % (float(dict_seq["S"])*(float(abundance_data[2]))))
						file_write.write("%s\t" % (float(dict_seq["T"])*(float(abundance_data[2]))))
						file_write.write("%s\t" % (float(dict_seq["V"])*(float(abundance_data[2]))))
						file_write.write("%s\t" % (float(dict_seq["W"])*(float(abundance_data[2]))))
						file_write.write("%s\t" % (float(dict_seq["Y"])*(float(abundance_data[2]))))
			
						#print the abundance values for each protein as a column 
						file_write.write("%s" % (abundance_data[2]))
						#file_write.write("%s\n" % (abundance_data[2])) # for some datasets we have to jump a line otherwise things get written one a single line
					#else :
						#print mapping[0]

	handle_mapping_table.close()				  
	handle_fasta.close()
	handle_dataset.close()
	


stop = timeit.default_timer()
print stop - start 