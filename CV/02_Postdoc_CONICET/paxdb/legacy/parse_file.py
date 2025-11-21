#!/usr/bin/python
# -*- coding: utf-8 -*-

import timeit
start = timeit.default_timer()
import os
import sys
import re
# from Bio import SeqIO # to parse the fasta file
# from Bio import SeqUtils # Miscellaneous functions for dealing with sequences.
# from Bio.SeqUtils import  ProtParam
# import numpy 

file_open = open("GtRNAdb-all-tRNAs.fa", 'rU')
trna_genes_list = file_open.readlines()

file_write = open("trna_from_gtrnadb.txt", 'a+')
file_write.write("species\ttRNA name\tamino-acid\tcodon\n")

for line in trna_genes_list :
	if line[0] == ">" :
		species = (line.split(".trna"))[0]
		species = re.sub(r"_chr.*", "", species)
		species = species.replace(">", "")
		trna_name = (((line.split("trna"))[1]).split(" "))[0]
		

		if "SeC(" in line :
			aa = (((line.split("  "))[1]).split(" "))[0]
			codon = (((line.split("("))[4]).split(")"))[0]
		else :
			aa = (((line.split(")"))[1]).split("("))[0]
			aa = aa.replace(" ", "")
			codon = (((line.split("("))[2]).split(")"))[0]
	file_write.write("%s\t%s\t%s\t%s\n" % (species, trna_name, aa, codon))
	




stop = timeit.default_timer()
print stop - start 