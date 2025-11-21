#!/usr/bin/python
# -*- coding: utf-8 -*-

#script to use to calculate the frequencies of amino acids in a fasta file


import timeit
start = timeit.default_timer()
import os
import sys
import re
from Bio import SeqIO # to parse the fasta file
from Bio import SeqUtils # Miscellaneous functions for dealing with sequences.
from Bio.SeqUtils import  ProtParam

file_write5 = open("/home/nina/scripts/elm_in_ank/fasta/amino_acid_count/ank_reviewed_canon.txt", 'a+')

handle_fasta = open("/home/nina/scripts/elm_in_ank/fasta/ank_reviewed_canon.fasta", "rU")  #for an obscur reason, handle_fasta had to be closed and reopened...


tot_a = 0
tot_c = 0
tot_d = 0
tot_e = 0
tot_f = 0
tot_g = 0
tot_h = 0
tot_i = 0
tot_k = 0
tot_l = 0
tot_m = 0
tot_n = 0
tot_p = 0
tot_q = 0
tot_r = 0
tot_s = 0
tot_t = 0
tot_v = 0
tot_w = 0
tot_y = 0


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
file_write.write("all amino-acids\n") # this column will be used to fill in the total counts of amino acids


for record in SeqIO.parse(handle_fasta, "fasta") : # for each uniprot in the fasta file : (for another obscur reason, the 2 for loops have to go in this order)
	my_seq = record.seq #get each sequence as a string
	X = ProtParam.ProteinAnalysis("%s" % my_seq)  # creates the variable containing the sequence, to be passed in the count_amino_acids function.
	dict_seq = X.count_amino_acids() # count the occurence of amino acid in each sequence and return as a dict
	#len_seq = len(my_seq)

	#print the list of proteins (containing uniprot id) in the first column
	file_write.write("%s\t" % record.id ) 
	
	#store the cumulated sums of each amino-acids
	tot_a = tot_a + float(dict_seq["A"])
	tot_c = tot_c + float(dict_seq["C"])
	tot_d = tot_d + float(dict_seq["D"])
	tot_e = tot_e + float(dict_seq["E"])
	tot_f = tot_f + float(dict_seq["F"])
	tot_g = tot_g + float(dict_seq["G"])
	tot_h = tot_h + float(dict_seq["H"])
	tot_i = tot_i + float(dict_seq["I"])
	tot_k = tot_k + float(dict_seq["K"])
	tot_l = tot_l + float(dict_seq["L"])
	tot_m = tot_m + float(dict_seq["M"])
	tot_n = tot_n + float(dict_seq["N"])
	tot_p = tot_p + float(dict_seq["P"])
	tot_q = tot_q + float(dict_seq["Q"])
	tot_r = tot_r + float(dict_seq["R"])
	tot_s = tot_s + float(dict_seq["S"])
	tot_t = tot_t + float(dict_seq["T"])
	tot_v = tot_v + float(dict_seq["V"])
	tot_w = tot_w + float(dict_seq["W"])
	tot_y = tot_y + float(dict_seq["Y"])
	
	#print the counts of each amino acid of each protein as a table
	tot_amino_acid_per_protein = float(dict_seq["A"]) + float(dict_seq["C"]) + float(dict_seq["D"]) + float(dict_seq["E"]) + float(dict_seq["F"]) + float(dict_seq["G"]) + float(dict_seq["H"]) + float(dict_seq["I"]) + float(dict_seq["K"]) + float(dict_seq["L"]) + float(dict_seq["M"]) + float(dict_seq["N"]) + float(dict_seq["P"]) + float(dict_seq["Q"]) + float(dict_seq["R"]) + float(dict_seq["S"]) + float(dict_seq["T"]) + float(dict_seq["V"]) + float(dict_seq["W"]) + float(dict_seq["Y"])
	# this tot_amino_acid_per_protein is to fill the last column. We could have used len_seq but since it can differ from tot_amino_acid_per_protein, it has been seen with nacho that we would use the total count instead of sequence length.
	file_write.write("%f\t" % float(dict_seq["A"]))
	file_write.write("%f\t" % float(dict_seq["C"]))
	file_write.write("%f\t" % float(dict_seq["D"]))
	file_write.write("%f\t" % float(dict_seq["E"]))
	file_write.write("%f\t" % float(dict_seq["F"]))
	file_write.write("%f\t" % float(dict_seq["G"]))
	file_write.write("%f\t" % float(dict_seq["H"]))
	file_write.write("%f\t" % float(dict_seq["I"]))
	file_write.write("%f\t" % float(dict_seq["K"]))
	file_write.write("%f\t" % float(dict_seq["L"]))
	file_write.write("%f\t" % float(dict_seq["M"]))
	file_write.write("%f\t" % float(dict_seq["N"]))
	file_write.write("%f\t" % float(dict_seq["P"]))
	file_write.write("%f\t" % float(dict_seq["Q"]))
	file_write.write("%f\t" % float(dict_seq["R"]))
	file_write.write("%f\t" % float(dict_seq["S"]))
	file_write.write("%f\t" % float(dict_seq["T"]))
	file_write.write("%f\t" % float(dict_seq["V"]))
	file_write.write("%f\t" % float(dict_seq["W"]))
	file_write.write("%f\t" % float(dict_seq["Y"]))
	file_write.write("%f\t" % float(tot_amino_acid_per_protein))
	
	
tot_amino_acid = tot_a + tot_c + tot_d + tot_e + tot_f + tot_g + tot_h + tot_i + tot_k + tot_l + tot_m + tot_n + tot_p + tot_q + tot_r + tot_s + tot_t + tot_v + tot_w + tot_y
	
file_write.write("\tTotal\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n" % (tot_a, tot_c, tot_d, tot_e, tot_f, tot_g, tot_h, tot_i, tot_k, tot_l, tot_m, tot_n, tot_p, tot_q, tot_r, tot_s, tot_t, tot_v, tot_w, tot_y, tot_amino_acid))
file_write.write("\tFraction\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n" % (tot_a/tot_amino_acid, tot_c/tot_amino_acid, tot_d/tot_amino_acid, tot_e/tot_amino_acid, tot_f/tot_amino_acid, tot_g/tot_amino_acid, tot_h/tot_amino_acid, tot_i/tot_amino_acid, tot_k/tot_amino_acid, tot_l/tot_amino_acid, tot_m/tot_amino_acid, tot_n/tot_amino_acid, tot_p/tot_amino_acid, tot_q/tot_amino_acid, tot_r/tot_amino_acid, tot_s/tot_amino_acid, tot_t/tot_amino_acid, tot_v/tot_amino_acid, tot_w/tot_amino_acid, tot_y/tot_amino_acid, tot_amino_acid/tot_amino_acid))	
	
handle_fasta.close()

for record in SeqIO.parse(handle_fasta, "fasta"):
	
	my_seq = record.seq #get each sequence as a string
	X = ProtParam.ProteinAnalysis("%s" % my_seq)  # creates the variable containing the sequence, to be passed in the count_amino_acids function.
	dict_seq = X.count_amino_acids() # will count frequency of amino acid in each sequence and return as a dict
	file_write5.write("%s\t" % record.id ) # will print the list of proteins in a column
	print record
	
	file_write5.write("%s\t" % dict_seq["A"]) # will print the counts of each amino acid of each protein as a table
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
	file_write5.write("\n")

handle_fasta.close()


stop = timeit.default_timer()
print stop - start 