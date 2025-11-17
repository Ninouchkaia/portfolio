#!/usr/bin/python
# -*- coding: utf-8 -*-

import timeit
start = timeit.default_timer()

from Bio import SeqIO
import sys


# if sys.argv[1] == "non_ankregions" : declare large linker size threshold as sys.argv[2]
# if sys.argv[1] == "ankregions" : declare large linker size threshold as sys.argv[2]
# if sys.argv[1] == "ank_linkers" : declare large linker size threshold as sys.argv[2]
# if sys.argv[1] == "all_ank_linkers"
# if sys.argv[1] == "ank_repeats"
# if sys.argv[1] == "ank_proteins"

#declare fasta ank_proteins according to uniprot
#and build sequence_dict :
ank_proteins_from_uniprot_fasta = "all-ank-20130926.fasta"
sequence_dict = {}
uniprot_list_total = []
with open(ank_proteins_from_uniprot_fasta, 'rU') as fasta_handle :
	for record in SeqIO.parse(fasta_handle, "fasta") :
		uniprot_id = record.id[3:9]
		uniprot_list_total.append(uniprot_id)
		uniprot_seq = str(record.seq)
		sequence_dict[uniprot_id] = uniprot_seq

#get the unp list of ank proteins for which rocio hmmer search detected repeats. (there are 1191)
ankyrin_list = [] 
with open("uniprot_hmmer1_results_for_all_ank_1238.txt", 'rU') as file_read :
	data = file_read.readlines()
	for line in data :	
		uniprot_id = line.replace("\n","")
		ankyrin_list.append(uniprot_id) 



if sys.argv[1] == "ank_proteins" :
	#build fasta ank_proteins according to rocio's detection:
	with open(ank_proteins_from_uniprot_fasta, 'rU') as fasta_handle :
		with open("ankproteins_1191.fasta", 'a+') as file_write :
			for uniprot_id in ankyrin_list :
				sequence = sequence_dict[uniprot_id]
				file_write.write(">sp|%s|%s\n" % (uniprot_id, len(sequence)))
				file_write.write("%s\n" % sequence)


if sys.argv[1] == "ank_repeats" :
	#build fasta for repeats:
	with open("ankrepeats_1191.fasta", 'a+') as file_write :
		for uniprot_id in ankyrin_list :
			with open("/Users/verstrat/temp_qb/paper_elm_2014/ank_repeats_analysis/compute_ank1238/%s_compute_ank.txt" % uniprot_id, 'rU') as file_read :
				data = file_read.readlines()
				for i in range(0, len(data)) :
					data[i] = data[i].split("\t")
					sequence_start = int(data[i][2])
					sequence_end = int(data[i][3])
					sequence = sequence_dict[uniprot_id][sequence_start-1:sequence_end]
					file_write.write(">sp|%s|%s|ank_repeat_%s|%s|%s|%s\n" % (uniprot_id, len(sequence_dict[uniprot_id]), i+1, sequence_start, sequence_end, len(sequence)))
					file_write.write("%s\n" % sequence)


if sys.argv[1] == "all_ank_linkers" :
	#build fasta for linkers (=between repeats OF ALL SIZES)):
	with open("anklinkers_1191.fasta", 'a+') as file_write :
		for uniprot_id in ankyrin_list :
			with open("/Users/verstrat/temp_qb/paper_elm_2014/ank_repeats_analysis/compute_linkers1191_correct/%s_anklinkers.txt" % uniprot_id, 'rU') as file_read :
				data = file_read.readlines()
				for i in range(0, len(data)) :
					data[i] = data[i].split("\t")
					anklinker_length = int(data[i][4].replace("\n",""))
					sequence_start = int(data[i][2])
					sequence_end = int(data[i][3])
					sequence = sequence_dict[uniprot_id][sequence_start-1:sequence_end]
					file_write.write(">sp|%s|%s|ank_linker_%s|%s|%s|%s\n" % (uniprot_id, len(sequence_dict[uniprot_id]), i+1, sequence_start, sequence_end, len(sequence)))
					file_write.write("%s\n" % sequence)


if sys.argv[1] == "ank_linkers" :
	threshold = int(sys.argv[2])
	#build fasta for linkers (=between repeats AND < 132aa)):
	with open("anklinkers_1191_threshold%s.fasta" % threshold, 'a+') as file_write :
		for uniprot_id in ankyrin_list :
			with open("/Users/verstrat/temp_qb/paper_elm_2014/ank_repeats_analysis/compute_linkers1191_correct/%s_anklinkers.txt" % uniprot_id, 'rU') as file_read :
				data = file_read.readlines()
				for i in range(0, len(data)) :
					data[i] = data[i].split("\t")
					anklinker_length = int(data[i][4].replace("\n",""))
					if anklinker_length < threshold :
						sequence_start = int(data[i][2])
						sequence_end = int(data[i][3])
						sequence = sequence_dict[uniprot_id][sequence_start-1:sequence_end]
						file_write.write(">sp|%s|%s|ank_linker_%s|%s|%s|%s\n" % (uniprot_id, len(sequence_dict[uniprot_id]), i+1, sequence_start, sequence_end, len(sequence)))
						file_write.write("%s\n" % sequence)

if sys.argv[1] == "ankregions" :
	threshold = int(sys.argv[2])
	#build fasta for regions:
	with open("ankregions_1191_threshold%s.fasta" % threshold, 'a+') as file_write :
		for uniprot_id in ankyrin_list :
			with open("/Users/verstrat/temp_qb/paper_elm_2014/ank_repeats_analysis/compute_ankregions1191_threshold%s/%s_ankregions.txt" % (threshold, uniprot_id), 'rU') as file_read :
				data = file_read.readlines()
				for i in range(0, len(data)) :
					data[i] = data[i].split("\t")
					sequence_start = int(data[i][2])
					sequence_end = int(data[i][3])
					sequence = sequence_dict[uniprot_id][sequence_start-1:sequence_end]
					file_write.write(">sp|%s|%s|ankregion_%s|%s|%s|%s\n" % (uniprot_id, len(sequence_dict[uniprot_id]), i+1, sequence_start, sequence_end, len(sequence)))
					file_write.write("%s\n" % sequence)

if sys.argv[1] == "non_ankregions" :
	threshold = int(sys.argv[2])
	#build fasta for regions:
	with open("non_ankregions_1191_threshold%s.fasta" % threshold, 'a+') as file_write :
		for uniprot_id in ankyrin_list :
			with open("/Users/verstrat/temp_qb/paper_elm_2014/ank_repeats_analysis/compute_non_ankregions1191_threshold%s/%s_non_ankregions.txt" % (threshold,uniprot_id), 'rU') as file_read :
				data = file_read.readlines()
				for i in range(0, len(data)) :
					data[i] = data[i].split("\t")
					sequence_start = int(data[i][2])
					sequence_end = int(data[i][3])
					sequence = sequence_dict[uniprot_id][sequence_start-1:sequence_end]
					file_write.write(">sp|%s|%s|non_ankregion_%s|%s|%s|%s\n" % (uniprot_id, len(sequence_dict[uniprot_id]), i+1, sequence_start, sequence_end, len(sequence)))
					file_write.write("%s\n" % sequence)


stop = timeit.default_timer()
print stop - start 
