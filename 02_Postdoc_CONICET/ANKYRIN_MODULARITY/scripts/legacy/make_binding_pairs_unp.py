#!/usr/bin/python
# -*- coding: utf-8 -*-

import timeit
start = timeit.default_timer()
import sys
import os
import glob
import re
from commands import getoutput # permet d'obtenir l'output d'une commande bash.
from numpy import prod
from Bio import SeqIO # to parse the fasta file
import collections
import math


# binding_partners_mapping_dict, ank_mapping_dict, merged_mapping_dict = {}, {}, {}

# with open("mapping_table_unp_string_uniref50_UPPERCASE.txt", 'rU') as file_open :
# 	data = file_open.readlines()
# 	for line in data :
# 		line = line.split("\t")
# 		string_id = line[0].upper()
# 		unp_id = line[1].replace("\n","")
# 		ank_mapping_dict[string_id] = unp_id
# 		merged_mapping_dict[string_id] = unp_id

# with open("binding_partners_mapping_fused12_UPPERCASE.txt", 'rU') as file_open :
# 	data = file_open.readlines()
# 	for line in data :
# 		line = line.split("\t")
# 		string_id = line[0].upper()
# 		unp_id = line[1].replace("\n","")
# 		binding_partners_mapping_dict[string_id] = unp_id
# 		if string_id not in merged_mapping_dict :
# 			merged_mapping_dict[string_id] = unp_id
# 		else :
# 			print string_id

# print len(ank_mapping_dict), len(binding_partners_mapping_dict), len(merged_mapping_dict)


# partnerA_list, partnerB_list, indep_proteins = [], [], []

# with open("binding_uniref50_score500_reduced.txt", 'rU') as file_open :
# 	data = file_open.readlines()
# 	with open("binding_uniref50_score500_reduced_unp_MAPPED.txt", 'a+') as file_write :
# 		with open("unmapped_binding_partners.txt", 'a+') as file_write2 :
# 			for line in data :
# 				line = line.split("\t")
# 				partnerA = line[0].upper()
# 				partnerA_list.append(partnerA)
# 				partnerB = line[1].upper()
# 				partnerB_list.append(partnerB)
# 				indep_proteins.append(partnerA)
# 				indep_proteins.append(partnerB)
# 				binding_score = line[5]
# 				if line[6] != '' :
# 					source = line[6].replace("\n","")
# 				else :
# 					source = line[7].replace("\n","")

# 				if partnerA in merged_mapping_dict :
# 					partnerA_unp = merged_mapping_dict[partnerA]
# 				else :
# 					file_write2.write("%s (unmapped)\t%s\t%s\t%s\n" % (partnerA, partnerB, binding_score, source))
				
# 				if partnerB in merged_mapping_dict :
# 					partnerB_unp = merged_mapping_dict[partnerB]
# 				else :
# 					file_write2.write("%s\t%s (unmapped)\t%s\t%s\n" % (partnerA, partnerB, binding_score, source))

# 				if partnerA not in merged_mapping_dict and partnerB not in merged_mapping_dict :
# 					print "REALLY?"

# 				if partnerA in merged_mapping_dict and partnerB in merged_mapping_dict :
# 					file_write.write("%s\t%s\t%s\t%s\n" % (partnerA_unp, partnerB_unp, binding_score, source))



# indep_proteins = list(set(indep_proteins))	

# print len(partnerA_list), len(set(partnerA_list))
# print len(partnerB_list), len(set(partnerB_list))
# print len(indep_proteins)


with open("binding_uniref50_score500_reduced_unp_MAPPED.txt", 'rU') as file_open :
	data = file_open.readlines()
	with open("interacting_pairs_list.txt", 'a+') as file_write :
		for line in data :
			line = line.split("\t")
			partnerA = line[0]
			partnerB = line[1]
			interacting_pair = "%s_%s" % (partnerA, partnerB)
			file_write.write("%s\n" % interacting_pair)



stop = timeit.default_timer()
print stop - start 