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
# from Bio import SeqIO # to parse the fasta file
import collections
import math

elm_domain_and_clan_mapping_dict = {}
domain_and_clan_elm_mapping_dict = {}

with open("elm_interaction_clans_update.txt", 'rU') as file_open :
	data = file_open.readlines()
	for line in data[1:] :
		line = line.split("\t")
		elm_name = line[0]
		if line[-1] == '\N\n' : 
			domain_or_clan = line[2]
		else : 
			domain_or_clan = line[-1].replace("\n","")

		if elm_name not in elm_domain_and_clan_mapping_dict :
			elm_domain_and_clan_mapping_dict[elm_name] = [domain_or_clan]
		else :
			elm_domain_and_clan_mapping_dict[elm_name].append(domain_or_clan)

		if domain_or_clan not in domain_and_clan_elm_mapping_dict :
			domain_and_clan_elm_mapping_dict[domain_or_clan] = [elm_name]
		else :
			domain_and_clan_elm_mapping_dict[domain_or_clan].append(elm_name)

print len(domain_and_clan_elm_mapping_dict), len(elm_domain_and_clan_mapping_dict)

# for i in domain_and_clan_elm_mapping_dict :
# 	print i, domain_and_clan_elm_mapping_dict[i]

# print domain_and_clan_elm_mapping_dict["CL0020"]

# with open("domains_clans_to_elms_mapping.txt", 'a+') as file_write :
# 	for i in domain_and_clan_elm_mapping_dict :
# 		print i, domain_and_clan_elm_mapping_dict[i]
# 		for j in range(0, len(domain_and_clan_elm_mapping_dict[i])) :
# 			file_write.write("%s\t%s\n" % (i, domain_and_clan_elm_mapping_dict[i][j]))

# with open("elms_to_domains_clans_mapping.txt", 'a+') as file_write :
# 	for i in elm_domain_and_clan_mapping_dict :
# 		print i, elm_domain_and_clan_mapping_dict[i]
# 		for j in range(0, len(elm_domain_and_clan_mapping_dict[i])) :
# 			file_write.write("%s\t%s\n" % (i, elm_domain_and_clan_mapping_dict[i][j]))


# filename = sys.argv[1]

# with open("counts/%s_binding_MAPPED.txt" % filename[:-4], 'a+') as file_write :
# 	with open("counts/%s" % filename, 'rU') as file_open :
# 		data = file_open.readlines()
# 		first_line_splitted = data[0].split("\t")
# 		file_write.write("%s\tbinding_domain_or_clan\t%s" % (first_line_splitted[0], '\t'.join(first_line_splitted[1:])))
# 		for line in data[1:] :
# 			line = line.split("\t")
# 			elm_name = line[0]
# 			if elm_name in elm_domain_and_clan_mapping_dict :
# 				if len(elm_domain_and_clan_mapping_dict[elm_name]) == 1 :
# 					binding_domain_or_clan = elm_domain_and_clan_mapping_dict[elm_name][0]
# 				else :
# 					binding_domain_or_clan = ', '.join(elm_domain_and_clan_mapping_dict[elm_name])
# 				print elm_name, binding_domain_or_clan
# 				file_write.write("%s\t%s\t%s" % (elm_name, binding_domain_or_clan, '\t'.join(line[1:])))
# 			else :
# 				print elm_name, "no binding domain or clan"
# 				file_write.write("%s\t%s\t%s" % (elm_name, "NULL", '\t'.join(line[1:])))



# filename = sys.argv[1]

# with open("counts/%s_binding_MAPPED.txt" % filename[:-4], 'a+') as file_write :
# 	with open("counts/%s" % filename, 'rU') as file_open :
# 		data = file_open.readlines()
# 		first_line_splitted = data[0].split("\t")
# 		file_write.write("%s\t%s\tbinding_elm\t%s" % ("domain", "domain_or_clan", '\t'.join(first_line_splitted[2:])))
# 		for line in data[1:] :
# 			line = line.split("\t")
# 			clan_name = line[1]
# 			domain_name = line[0]
# 			if clan_name == 'NULL' :
# 				domain_or_clan_name = line[0]
# 			else : 
# 				domain_or_clan_name = line[1]

# 			if domain_or_clan_name in domain_and_clan_elm_mapping_dict :
# 				if len(domain_and_clan_elm_mapping_dict[domain_or_clan_name]) == 1 :
# 					binding_elm = domain_and_clan_elm_mapping_dict[domain_or_clan_name][0]
# 				else :
# 					binding_elm = ', '.join(domain_and_clan_elm_mapping_dict[domain_or_clan_name])


# 				file_write.write("%s\t%s\t%s\t%s" % (domain_name, domain_or_clan_name, binding_elm, '\t'.join(line[2:])))

# 			else :
# 				#print domain_name, "no binding elm"
# 				file_write.write("%s\t%s\t%s\t%s" % (domain_name, domain_or_clan_name, "NULL", '\t'.join(line[2:])))








stop = timeit.default_timer()
print stop - start 