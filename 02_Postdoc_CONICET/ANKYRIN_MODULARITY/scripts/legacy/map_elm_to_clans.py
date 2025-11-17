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


pfam_family_with_no_clan = []
pfam_family_list = []
pfam_clan_list = []
family_clan_mapping_dict = {}
family_clan_mapping_dict_update = {}
with open("Pfam-A.clans.txt", 'rU') as file_open :
	data = file_open.readlines()
	for line in data[1:] :
		line = line.split("\t")
		pfam_family = line[0]
		pfam_family_list.append(pfam_family)
		pfam_clan = line[2]
		pfam_clan_ID = line[1]
		if pfam_clan != '\N' :
			pfam_clan_list.append(pfam_clan)
		if pfam_clan == '\N' :
			pfam_family_with_no_clan.append(pfam_family)
		family_clan_mapping_dict[pfam_family] = pfam_clan
		family_clan_mapping_dict_update[pfam_family] = pfam_clan_ID

print len(pfam_family_list)
print len(family_clan_mapping_dict)
# print family_clan_mapping_dict
print len(pfam_clan_list)
print len(set(pfam_clan_list))

print len(pfam_family_with_no_clan)

# with open("elm_interaction_domains_modified.txt", 'rU') as file_open :
# 	data = file_open.readlines()
# 	with open("elm_interaction_clans.txt", 'a+') as file_write :
# 		file_write.write("%s\tPfam_clan\n" % data[0].replace("\n",""))
# 		for line in data[1:] :
# 			line = line.split("\t")
# 			elm_name = line[0]
# 			pfam_family = line[1]
# 			pfam_clan = family_clan_mapping_dict[pfam_family]
# 			file_write.write("%s\t%s\t%s\n" % ('\t'.join(line[:-1]), line[-1].replace("\n",""), pfam_clan))

with open("elm_interaction_domains_modified.txt", 'rU') as file_open :
	data = file_open.readlines()
	with open("elm_interaction_clans_update.txt", 'a+') as file_write :
		file_write.write("%s\tPfam_clan_ID\n" % data[0].replace("\n",""))
		for line in data[1:] :
			line = line.split("\t")
			elm_name = line[0]
			pfam_family = line[1]
			pfam_clan_ID = family_clan_mapping_dict_update[pfam_family]
			file_write.write("%s\t%s\t%s\n" % ('\t'.join(line[:-1]), line[-1].replace("\n",""), pfam_clan_ID))





stop = timeit.default_timer()
print stop - start 