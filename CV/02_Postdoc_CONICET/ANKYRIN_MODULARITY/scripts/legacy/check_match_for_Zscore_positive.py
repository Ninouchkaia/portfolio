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





elm_domain_mapping_dict = {}
domain_elm_mapping_dict = {}
with open("elm_interaction_domains_modified.txt", 'rU') as file_open :
	data = file_open.readlines()
	for line in data[1:] :
		line = line.split("\t")
		elm_name = line[0]
		pfam_family = line[2]
		elm_domain_mapping_dict[elm_name] = pfam_family		
		domain_elm_mapping_dict[pfam_family] = elm_name

# elm_clan_mapping_dict = {}
# clan_elm_mapping_dict = {}
# with open("elm_interaction_domains_modified.txt", 'rU') as file_open :
# 	data = file_open.readlines()
# 	for line in data[1:] :
# 		line = line.split("\t")
# 		elm_name = line[0]
# 		pfam_family = line[2]
# 		elm_domain_mapping_dict[elm_name] = pfam_family
# 		elm_domain_mapping_dict[elm_name] = pfam_family


elms_with_Zscore_positive = []
elms_with_Zscore_negative = []
elms_in_ank = []
with open("elms_counts_in_ank_proteins_with_conservation_Zscores_colored.txt", 'rU') as file_open :
	data = file_open.readlines()
	for line in data[1:] :
		line = line.split("\t")
		Zscore = float(line[-2])
		if Zscore > 0 :
			elm_name = line[0]
			elms_with_Zscore_positive.append(elm_name)
		if Zscore <= 0 :
			elm_name = line[0]
			elms_with_Zscore_negative.append(elm_name)
		elm_name = line[0]
		elms_in_ank.append(elm_name)






interacting_domains_corresponding_to_elms_with_Zscore_positive = []
# with open("elms_in_ank_with_Zscore_positive.txt", 'a+') as file_write :
for elm_name in elms_with_Zscore_positive :
	if elm_name in elm_domain_mapping_dict :
		print elm_name, elm_domain_mapping_dict[elm_name]
		interacting_domains_corresponding_to_elms_with_Zscore_positive.append(elm_domain_mapping_dict[elm_name])
	else : 
		print elm_name, "doesnt have an interacting domain"
interacting_domains_corresponding_to_elms_with_Zscore_positive = list(set(interacting_domains_corresponding_to_elms_with_Zscore_positive))
print "len(interacting_domains_corresponding_to_elms_with_Zscore_positive)", len(interacting_domains_corresponding_to_elms_with_Zscore_positive)

# interacting_domains_corresponding_to_elms_with_Zscore_negative = []
# for elm_name in elms_with_Zscore_negative :
# 	if elm_name in elm_domain_mapping_dict :
# 		print elm_name, elm_domain_mapping_dict[elm_name]
# 		interacting_domains_corresponding_to_elms_with_Zscore_negative.append(elm_domain_mapping_dict[elm_name])
# 	else : 
# 		print elm_name, "doesnt have an interacting domain"
# interacting_domains_corresponding_to_elms_with_Zscore_negative = list(set(interacting_domains_corresponding_to_elms_with_Zscore_negative))
# print "len(interacting_domains_corresponding_to_elms_with_Zscore_negative)", len(interacting_domains_corresponding_to_elms_with_Zscore_negative)


# interacting_domains_corresponding_to_elms = []
# for elm_name in elms_in_ank :
# 	if elm_name in elm_domain_mapping_dict :
# 		print elm_name, elm_domain_mapping_dict[elm_name]
# 		interacting_domains_corresponding_to_elms.append(elm_domain_mapping_dict[elm_name])
# 	else : 
# 		print elm_name, "doesnt have an interacting domain"
# interacting_domains_corresponding_to_elms = list(set(interacting_domains_corresponding_to_elms))
# print "len(interacting_domains_corresponding_to_elms)", len(interacting_domains_corresponding_to_elms)





domains_with_Zscore_positive = []
with open("Pfam_domains_in_binding_partners_2038.fasta_MaxHomologs_1000_Zscores_color.txt", 'rU') as file_open :
	data = file_open.readlines()
	for line in data[1:] :
		line = line.split("\t")
		Zscore = float(line[-3])
		if Zscore > 0 :
			domain_name = line[0]
			domains_with_Zscore_positive.append(domain_name)
print "len(domains_with_Zscore_positive)", len(domains_with_Zscore_positive)

# domains_with_Zscore_negative = []
# with open("Pfam_domains_in_binding_partners_2038.fasta_MaxHomologs_1000_Zscores_color.txt", 'rU') as file_open :
# 	data = file_open.readlines()
# 	for line in data[1:] :
# 		line = line.split("\t")
# 		Zscore = float(line[-3])
# 		if Zscore <= 0 :
# 			domain_name = line[0]
# 			domains_with_Zscore_negative.append(domain_name)
# print "len(domains_with_Zscore_negative)", len(domains_with_Zscore_negative)


# domains_in_partners = []
# with open("Pfam_domains_in_binding_partners_2038.fasta_MaxHomologs_1000_Zscores_color.txt", 'rU') as file_open :
# 	data = file_open.readlines()
# 	for line in data[1:] :
# 		line = line.split("\t")
# 		domain_name = line[0]
# 		domains_in_partners.append(domain_name)
# print "len(domains_in_partners)", len(domains_in_partners)





domains_that_are_enriched_in_partners_and_match_enriched_elms_in_ank = []
for domain_name in domains_with_Zscore_positive :
	if domain_name in interacting_domains_corresponding_to_elms_with_Zscore_positive :
		domains_that_are_enriched_in_partners_and_match_enriched_elms_in_ank.append(domain_name)
		print domain_name, " : ", domain_elm_mapping_dict[domain_name]
print len(domains_that_are_enriched_in_partners_and_match_enriched_elms_in_ank)
			 
# domains_that_are_depleted_in_partners_and_match_depleted_elms_in_ank = []
# for domain_name in domains_with_Zscore_negative :
# 	if domain_name in interacting_domains_corresponding_to_elms_with_Zscore_negative :
# 		domains_that_are_depleted_in_partners_and_match_depleted_elms_in_ank.append(domain_name)
# 		print domain_name, " : ", domain_elm_mapping_dict[domain_name]
# print len(domains_that_are_depleted_in_partners_and_match_depleted_elms_in_ank)

# domains_in_partners_that_match_elms_in_ank = []
# for domain_name in domains_in_partners :
# 	if domain_name in interacting_domains_corresponding_to_elms :
# 		domains_in_partners_that_match_elms_in_ank.append(domain_name)
# 		print domain_name, " : ", domain_elm_mapping_dict[domain_name]
# print len(domains_in_partners_that_match_elms_in_ank)










stop = timeit.default_timer()
print stop - start 