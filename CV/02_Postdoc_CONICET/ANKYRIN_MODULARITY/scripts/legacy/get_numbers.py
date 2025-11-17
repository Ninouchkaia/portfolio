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


# filename = sys.argv[1]
# # filename2 = sys.argv[2]

# elm_list = []

# # general counting
# with open(filename, 'rU') as file_open :
# 	data = file_open.readlines()
# 	for line in data[1:] :
# 		line = line.split("\t")
# 		elm_name = line[0]
# 		elm_list.append(elm_name)
# elm_list = list(set(elm_list))
# print len(elm_list)

# clv, deg, doc, lig, mod, trg = 0,0,0,0,0,0
# for elm_name in elm_list :
# 	if elm_name[:3] == 'CLV' :
# 		clv+=1
# 	elif elm_name[:3] == 'DEG' :
# 		deg+=1
# 	elif elm_name[:3] == 'DOC' :
# 		doc+=1
# 	elif elm_name[:3] == 'LIG' :
# 		lig+=1
# 	elif elm_name[:3] == 'MOD' :
# 		mod+=1
# 	elif elm_name[:3] == 'TRG' :
# 		trg+=1

# # enriched elms
# with open(filename, 'rU') as file_open :
# 	data = file_open.readlines()
# 	for line in data[1:] :
# 		line = line.split("\t")
# 		elm_name = line[0]
# 		enrichment = float(line[5])
# 		if enrichment >= 1 :
# 			elm_list.append(elm_name)
# elm_list = list(set(elm_list))
# print len(elm_list)

# clv, deg, doc, lig, mod, trg = 0,0,0,0,0,0
# for elm_name in elm_list :
# 	if elm_name[:3] == 'CLV' :
# 		clv+=1
# 	elif elm_name[:3] == 'DEG' :
# 		deg+=1
# 	elif elm_name[:3] == 'DOC' :
# 		doc+=1
# 	elif elm_name[:3] == 'LIG' :
# 		lig+=1
# 	elif elm_name[:3] == 'MOD' :
# 		mod+=1
# 	elif elm_name[:3] == 'TRG' :
# 		trg+=1

# # enriched and conserved elms
# with open(filename, 'rU') as file_open :
# 	data = file_open.readlines()
# 	for line in data[1:] :
# 		line = line.split("\t")
# 		elm_name = line[0]
# 		enrichment = float(line[5]) 
# 		if enrichment >= 1 :
# 			color = line[-1].replace("\n","")
# 			if color == '\"red\"' :
# 				elm_list.append(elm_name)

# elm_list = list(set(elm_list))
# print len(elm_list)
# print "%s" % (", ". join(elm_list))

# clv, deg, doc, lig, mod, trg = 0,0,0,0,0,0
# for elm_name in elm_list :
# 	if elm_name[:3] == 'CLV' :
# 		clv+=1
# 	elif elm_name[:3] == 'DEG' :
# 		deg+=1
# 	elif elm_name[:3] == 'DOC' :
# 		doc+=1
# 	elif elm_name[:3] == 'LIG' :
# 		lig+=1
# 	elif elm_name[:3] == 'MOD' :
# 		mod+=1
# 	elif elm_name[:3] == 'TRG' :
# 		trg+=1

# print "clv, deg, doc, lig, mod, trg"
# print clv, deg, doc, lig, mod, trg


# # extract some elm domain mapping and store binding domains
# elm_list = ["LIG_WRPW_1", "MOD_SPalmitoyl_4", "LIG_TPR", "MOD_ASX_betaOH_EGF", "LIG_WH1", "LIG_PCNA_PIPBox_1", "MOD_TYR_CSK", "LIG_Sin3_1", "DOC_AGCK_PIF_1"]
# binding_domain_list = []
# with open(filename, 'rU') as file_open :
# 	data = file_open.readlines()
# 	for line in data[1:] :
# 		line = line.split("\t")
# 		elm_name = line[0]
# 		binding_domain = line[2]
# 		if elm_name in elm_list :
# 			print "%s" % ("\t".join(line))
# 			binding_domain_list.append(binding_domain)

# binding_domain_list = list(set(binding_domain_list))
# print len(binding_domain_list)
# print "\n"
# # extract some domains binding certain elms
# with open(filename2, 'rU') as file_open :
# 	data = file_open.readlines()
# 	for line in data[1:] :
# 		line = line.split("\t")
# 		domain_name = line[0]
# 		if domain_name in binding_domain_list : 
# 			print "%s" % ("\t".join(line))

# # extract some elm clan/domain mapping and store binding clan/domain
# elm_list = ["LIG_WRPW_1", "MOD_SPalmitoyl_4", "LIG_TPR", "MOD_ASX_betaOH_EGF", "LIG_WH1", "LIG_PCNA_PIPBox_1", "MOD_TYR_CSK", "LIG_Sin3_1", "DOC_AGCK_PIF_1"]
# binding_clan_or_domain_list = []
# with open(filename, 'rU') as file_open :
# 	data = file_open.readlines()
# 	for line in data[1:] :
# 		line = line.split("\t")
# 		elm_name = line[0]
# 		binding_domain_or_clan = line[1].replace("\n","")
# 		if elm_name in elm_list :
# 			print "%s" % ("\t".join(line))
# 			binding_clan_or_domain_list.append(binding_domain_or_clan)

# binding_clan_or_domain_list = list(set(binding_clan_or_domain_list))
# print len(binding_clan_or_domain_list)
# print "\n"
# # extract some domains binding certain elms
# with open(filename2, 'rU') as file_open :
# 	data = file_open.readlines()
# 	for line in data[1:] :
# 		line = line.split("\t")
# 		domain_or_clan_name = line[1].replace("\n","")
# 		if domain_or_clan_name in binding_clan_or_domain_list : 
# 			print "%s" % ("\t".join(line[:4]))


# # get Zscore for a list of specific domains :
# domains_data = sys.argv[1]
# domains_list = ["zf-DHHC",
# "DZR",
# "zf-RanBP",
# "zf-NADH-PPase",
# "PH",
# "GRAM",
# "PH_8",
# "PID",
# "Pkinase_Tyr",
# "Pkinase",
# "TPR_8",
# "TPR_2",
# "Arm",
# "TPR_11",
# "TPR_7",
# "RCC1"]
# with open(domains_data, 'rU') as file_open :
# 	data = file_open.readlines()
# 	for domain_name_from_list in domains_list :

# 		for line in data[1:] :
# 			line = line.split("\t")
# 			domain_name = line[0]
# 			Zscore = float(line[4])
# 			if domain_name == domain_name_from_list :
# 				print line[0], line[1], line[3], line[4]


# get enrichment score for a list of specific domains :
domains_data = sys.argv[1]
domains_list = ['ZU5', 'DBB', 'TRP_2', 'Pre-SET', 'Death', 'Ion_trans', 'EGF', 'Ank_3', 'HECT', 'ZZ', 'PH', 'DUF3354', 'hEGF', 'SH3_2', 'NODP', 'SAM_2', 'Notch', 'PID', 'CAP_GLY', 'DUF3454', 'zf-RanBP', 'RHD']

with open(domains_data, 'rU') as file_open :
	data = file_open.readlines()
	for domain_name_from_list in domains_list :

		for line in data[1:] :
			line = line.split("\t")
			domain_name = line[0]
			enrichment = float(line[3])
			if domain_name == domain_name_from_list :
				print "\t".join(line[0:-2])

# get mapping clans from ELM list : 
# elms_data = sys.argv[1]
# elms_list = ["LIG_EH1_1",
# "LIG_FAT_LD_1",
# "DEG_MDM2_1",
# "LIG_NRBOX",
# "MOD_OGLYCOS",
# "TRG_ER_FFAT_1",
# "DOC_AGCK_PIF_1",
# "MOD_ASX_betaOH_EGF",
# "LIG_SPRY_1"]
# with open(elms_data, 'rU') as file_open :
# 	data = file_open.readlines()
# 	for line in data[1:] :
# 		line = line.split("\t")
# 		elm_name = line[0]
# 		clan_or_domain_name = line[1].replace("\n","")
# 		if elm_name in elms_list :
# 			print elm_name, clan_or_domain_name



stop = timeit.default_timer()
print stop - start 