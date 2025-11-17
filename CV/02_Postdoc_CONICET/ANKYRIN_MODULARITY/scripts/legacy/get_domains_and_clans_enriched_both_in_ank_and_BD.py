#!/usr/bin/python
# -*- coding: utf-8 -*-

import timeit
start = timeit.default_timer()
import sys
import os
import glob
import re
from commands import getoutput # permet d'obtenir l'output d'une commande bash.
import numpy
# from Bio import SeqIO # to parse the fasta file
import collections
import math

domains_in_ank = sys.argv[1]
domains_in_BD = sys.argv[2]




#####################################################################################################
###########################        #based on enrichment > average_value     #########################
#####################################################################################################


# # #average value for enrichment : 
# enrichment_values_in_ank = []
# with open(domains_in_ank, 'rU') as file_open :
# 	data = file_open.readlines()
# 	for line in data[1:] :
# 		line = line.split("\t")
# 		enrichment = float(line[3])
# 		enrichment_values_in_ank.append(enrichment)
# average_enrichment_in_ank = numpy.mean(enrichment_values_in_ank)
# print average_enrichment_in_ank #3.39631611409


# enrichment_values_in_BD = []
# with open(domains_in_BD, 'rU') as file_open :
# 	data = file_open.readlines()
# 	for line in data[1:] :
# 		line = line.split("\t")
# 		enrichment = float(line[3])
# 		enrichment_values_in_BD.append(enrichment)
# average_enrichment_in_BD = numpy.mean(enrichment_values_in_BD)
# print average_enrichment_in_BD #2.79198264968



# enriched_domains_in_ank, enriched_clans_in_ank = [], []
# with open(domains_in_ank, 'rU') as file_open :
# 	data = file_open.readlines()
# 	for line in data[1:] :
# 		line = line.split("\t")
# 		domain_name = line[0]
# 		clan_name = line[1]
# 		enrichment = float(line[3])
# 		if enrichment > average_enrichment_in_ank :
# 			enriched_domains_in_ank.append(domain_name)
# 			enriched_clans_in_ank.append(clan_name)
# enriched_domains_in_ank = list(set(enriched_domains_in_ank))
# enriched_clans_in_ank = list(set(enriched_clans_in_ank))

# enriched_domains_in_BD, enriched_clans_in_BD = [], []
# with open(domains_in_BD, 'rU') as file_open :
# 	data = file_open.readlines()
# 	for line in data[1:] :
# 		line = line.split("\t")
# 		domain_name = line[0]
# 		clan_name = line[1]
# 		enrichment = float(line[3])
# 		if enrichment > average_enrichment_in_BD :
# 			enriched_domains_in_BD.append(domain_name)
# 			enriched_clans_in_BD.append(clan_name)
# enriched_domains_in_BD = list(set(enriched_domains_in_BD))
# enriched_clans_in_BD = list(set(enriched_clans_in_BD))


# enriched_domains_in_ank_and_BD, enriched_clans_in_ank_and_BD = [], []
# for domain_name in enriched_domains_in_ank :
# 	if domain_name in enriched_domains_in_BD :
# 		enriched_domains_in_ank_and_BD.append(domain_name)
# for clan_name in enriched_clans_in_ank :
# 	if clan_name in enriched_clans_in_BD :
# 		enriched_clans_in_ank_and_BD.append(clan_name)

# print "len(enriched_domains_in_ank_and_BD)", len(enriched_domains_in_ank_and_BD) 
# print "len(enriched_clans_in_ank_and_BD)", len(enriched_clans_in_ank_and_BD)

# print enriched_domains_in_ank_and_BD
# print enriched_clans_in_ank_and_BD

# len(enriched_domains_in_ank_and_BD) 22
# len(enriched_clans_in_ank_and_BD) 26
# ['ZU5', 'DBB', 'TRP_2', 'Pre-SET', 'Death', 'Ion_trans', 'EGF', 'Ank_3', 'HECT', 'ZZ', 'PH', 'DUF3354', 'hEGF', 'SH3_2', 'NODP', 'SAM_2', 'Notch', 'PID', 'CAP_GLY', 'DUF3454', 'zf-RanBP', 'RHD']
# ['DBB', 'TRP_2', 'CL0073', 'Pre-SET', 'CL0159', 'ZU5', 'CL0229', 'CL0202', 'CL0125', 'CL0001', 'CL0552', 'CL0006', 'CL0465', 'DUF3354', 'CL0030', 'CL0266', 'CL0033', 'CL0167', 'CL0023', 'CL0041', 'CL0003', 'NODP', 'CL0010', 'Notch', 'CAP_GLY', 'DUF3454']

######################################################################################################
############################        #based on enrichment > 1      ####################################
######################################################################################################

# 
# enriched_domains_in_ank, enriched_clans_in_ank = [], []
# with open(domains_in_ank, 'rU') as file_open :
# 	data = file_open.readlines()
# 	for line in data[1:] :
# 		line = line.split("\t")
# 		domain_name = line[0]
# 		clan_name = line[1]
# 		enrichment = float(line[3])
# 		if enrichment > 1 :
# 			enriched_domains_in_ank.append(domain_name)
# 			enriched_clans_in_ank.append(clan_name)
# enriched_domains_in_ank = list(set(enriched_domains_in_ank))
# enriched_clans_in_ank = list(set(enriched_clans_in_ank))

# enriched_domains_in_BD, enriched_clans_in_BD = [], []
# with open(domains_in_BD, 'rU') as file_open :
# 	data = file_open.readlines()
# 	for line in data[1:] :
# 		line = line.split("\t")
# 		domain_name = line[0]
# 		clan_name = line[1]
# 		enrichment = float(line[3])
# 		if enrichment > 1 :
# 			enriched_domains_in_BD.append(domain_name)
# 			enriched_clans_in_BD.append(clan_name)
# enriched_domains_in_BD = list(set(enriched_domains_in_BD))
# enriched_clans_in_BD = list(set(enriched_clans_in_BD))


# enriched_domains_in_ank_and_BD, enriched_clans_in_ank_and_BD = [], []
# for domain_name in enriched_domains_in_ank :
# 	if domain_name in enriched_domains_in_BD :
# 		enriched_domains_in_ank_and_BD.append(domain_name)
# for clan_name in enriched_clans_in_ank :
# 	if clan_name in enriched_clans_in_BD :
# 		enriched_clans_in_ank_and_BD.append(clan_name)

# print "len(enriched_domains_in_ank_and_BD)", len(enriched_domains_in_ank_and_BD) 
# print "len(enriched_clans_in_ank_and_BD)", len(enriched_clans_in_ank_and_BD)

# print enriched_domains_in_ank_and_BD
# print enriched_clans_in_ank_and_BD

# len(enriched_domains_in_ank_and_BD) 60
# len(enriched_clans_in_ank_and_BD) 60

# ['TRP_2', 'VPS9', 'Ion_trans', 'PARP', 'TIG', 'EGF', 'Ank_3', 'PID', 'R3H', 'FERM_M', 'cNMP_binding', 'DBB', 'BRCT', 'CAP_GLY', 'SH3_9', 'ArfGap', 'SH2', 'HECT', 'Pkinase_Tyr', 'SH3_1', 'PH', 'G-patch', 'SH3_2', 'SET', 'WGR', 'RHD', 'BTB', 'zf-RanBP', 'PARP_reg', 'Death', 'ZZ', 'Actin', 'zf-CCCH', 'DAGK_cat', 'PDZ', 'WH2', 'Pre-SET', 'LRR_7', 'Myosin_head', 'ZU5', 'Ank', 'Chromo', 'BAR', 'zf-C3HC4_3', 'C1_1', 'hEGF', 'Glutaminase', 'IBR', 'KilA-N', 'DUF3354', 'IQ', 'NODP', 'Notch', 'ELMO_CED12', 'EGF_CA', 'SAM_1', 'TPR_7', 'SAM_2', 'DUF3454', 'OTU']
# ['TRP_2', 'CL0033', 'CL0030', 'VPS9', 'CL0202', 'CL0125', 'CL0459', 'CL0145', 'R3H', 'CL0021', 'CL0020', 'CL0023', 'CL0022', 'FERM_M', 'cNMP_binding', 'DBB', 'CAP_GLY', 'CL0449', 'CL0159', 'ArfGap', 'CL0390', 'CL0229', 'CL0537', 'CL0306', 'CL0167', 'WGR', 'CL0093', 'CL0016', 'CL0541', 'CL0010', 'CL0013', 'PARP_reg', 'CL0073', 'CL0003', 'SET', 'CL0084', 'CL0001', 'CL0552', 'CL0006', 'CL0240', 'CL0108', 'WH2', 'Pre-SET', 'ZU5', 'CL0271', 'CL0384', 'IBR', 'CL0266', 'CL0186', 'DUF3454', 'DUF3354', 'CL0049', 'IQ', 'CL0041', 'NODP', 'Notch', 'ELMO_CED12', 'KilA-N', 'CL0465', 'CL0466']

######################################################################################################
############################        #based on enrichment > 2      ####################################
######################################################################################################


# enriched_domains_in_ank, enriched_clans_in_ank = [], []
# with open(domains_in_ank, 'rU') as file_open :
# 	data = file_open.readlines()
# 	for line in data[1:] :
# 		line = line.split("\t")
# 		domain_name = line[0]
# 		clan_name = line[1]
# 		enrichment = float(line[3])
# 		if enrichment > 2 :
# 			enriched_domains_in_ank.append(domain_name)
# 			enriched_clans_in_ank.append(clan_name)
# enriched_domains_in_ank = list(set(enriched_domains_in_ank))
# enriched_clans_in_ank = list(set(enriched_clans_in_ank))

# enriched_domains_in_BD, enriched_clans_in_BD = [], []
# with open(domains_in_BD, 'rU') as file_open :
# 	data = file_open.readlines()
# 	for line in data[1:] :
# 		line = line.split("\t")
# 		domain_name = line[0]
# 		clan_name = line[1]
# 		enrichment = float(line[3])
# 		if enrichment > 2 :
# 			enriched_domains_in_BD.append(domain_name)
# 			enriched_clans_in_BD.append(clan_name)
# enriched_domains_in_BD = list(set(enriched_domains_in_BD))
# enriched_clans_in_BD = list(set(enriched_clans_in_BD))


# enriched_domains_in_ank_and_BD, enriched_clans_in_ank_and_BD = [], []
# for domain_name in enriched_domains_in_ank :
# 	if domain_name in enriched_domains_in_BD :
# 		enriched_domains_in_ank_and_BD.append(domain_name)
# for clan_name in enriched_clans_in_ank :
# 	if clan_name in enriched_clans_in_BD :
# 		enriched_clans_in_ank_and_BD.append(clan_name)

# print "len(enriched_domains_in_ank_and_BD)", len(enriched_domains_in_ank_and_BD) 
# print "len(enriched_clans_in_ank_and_BD)", len(enriched_clans_in_ank_and_BD)

# print enriched_domains_in_ank_and_BD
# print enriched_clans_in_ank_and_BD

# len(enriched_domains_in_ank_and_BD) 40
# len(enriched_clans_in_ank_and_BD) 47
# ['TRP_2', 'Ion_trans', 'PARP', 'EGF', 'Ank_3', 'PID', 'FERM_M', 'DBB', 'CAP_GLY', 'SH3_9', 'ArfGap', 'HECT', 'Pkinase_Tyr', 'SH3_1', 'PH', 'SH3_2', 'WGR', 'RHD', 'zf-RanBP', 'PARP_reg', 'Death', 'ZZ', 'DAGK_cat', 'WH2', 'Pre-SET', 'ZU5', 'Chromo', 'BAR', 'zf-C3HC4_3', 'C1_1', 'hEGF', 'IBR', 'DUF3354', 'IQ', 'NODP', 'Notch', 'ELMO_CED12', 'SAM_1', 'SAM_2', 'DUF3454']
# ['TRP_2', 'CL0033', 'CL0030', 'CL0202', 'CL0125', 'CL0145', 'CL0021', 'CL0020', 'CL0023', 'CL0022', 'FERM_M', 'DBB', 'CAP_GLY', 'CL0159', 'ArfGap', 'CL0229', 'CL0537', 'CL0306', 'CL0167', 'WGR', 'CL0093', 'CL0016', 'CL0010', 'PARP_reg', 'CL0073', 'CL0003', 'CL0084', 'CL0001', 'CL0552', 'CL0006', 'CL0240', 'WH2', 'Pre-SET', 'ZU5', 'CL0271', 'IBR', 'CL0266', 'CL0186', 'DUF3454', 'DUF3354', 'CL0049', 'IQ', 'CL0041', 'NODP', 'Notch', 'ELMO_CED12', 'CL0465']


######################################################################################################
############################        #based on enrichment > 3      ####################################
######################################################################################################


# enriched_domains_in_ank, enriched_clans_in_ank = [], []
# with open(domains_in_ank, 'rU') as file_open :
# 	data = file_open.readlines()
# 	for line in data[1:] :
# 		line = line.split("\t")
# 		domain_name = line[0]
# 		clan_name = line[1]
# 		enrichment = float(line[3])
# 		if enrichment > 3 :
# 			enriched_domains_in_ank.append(domain_name)
# 			enriched_clans_in_ank.append(clan_name)
# enriched_domains_in_ank = list(set(enriched_domains_in_ank))
# enriched_clans_in_ank = list(set(enriched_clans_in_ank))

# enriched_domains_in_BD, enriched_clans_in_BD = [], []
# with open(domains_in_BD, 'rU') as file_open :
# 	data = file_open.readlines()
# 	for line in data[1:] :
# 		line = line.split("\t")
# 		domain_name = line[0]
# 		clan_name = line[1]
# 		enrichment = float(line[3])
# 		if enrichment > 3 :
# 			enriched_domains_in_BD.append(domain_name)
# 			enriched_clans_in_BD.append(clan_name)
# enriched_domains_in_BD = list(set(enriched_domains_in_BD))
# enriched_clans_in_BD = list(set(enriched_clans_in_BD))


# enriched_domains_in_ank_and_BD, enriched_clans_in_ank_and_BD = [], []
# for domain_name in enriched_domains_in_ank :
# 	if domain_name in enriched_domains_in_BD :
# 		enriched_domains_in_ank_and_BD.append(domain_name)
# for clan_name in enriched_clans_in_ank :
# 	if clan_name in enriched_clans_in_BD :
# 		enriched_clans_in_ank_and_BD.append(clan_name)

# print "len(enriched_domains_in_ank_and_BD)", len(enriched_domains_in_ank_and_BD) 
# print "len(enriched_clans_in_ank_and_BD)", len(enriched_clans_in_ank_and_BD)

# print enriched_domains_in_ank_and_BD
# print enriched_clans_in_ank_and_BD

# len(enriched_domains_in_ank_and_BD) 18
# len(enriched_clans_in_ank_and_BD) 27
# ['TRP_2', 'Ion_trans', 'EGF', 'Ank_3', 'PID', 'DBB', 'CAP_GLY', 'HECT', 'SH3_1', 'RHD', 'PARP_reg', 'Death', 'ZZ', 'ZU5', 'DUF3354', 'NODP', 'SAM_2', 'DUF3454']
# ['PARP_reg', 'DBB', 'TRP_2', 'CL0073', 'CL0159', 'ZU5', 'CL0020', 'CL0229', 'CL0202', 'CL0125', 'CL0001', 'CL0552', 'CL0006', 'CL0465', 'DUF3354', 'CL0030', 'CL0266', 'CL0145', 'CL0033', 'CL0167', 'CL0023', 'CL0041', 'CL0003', 'NODP', 'CL0010', 'CAP_GLY', 'DUF3454']


######################################################################################################
############################    based on Zscore enrichment > 1    ####################################
######################################################################################################
# # 
# enriched_domains_in_ank, enriched_clans_in_ank = [], []
# with open(domains_in_ank, 'rU') as file_open :
# 	data = file_open.readlines()
# 	for line in data[1:] :
# 		line = line.split("\t")
# 		domain_name = line[0]
# 		clan_name = line[1]
# 		Zscore = float(line[4])
# 		if Zscore > 1 :
# 			enriched_domains_in_ank.append(domain_name)
# 			enriched_clans_in_ank.append(clan_name)
# enriched_domains_in_ank = list(set(enriched_domains_in_ank))
# enriched_clans_in_ank = list(set(enriched_clans_in_ank))

# enriched_domains_in_BD, enriched_clans_in_BD = [], []
# with open(domains_in_BD, 'rU') as file_open :
# 	data = file_open.readlines()
# 	for line in data[1:] :
# 		line = line.split("\t")
# 		domain_name = line[0]
# 		clan_name = line[1]
# 		Zscore = float(line[4])
# 		if Zscore > 1 :
# 			enriched_domains_in_BD.append(domain_name)
# 			enriched_clans_in_BD.append(clan_name)
# enriched_domains_in_BD = list(set(enriched_domains_in_BD))
# enriched_clans_in_BD = list(set(enriched_clans_in_BD))


# enriched_domains_in_ank_and_BD, enriched_clans_in_ank_and_BD = [], []
# for domain_name in enriched_domains_in_ank :
# 	if domain_name in enriched_domains_in_BD :
# 		enriched_domains_in_ank_and_BD.append(domain_name)
# for clan_name in enriched_clans_in_ank :
# 	if clan_name in enriched_clans_in_BD :
# 		enriched_clans_in_ank_and_BD.append(clan_name)

# print "len(enriched_domains_in_ank_and_BD)", len(enriched_domains_in_ank_and_BD) 
# print "len(enriched_clans_in_ank_and_BD)", len(enriched_clans_in_ank_and_BD)

# print enriched_domains_in_ank_and_BD
# print enriched_clans_in_ank_and_BD

# len(enriched_domains_in_ank_and_BD) 2
# len(enriched_clans_in_ank_and_BD) 3
# ['DBB', 'TRP_2']
# ['DBB', 'TRP_2', 'CL0229']


######################################################################################################
############################    based on Zscore enrichment > 0.5  ####################################
######################################################################################################
# enriched_domains_in_ank, enriched_clans_in_ank = [], []
# with open(domains_in_ank, 'rU') as file_open :
# 	data = file_open.readlines()
# 	for line in data[1:] :
# 		line = line.split("\t")
# 		domain_name = line[0]
# 		clan_name = line[1]
# 		Zscore = float(line[4])
# 		if Zscore > 0.5 :
# 			enriched_domains_in_ank.append(domain_name)
# 			enriched_clans_in_ank.append(clan_name)
# enriched_domains_in_ank = list(set(enriched_domains_in_ank))
# enriched_clans_in_ank = list(set(enriched_clans_in_ank))

# enriched_domains_in_BD, enriched_clans_in_BD = [], []
# with open(domains_in_BD, 'rU') as file_open :
# 	data = file_open.readlines()
# 	for line in data[1:] :
# 		line = line.split("\t")
# 		domain_name = line[0]
# 		clan_name = line[1]
# 		Zscore = float(line[4])
# 		if Zscore > 0.5 :
# 			enriched_domains_in_BD.append(domain_name)
# 			enriched_clans_in_BD.append(clan_name)
# enriched_domains_in_BD = list(set(enriched_domains_in_BD))
# enriched_clans_in_BD = list(set(enriched_clans_in_BD))


# enriched_domains_in_ank_and_BD, enriched_clans_in_ank_and_BD = [], []
# for domain_name in enriched_domains_in_ank :
# 	if domain_name in enriched_domains_in_BD :
# 		enriched_domains_in_ank_and_BD.append(domain_name)
# for clan_name in enriched_clans_in_ank :
# 	if clan_name in enriched_clans_in_BD :
# 		enriched_clans_in_ank_and_BD.append(clan_name)

# print "len(enriched_domains_in_ank_and_BD)", len(enriched_domains_in_ank_and_BD) 
# print "len(enriched_clans_in_ank_and_BD)", len(enriched_clans_in_ank_and_BD)

# print enriched_domains_in_ank_and_BD
# print enriched_clans_in_ank_and_BD

# len(enriched_domains_in_ank_and_BD) 5
# len(enriched_clans_in_ank_and_BD) 10
# ['DBB', 'TRP_2', 'EGF', 'DUF3454', 'RHD']
# ['DBB', 'TRP_2', 'CL0073', 'CL0159', 'CL0229', 'CL0001', 'CL0023', 'CL0041', 'CL0003', 'DUF3454']







stop = timeit.default_timer()
print stop - start 