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

elms_in_ank = sys.argv[1]
elms_in_BD = sys.argv[2]




# #####################################################################################################
# ###########################        #based on enrichment > average_value     #########################
# #####################################################################################################


# # #average value for enrichment : 
# enrichment_values_in_ank = []
# with open(elms_in_ank, 'rU') as file_open :
# 	data = file_open.readlines()
# 	for line in data[1:] :
# 		line = line.split("\t")
# 		enrichment = float(line[6])
# 		enrichment_values_in_ank.append(enrichment)
# average_enrichment_in_ank = numpy.mean(enrichment_values_in_ank)
# print average_enrichment_in_ank #0.320082418061


# enrichment_values_in_BD = []
# with open(elms_in_BD, 'rU') as file_open :
# 	data = file_open.readlines()
# 	for line in data[1:] :
# 		line = line.split("\t")
# 		enrichment = float(line[6])
# 		enrichment_values_in_BD.append(enrichment)
# average_enrichment_in_BD = numpy.mean(enrichment_values_in_BD)
# print average_enrichment_in_BD #0.663061310126



# enriched_elms_in_ank = []
# with open(elms_in_ank, 'rU') as file_open :
# 	data = file_open.readlines()
# 	for line in data[1:] :
# 		line = line.split("\t")
# 		elm_name = line[0]
# 		enrichment = float(line[6])
# 		if enrichment > average_enrichment_in_ank :
# 			enriched_elms_in_ank.append(elm_name)
# enriched_elms_in_ank = list(set(enriched_elms_in_ank))

# enriched_elms_in_BD = []
# with open(elms_in_BD, 'rU') as file_open :
# 	data = file_open.readlines()
# 	for line in data[1:] :
# 		line = line.split("\t")
# 		elm_name = line[0]
# 		enrichment = float(line[6])
# 		if enrichment > average_enrichment_in_BD :
# 			enriched_elms_in_BD.append(elm_name)
# enriched_elms_in_BD = list(set(enriched_elms_in_BD))


# enriched_elms_in_ank_and_BD = []
# for elm_name in enriched_elms_in_ank :
# 	if elm_name in enriched_elms_in_BD :
# 		enriched_elms_in_ank_and_BD.append(elm_name)

# print "len(enriched_elms_in_ank_and_BD)", len(enriched_elms_in_ank_and_BD) 

# print enriched_elms_in_ank_and_BD

# len(enriched_elms_in_ank_and_BD) 46
# ['LIG_GLEBS_BUB3_1', 'LIG_Rb_LxCxE_1', 'TRG_ENDOCYTIC_2', 'DOC_PIKK_1', 'MOD_ASX_betaOH_EGF', 'LIG_NRBOX', 'LIG_CORNRBOX', 'DEG_MDM2_1', 'LIG_Rb_pABgroove_1', 'LIG_SH2_PTP2', 'LIG_SPRY_1', 'LIG_Actin_WH2_2', 'LIG_TYR_ITIM', 'LIG_Actin_WH2_1', 'LIG_AP2alpha_1', 'LIG_OCRL_FandH_1', 'LIG_PTB_Phospho_1', 'LIG_BIR_III_3', 'LIG_TYR_ITSM', 'LIG_Actin_RPEL_3', 'DOC_SPAK_OSR1_1', 'LIG_SH2_STAT5', 'DOC_AGCK_PIF_2', 'LIG_FAT_LD_1', 'TRG_PEX_2', 'LIG_PTB_Apo_2', 'LIG_Clathr_ClatBox_2', 'MOD_OFUCOSY', 'LIG_HOMEOBOX', 'TRG_AP2beta_CARGO_1', 'MOD_CMANNOS', 'LIG_EH1_1', 'TRG_NLS_MonoCore_2', 'LIG_TYR_ITAM', 'LIG_PCNA_PIPBox_1', 'DOC_AGCK_PIF_1', 'LIG_Sin3_3', 'LIG_SH2_GRB2', 'MOD_N-GLC_2', 'LIG_CtBP_PxDLS_1', 'MOD_OGLYCOS', 'LIG_eIF4E_1', 'TRG_PEX_1', 'LIG_BIR_III_1', 'LIG_CAP-Gly_1', 'MOD_CAAXbox']




# #####################################################################################################
# ###########################        #based on enrichment > average_value     #########################
# ###########################        and on conservation > average_value     #########################
# #####################################################################################################


# #average value for enrichment : 
enrichment_values_in_ank = []
with open(elms_in_ank, 'rU') as file_open :
	data = file_open.readlines()
	for line in data[1:] :
		line = line.split("\t")
		enrichment = float(line[6])
		enrichment_values_in_ank.append(enrichment)
average_enrichment_in_ank = numpy.mean(enrichment_values_in_ank)
print average_enrichment_in_ank #0.320082418061


enrichment_values_in_BD = []
with open(elms_in_BD, 'rU') as file_open :
	data = file_open.readlines()
	for line in data[1:] :
		line = line.split("\t")
		enrichment = float(line[6])
		enrichment_values_in_BD.append(enrichment)
average_enrichment_in_BD = numpy.mean(enrichment_values_in_BD)
print average_enrichment_in_BD #0.663061310126


# #average value for conservation : 
conservation_values_in_ank = []
with open(elms_in_ank, 'rU') as file_open :
	data = file_open.readlines()
	for line in data[1:] :
		line = line.split("\t")
		conservation = float(line[7])
		conservation_values_in_ank.append(conservation)
average_conservation_in_ank = numpy.mean(conservation_values_in_ank)
print average_conservation_in_ank #0.497151327904



conservation_values_in_BD = []
with open(elms_in_BD, 'rU') as file_open :
	data = file_open.readlines()
	for line in data[1:] :
		line = line.split("\t")
		conservation = float(line[7])
		conservation_values_in_BD.append(conservation)
average_conservation_in_BD = numpy.mean(conservation_values_in_BD)
print average_conservation_in_BD #0.40097531464



enriched_elms_in_ank = []
with open(elms_in_ank, 'rU') as file_open :
	data = file_open.readlines()
	for line in data[1:] :
		line = line.split("\t")
		elm_name = line[0]
		enrichment = float(line[6])
		conservation = float(line[7])
		if enrichment > average_enrichment_in_ank and conservation > average_conservation_in_ank :
			enriched_elms_in_ank.append(elm_name)
enriched_elms_in_ank = list(set(enriched_elms_in_ank))

enriched_elms_in_BD = []
with open(elms_in_BD, 'rU') as file_open :
	data = file_open.readlines()
	for line in data[1:] :
		line = line.split("\t")
		elm_name = line[0]
		enrichment = float(line[6])
		conservation = float(line[7])
		if enrichment > average_enrichment_in_BD and conservation > average_conservation_in_BD :
			enriched_elms_in_BD.append(elm_name)
enriched_elms_in_BD = list(set(enriched_elms_in_BD))


enriched_elms_in_ank_and_BD = []
for elm_name in enriched_elms_in_ank :
	if elm_name in enriched_elms_in_BD :
		enriched_elms_in_ank_and_BD.append(elm_name)

print "len(enriched_elms_in_ank_and_BD)", len(enriched_elms_in_ank_and_BD) 

print enriched_elms_in_ank_and_BD








######################################################################################################
############################        #based on enrichment > 1      ####################################
######################################################################################################


# enriched_elms_in_ank = []
# with open(elms_in_ank, 'rU') as file_open :
# 	data = file_open.readlines()
# 	for line in data[1:] :
# 		line = line.split("\t")
# 		elm_name = line[0]
# 		enrichment = float(line[6])
# 		if enrichment > 1 :
# 			enriched_elms_in_ank.append(elm_name)
# enriched_elms_in_ank = list(set(enriched_elms_in_ank))

# enriched_elms_in_BD = []
# with open(elms_in_BD, 'rU') as file_open :
# 	data = file_open.readlines()
# 	for line in data[1:] :
# 		line = line.split("\t")
# 		elm_name = line[0]
# 		enrichment = float(line[6])
# 		if enrichment > 1 :
# 			enriched_elms_in_BD.append(elm_name)
# enriched_elms_in_BD = list(set(enriched_elms_in_BD))


# enriched_elms_in_ank_and_BD = []
# for elm_name in enriched_elms_in_ank :
# 	if elm_name in enriched_elms_in_BD :
# 		enriched_elms_in_ank_and_BD.append(elm_name)

# print "len(enriched_elms_in_ank_and_BD)", len(enriched_elms_in_ank_and_BD) 

# print enriched_elms_in_ank_and_BD

# # en(enriched_elms_in_ank_and_BD) 28
# # ['LIG_TYR_ITAM', 'TRG_ENDOCYTIC_2', 'MOD_ASX_betaOH_EGF', 'LIG_eIF4E_1', 'LIG_CORNRBOX', 'DEG_MDM2_1', 'LIG_EH1_1', 'LIG_AP2alpha_1', 'LIG_OCRL_FandH_1', 'LIG_PTB_Phospho_1', 'LIG_GLEBS_BUB3_1', 'LIG_Actin_RPEL_3', 'LIG_SPRY_1', 'LIG_FAT_LD_1', 'TRG_PEX_2', 'LIG_PTB_Apo_2', 'LIG_Clathr_ClatBox_2', 'LIG_NRBOX', 'LIG_HOMEOBOX', 'LIG_PCNA_PIPBox_1', 'DOC_AGCK_PIF_1', 'DOC_AGCK_PIF_2', 'MOD_OGLYCOS', 'LIG_Rb_pABgroove_1', 'LIG_TYR_ITIM', 'MOD_OFUCOSY', 'LIG_CAP-Gly_1', 'MOD_CAAXbox']

######################################################################################################
############################        #based on enrichment > 2      ####################################
######################################################################################################
# enriched_elms_in_ank = []
# with open(elms_in_ank, 'rU') as file_open :
# 	data = file_open.readlines()
# 	for line in data[1:] :
# 		line = line.split("\t")
# 		elm_name = line[0]
# 		enrichment = float(line[6])
# 		if enrichment > 2 :
# 			enriched_elms_in_ank.append(elm_name)
# enriched_elms_in_ank = list(set(enriched_elms_in_ank))

# enriched_elms_in_BD = []
# with open(elms_in_BD, 'rU') as file_open :
# 	data = file_open.readlines()
# 	for line in data[1:] :
# 		line = line.split("\t")
# 		elm_name = line[0]
# 		enrichment = float(line[6])
# 		if enrichment > 2 :
# 			enriched_elms_in_BD.append(elm_name)
# enriched_elms_in_BD = list(set(enriched_elms_in_BD))


# enriched_elms_in_ank_and_BD = []
# for elm_name in enriched_elms_in_ank :
# 	if elm_name in enriched_elms_in_BD :
# 		enriched_elms_in_ank_and_BD.append(elm_name)

# print "len(enriched_elms_in_ank_and_BD)", len(enriched_elms_in_ank_and_BD) 

# print enriched_elms_in_ank_and_BD

# len(enriched_elms_in_ank_and_BD) 4
# ['MOD_ASX_betaOH_EGF', 'MOD_OGLYCOS', 'LIG_TYR_ITAM', 'LIG_PCNA_PIPBox_1']

######################################################################################################
############################        #based on enrichment > 3      ####################################
######################################################################################################

# enriched_elms_in_ank = []
# with open(elms_in_ank, 'rU') as file_open :
# 	data = file_open.readlines()
# 	for line in data[1:] :
# 		line = line.split("\t")
# 		elm_name = line[0]
# 		enrichment = float(line[6])
# 		if enrichment > 3 :
# 			enriched_elms_in_ank.append(elm_name)
# enriched_elms_in_ank = list(set(enriched_elms_in_ank))

# enriched_elms_in_BD = []
# with open(elms_in_BD, 'rU') as file_open :
# 	data = file_open.readlines()
# 	for line in data[1:] :
# 		line = line.split("\t")
# 		elm_name = line[0]
# 		enrichment = float(line[6])
# 		if enrichment > 3 :
# 			enriched_elms_in_BD.append(elm_name)
# enriched_elms_in_BD = list(set(enriched_elms_in_BD))


# enriched_elms_in_ank_and_BD = []
# for elm_name in enriched_elms_in_ank :
# 	if elm_name in enriched_elms_in_BD :
# 		enriched_elms_in_ank_and_BD.append(elm_name)

# print "len(enriched_elms_in_ank_and_BD)", len(enriched_elms_in_ank_and_BD) 

# print enriched_elms_in_ank_and_BD

# len(enriched_elms_in_ank_and_BD) 2
# ['MOD_OGLYCOS', 'MOD_ASX_betaOH_EGF']

######################################################################################################
############################    based on Zscore enrichment > 1    ####################################
######################################################################################################
# #
# enriched_elms_in_ank = []
# with open(elms_in_ank, 'rU') as file_open :
# 	data = file_open.readlines()
# 	for line in data[1:] :
# 		line = line.split("\t")
# 		elm_name = line[0]
# 		Zscore = float(line[8])
# 		if Zscore > 1 :
# 			enriched_elms_in_ank.append(elm_name)
# enriched_elms_in_ank = list(set(enriched_elms_in_ank))

# enriched_elms_in_BD = []
# with open(elms_in_BD, 'rU') as file_open :
# 	data = file_open.readlines()
# 	for line in data[1:] :
# 		line = line.split("\t")
# 		elm_name = line[0]
# 		Zscore = float(line[8])
# 		if Zscore > 1 :
# 			enriched_elms_in_BD.append(elm_name)
# enriched_elms_in_BD = list(set(enriched_elms_in_BD))


# enriched_elms_in_ank_and_BD = []
# for elm_name in enriched_elms_in_ank :
# 	if elm_name in enriched_elms_in_BD :
# 		enriched_elms_in_ank_and_BD.append(elm_name)

# print "len(enriched_elms_in_ank_and_BD)", len(enriched_elms_in_ank_and_BD) 

# print enriched_elms_in_ank_and_BD

# len(enriched_elms_in_ank_and_BD) 6
# ['MOD_ASX_betaOH_EGF', 'MOD_OGLYCOS', 'LIG_TYR_ITAM', 'LIG_CORNRBOX', 'LIG_PCNA_PIPBox_1', 'DOC_AGCK_PIF_1']


######################################################################################################
############################    based on Zscore enrichment > 0.5  ####################################
######################################################################################################

# enriched_elms_in_ank = []
# with open(elms_in_ank, 'rU') as file_open :
# 	data = file_open.readlines()
# 	for line in data[1:] :
# 		line = line.split("\t")
# 		elm_name = line[0]
# 		Zscore = float(line[8])
# 		if Zscore > 0.5 :
# 			enriched_elms_in_ank.append(elm_name)
# enriched_elms_in_ank = list(set(enriched_elms_in_ank))

# enriched_elms_in_BD = []
# with open(elms_in_BD, 'rU') as file_open :
# 	data = file_open.readlines()
# 	for line in data[1:] :
# 		line = line.split("\t")
# 		elm_name = line[0]
# 		Zscore = float(line[8])
# 		if Zscore > 0.5 :
# 			enriched_elms_in_BD.append(elm_name)
# enriched_elms_in_BD = list(set(enriched_elms_in_BD))


# enriched_elms_in_ank_and_BD = []
# for elm_name in enriched_elms_in_ank :
# 	if elm_name in enriched_elms_in_BD :
# 		enriched_elms_in_ank_and_BD.append(elm_name)

# print "len(enriched_elms_in_ank_and_BD)", len(enriched_elms_in_ank_and_BD) 

# print enriched_elms_in_ank_and_BD

# len(enriched_elms_in_ank_and_BD) 18
# ['LIG_TYR_ITAM', 'MOD_ASX_betaOH_EGF', 'LIG_CORNRBOX', 'DEG_MDM2_1', 'LIG_EH1_1', 'LIG_AP2alpha_1', 'LIG_Actin_RPEL_3', 'LIG_SPRY_1', 'LIG_FAT_LD_1', 'TRG_PEX_2', 'LIG_HOMEOBOX', 'LIG_PCNA_PIPBox_1', 'DOC_AGCK_PIF_1', 'DOC_AGCK_PIF_2', 'MOD_OGLYCOS', 'LIG_Rb_pABgroove_1', 'MOD_OFUCOSY', 'MOD_CAAXbox']

stop = timeit.default_timer()
print stop - start 