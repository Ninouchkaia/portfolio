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


elm_clan_and_family_mapping, clan_and_family_elm_mapping = {}, {}
elm_names_mapped_to_family_only = []
with open("elm_interaction_clans_update.txt", 'rU') as file_open :
	data = file_open.readlines()
	for line in data[1:] :
		line = line.split("\t")
		elm_name = line[0]
		pfam_clan = line[-1].replace("\n","")
		pfam_family = line[1]


		if pfam_clan != "\N" :

			if elm_name not in elm_clan_and_family_mapping :
				elm_clan_and_family_mapping[elm_name] = [pfam_clan]
			else :
				elm_clan_and_family_mapping[elm_name].append(pfam_clan)

			if pfam_clan not in clan_and_family_elm_mapping :
				clan_and_family_elm_mapping[pfam_clan] = [elm_name]
			else :
				clan_and_family_elm_mapping[pfam_clan].append(elm_name)

		else :
			if elm_name not in elm_clan_and_family_mapping :
				elm_clan_and_family_mapping[elm_name] = [pfam_family]
			else :
				elm_clan_and_family_mapping[elm_name].append(pfam_family)

			if pfam_family not in clan_and_family_elm_mapping :
				clan_and_family_elm_mapping[pfam_family] = [elm_name]
			else :
				clan_and_family_elm_mapping[pfam_family].append(elm_name)

families_concerned = []
for elm_name in elm_clan_and_family_mapping :
	if len(elm_clan_and_family_mapping[elm_name]) == 1 :
		if elm_clan_and_family_mapping[elm_name][0][:2] == "PF" :
			print elm_name, elm_clan_and_family_mapping[elm_name][0]
			elm_names_mapped_to_family_only.append(elm_name)
			families_concerned.append(elm_clan_and_family_mapping[elm_name][0])
families_concerned = list(set(families_concerned))

families_concerned = list(set(families_concerned))
print len(clan_and_family_elm_mapping), len(elm_clan_and_family_mapping), len(elm_names_mapped_to_family_only), len(families_concerned)



# for i in clan_elm_mapping :
# 	print i, clan_elm_mapping[i]
# print "\n"
# for i in elm_clan_mapping :
# 	print i, elm_clan_mapping[i]



counter1, counter2 = 0,0
counter1bis, counter2bis = 0,0

ank_and_BD_list = []
with open("interacting_pairs_list_REORDED.txt", 'rU') as file_open :
	data = file_open.readlines()
	for line in data :
		partnerA = line[0:6]
		partnerB = line[7:].replace("\n","")
		ank_and_BD_list.append(partnerA)
		ank_and_BD_list.append(partnerB)

ank_and_BD_list = list(set(ank_and_BD_list))
print len(ank_and_BD_list) #2321

########################################################################        IN ANK AND BD    ######################################################################		

# for protein in ank_and_BD_list :
# 	print protein

# 	pfam_clan_search_protein = []
# 	domain_search_protein = []
# 	domains_in_protein = []
# 	pfam_structures_in_protein = []
# 	elms_in_protein = []

# 	if os.path.isfile("pfam_in_all_ank1234_and_BD/%s_pfam.txt" % protein) == True :
# 		pfam_clan_search_protein = (getoutput("awk '{ print $15 }' pfam_in_all_ank1234_and_BD/%s_pfam.txt" % protein)).split("\n")
# 		domain_search_protein = (getoutput("awk '{ print $6 }' pfam_in_all_ank1234_and_BD/%s_pfam.txt" % protein)).split("\n")
		
# 		for domain in domain_search_protein :
# 			domain_in_protein = (domain.split("."))[0]
# 			domains_in_protein.append(domain_in_protein)

# 		# print pfam_clan_search_protein, domains_in_protein
# 		# print len(pfam_clan_search_protein), len(domains_in_protein)

# 		if len(domains_in_protein) > 0 :
# 			for i in range(0,len(domains_in_protein)) :
# 				if pfam_clan_search_protein[i][:2] == 'CL' :
# 					pfam_structures_in_protein.append(pfam_clan_search_protein[i])
# 				else :
# 					pfam_structures_in_protein.append(domains_in_protein[i])
# 		else :
# 			pass		

# 		pfam_structures_in_protein = list(set(pfam_structures_in_protein))
# 		# print pfam_structures_in_protein


# 	else :
# 		#print "%s n'a pas de pfam_structures" % protein
# 		pass



# 	elms_in_protein = list(set((getoutput("grep \"%s\" elms_search_in_interacting_pairs.txt | awk '{ print $2 }'" % protein).split("\n"))))

# 	#print elms_in_partnerB


# 	for pfam_structure in pfam_structures_in_protein :
# 		#print pfam_structure
# 		counter1bis = counter1bis + 1
# 		if pfam_structure in clan_and_family_elm_mapping :
# 			binding_elms = clan_and_family_elm_mapping[pfam_structure]
# 			for elm in binding_elms :
# 				if elm in elms_in_protein :
# 					counter1 = counter1 + 1
# 				#	print elm
# 		else :
# 			#print "This Pfam structure %s was not shown to bind to any elm" % pfam_structure
# 			pass


# print counter1
# print counter1bis
# print "\n"




##########################################################      IN ANKYRINS ONLY        ##########################################################

# ank_list = []

# with open("mapping_table_unp_string_uniref50_UPPERCASE.txt", 'rU') as file_open :
# 	data = file_open.readlines()
# 	for line in data :
# 		line = line.split("\t")
# 		ank_id = line[1].replace("\n","")
# 		ank_list.append(ank_id)
# ank_list = list(set(ank_list))

# for protein in ank_list :
# 	print protein

# 	pfam_clan_search_protein = []
# 	domain_search_protein = []
# 	domains_in_protein = []
# 	pfam_structures_in_protein = []
# 	elms_in_protein = []

# 	if os.path.isfile("pfam_in_all_ank1234_and_BD/%s_pfam.txt" % protein) == True :
# 		pfam_clan_search_protein = (getoutput("awk '{ print $15 }' pfam_in_all_ank1234_and_BD/%s_pfam.txt" % protein)).split("\n")
# 		domain_search_protein = (getoutput("awk '{ print $6 }' pfam_in_all_ank1234_and_BD/%s_pfam.txt" % protein)).split("\n")
		
# 		for domain in domain_search_protein :
# 			domain_in_protein = (domain.split("."))[0]
# 			domains_in_protein.append(domain_in_protein)

# 		# print pfam_clan_search_protein, domains_in_protein
# 		# print len(pfam_clan_search_protein), len(domains_in_protein)

# 		if len(domains_in_protein) > 0 :
# 			for i in range(0,len(domains_in_protein)) :
# 				if pfam_clan_search_protein[i][:2] == 'CL' :
# 					pfam_structures_in_protein.append(pfam_clan_search_protein[i])
# 				else :
# 					pfam_structures_in_protein.append(domains_in_protein[i])
# 		else :
# 			pass		

# 		pfam_structures_in_protein = list(set(pfam_structures_in_protein))
# 		# print pfam_structures_in_protein


# 	else :
# 		#print "%s n'a pas de pfam_structures" % protein
# 		pass



# 	elms_in_protein = list(set((getoutput("grep \"%s\" elms_search_in_interacting_pairs.txt | awk '{ print $2 }'" % protein).split("\n"))))

# 	#print elms_in_partnerB


# 	for pfam_structure in pfam_structures_in_protein :
# 		#print pfam_structure
# 		counter1bis = counter1bis + 1
# 		if pfam_structure in clan_and_family_elm_mapping :
# 			binding_elms = clan_and_family_elm_mapping[pfam_structure]
# 			for elm in binding_elms :
# 				if elm in elms_in_protein :
# 					counter1 = counter1 + 1
# 				#	print elm
# 		else :
# 			#print "This Pfam structure %s was not shown to bind to any elm" % pfam_structure
# 			pass


# print counter1
# print counter1bis
# print "\n"


##########################################################      IN BD ONLY        ##########################################################

# ank_list = []
# BD_list = []

# with open("mapping_table_unp_string_uniref50_UPPERCASE.txt", 'rU') as file_open :
# 	data = file_open.readlines()
# 	for line in data :
# 		line = line.split("\t")
# 		ank_id = line[1].replace("\n","")
# 		ank_list.append(ank_id)
# ank_list = list(set(ank_list))

# BD_list = list(set(ank_and_BD_list)-set(ank_list))
# print len(BD_list)
# # 
# for protein in BD_list :
# 	print protein

# 	pfam_clan_search_protein = []
# 	domain_search_protein = []
# 	domains_in_protein = []
# 	pfam_structures_in_protein = []
# 	elms_in_protein = []

# 	if os.path.isfile("pfam_in_all_ank1234_and_BD/%s_pfam.txt" % protein) == True :
# 		pfam_clan_search_protein = (getoutput("awk '{ print $15 }' pfam_in_all_ank1234_and_BD/%s_pfam.txt" % protein)).split("\n")
# 		domain_search_protein = (getoutput("awk '{ print $6 }' pfam_in_all_ank1234_and_BD/%s_pfam.txt" % protein)).split("\n")
		
# 		for domain in domain_search_protein :
# 			domain_in_protein = (domain.split("."))[0]
# 			domains_in_protein.append(domain_in_protein)

# 		# print pfam_clan_search_protein, domains_in_protein
# 		# print len(pfam_clan_search_protein), len(domains_in_protein)

# 		if len(domains_in_protein) > 0 :
# 			for i in range(0,len(domains_in_protein)) :
# 				if pfam_clan_search_protein[i][:2] == 'CL' :
# 					pfam_structures_in_protein.append(pfam_clan_search_protein[i])
# 				else :
# 					pfam_structures_in_protein.append(domains_in_protein[i])
# 		else :
# 			pass		

# 		pfam_structures_in_protein = list(set(pfam_structures_in_protein))
# 		# print pfam_structures_in_protein


# 	else :
# 		#print "%s n'a pas de pfam_structures" % protein
# 		pass



# 	elms_in_protein = list(set((getoutput("grep \"%s\" elms_search_in_interacting_pairs.txt | awk '{ print $2 }'" % protein).split("\n"))))

# 	#print elms_in_partnerB


# 	for pfam_structure in pfam_structures_in_protein :
# 		#print pfam_structure
# 		counter1bis = counter1bis + 1
# 		if pfam_structure in clan_and_family_elm_mapping :
# 			binding_elms = clan_and_family_elm_mapping[pfam_structure]
# 			for elm in binding_elms :
# 				if elm in elms_in_protein :
# 					counter1 = counter1 + 1
# 				#	print elm
# 		else :
# 			#print "This Pfam structure %s was not shown to bind to any elm" % pfam_structure
# 			pass


# print counter1
# print counter1bis
# print "\n"


# ########################################################################        IN ANK AND BD    ######################################################################		
# for protein in ank_and_BD_list :
# 	print protein

# 	elms_in_protein = list(set((getoutput("grep \"%s\" elms_search_in_interacting_pairs.txt | awk '{ print $2 }'" % protein).split("\n"))))


# 	pfam_clan_search_protein = []
# 	domain_search_protein = []
# 	domains_in_protein = []
# 	pfam_structures_in_protein = []
# 	elms_in_protein = []

# 	if os.path.isfile("pfam_in_all_ank1234_and_BD/%s_pfam.txt" % protein) == True :
# 		pfam_clan_search_protein = (getoutput("awk '{ print $15 }' pfam_in_all_ank1234_and_BD/%s_pfam.txt" % protein)).split("\n")
# 		domain_search_protein = (getoutput("awk '{ print $6 }' pfam_in_all_ank1234_and_BD/%s_pfam.txt" % protein)).split("\n")
		
# 		for domain in domain_search_protein :
# 			domain_in_protein = (domain.split("."))[0]
# 			domains_in_protein.append(domain_in_protein)

# 		# print pfam_clan_search_protein, domains_in_protein
# 		# print len(pfam_clan_search_protein), len(domains_in_protein)

# 		if len(domains_in_protein) > 0 :
# 			for i in range(0,len(domains_in_protein)) :
# 				if pfam_clan_search_protein[i][:2] == 'CL' :
# 					pfam_structures_in_protein.append(pfam_clan_search_protein[i])
# 				else :
# 					pfam_structures_in_protein.append(domains_in_protein[i])
# 		else :
# 			pass		

# 		pfam_structures_in_protein = list(set(pfam_structures_in_protein))
# 		# print pfam_structures_in_protein


# 	else :
# 		#print "%s n'a pas de pfam_structures" % protein
# 		pass



# 	elms_in_protein = list(set((getoutput("grep \"%s\" elms_search_in_interacting_pairs.txt | awk '{ print $2 }'" % protein).split("\n"))))

# 	#print elms_in_partnerB


# 	for elm in elms_in_protein :
# 		#print elm
# 		counter2bis = counter2bis + 1
# 		if elm in elm_clan_and_family_mapping :
# 			binding_pfam_structures = elm_clan_and_family_mapping[elm]
# 			for pfam_structure in binding_pfam_structures :
# 				if pfam_structure in pfam_structures_in_protein :
# 					counter2 = counter2 + 1
# 			#		print clan
# 		else :
# 			#print "This elm %s was not shown to bind to any clan" % elm
# 			pass

# print counter2
# print counter2bis
# print "\n"

##########################################################      IN ANKYRINS ONLY        ##########################################################
# ank_list = []

# with open("mapping_table_unp_string_uniref50_UPPERCASE.txt", 'rU') as file_open :
# 	data = file_open.readlines()
# 	for line in data :
# 		line = line.split("\t")
# 		ank_id = line[1].replace("\n","")
# 		ank_list.append(ank_id)
# ank_list = list(set(ank_list))

# for protein in ank_list :
# 	print protein

# 	elms_in_protein = list(set((getoutput("grep \"%s\" elms_search_in_interacting_pairs.txt | awk '{ print $2 }'" % protein).split("\n"))))


# 	pfam_clan_search_protein = []
# 	domain_search_protein = []
# 	domains_in_protein = []
# 	pfam_structures_in_protein = []
# 	elms_in_protein = []

# 	if os.path.isfile("pfam_in_all_ank1234_and_BD/%s_pfam.txt" % protein) == True :
# 		pfam_clan_search_protein = (getoutput("awk '{ print $15 }' pfam_in_all_ank1234_and_BD/%s_pfam.txt" % protein)).split("\n")
# 		domain_search_protein = (getoutput("awk '{ print $6 }' pfam_in_all_ank1234_and_BD/%s_pfam.txt" % protein)).split("\n")
		
# 		for domain in domain_search_protein :
# 			domain_in_protein = (domain.split("."))[0]
# 			domains_in_protein.append(domain_in_protein)

# 		# print pfam_clan_search_protein, domains_in_protein
# 		# print len(pfam_clan_search_protein), len(domains_in_protein)

# 		if len(domains_in_protein) > 0 :
# 			for i in range(0,len(domains_in_protein)) :
# 				if pfam_clan_search_protein[i][:2] == 'CL' :
# 					pfam_structures_in_protein.append(pfam_clan_search_protein[i])
# 				else :
# 					pfam_structures_in_protein.append(domains_in_protein[i])
# 		else :
# 			pass		

# 		pfam_structures_in_protein = list(set(pfam_structures_in_protein))
# 		# print pfam_structures_in_protein


# 	else :
# 		#print "%s n'a pas de pfam_structures" % protein
# 		pass



# 	elms_in_protein = list(set((getoutput("grep \"%s\" elms_search_in_interacting_pairs.txt | awk '{ print $2 }'" % protein).split("\n"))))

# 	#print elms_in_partnerB


# 	for elm in elms_in_protein :
# 		#print elm
# 		counter2bis = counter2bis + 1
# 		if elm in elm_clan_and_family_mapping :
# 			binding_pfam_structures = elm_clan_and_family_mapping[elm]
# 			for pfam_structure in binding_pfam_structures :
# 				if pfam_structure in pfam_structures_in_protein :
# 					counter2 = counter2 + 1
# 			#		print clan
# 		else :
# 			#print "This elm %s was not shown to bind to any clan" % elm
# 			pass

# print counter2
# print counter2bis
# print "\n"

##########################################################      IN BD ONLY        ##########################################################
ank_list = []
BD_list = []

with open("mapping_table_unp_string_uniref50_UPPERCASE.txt", 'rU') as file_open :
	data = file_open.readlines()
	for line in data :
		line = line.split("\t")
		ank_id = line[1].replace("\n","")
		ank_list.append(ank_id)
ank_list = list(set(ank_list))

BD_list = list(set(ank_and_BD_list)-set(ank_list))
print len(BD_list)
# 
for protein in BD_list :
	print protein

	elms_in_protein = list(set((getoutput("grep \"%s\" elms_search_in_interacting_pairs.txt | awk '{ print $2 }'" % protein).split("\n"))))


	pfam_clan_search_protein = []
	domain_search_protein = []
	domains_in_protein = []
	pfam_structures_in_protein = []
	elms_in_protein = []

	if os.path.isfile("pfam_in_all_ank1234_and_BD/%s_pfam.txt" % protein) == True :
		pfam_clan_search_protein = (getoutput("awk '{ print $15 }' pfam_in_all_ank1234_and_BD/%s_pfam.txt" % protein)).split("\n")
		domain_search_protein = (getoutput("awk '{ print $6 }' pfam_in_all_ank1234_and_BD/%s_pfam.txt" % protein)).split("\n")
		
		for domain in domain_search_protein :
			domain_in_protein = (domain.split("."))[0]
			domains_in_protein.append(domain_in_protein)

		# print pfam_clan_search_protein, domains_in_protein
		# print len(pfam_clan_search_protein), len(domains_in_protein)

		if len(domains_in_protein) > 0 :
			for i in range(0,len(domains_in_protein)) :
				if pfam_clan_search_protein[i][:2] == 'CL' :
					pfam_structures_in_protein.append(pfam_clan_search_protein[i])
				else :
					pfam_structures_in_protein.append(domains_in_protein[i])
		else :
			pass		

		pfam_structures_in_protein = list(set(pfam_structures_in_protein))
		# print pfam_structures_in_protein


	else :
		#print "%s n'a pas de pfam_structures" % protein
		pass



	elms_in_protein = list(set((getoutput("grep \"%s\" elms_search_in_interacting_pairs.txt | awk '{ print $2 }'" % protein).split("\n"))))

	#print elms_in_partnerB


	for elm in elms_in_protein :
		#print elm
		counter2bis = counter2bis + 1
		if elm in elm_clan_and_family_mapping :
			binding_pfam_structures = elm_clan_and_family_mapping[elm]
			for pfam_structure in binding_pfam_structures :
				if pfam_structure in pfam_structures_in_protein :
					counter2 = counter2 + 1
			#		print clan
		else :
			#print "This elm %s was not shown to bind to any clan" % elm
			pass

print counter2
print counter2bis
print "\n"



stop = timeit.default_timer()
print stop - start 