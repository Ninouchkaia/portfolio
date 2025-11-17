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



top_elms_in_ank_proteins = []
top_elms_in_binding_partners = []
top_elms_in_interacting_pairs = []

top_families_in_ank_proteins = []
top_families_in_binding_partners = []
top_families_in_interacting_pairs = []

with open("top_elms_counts_in_ank_proteins_with_conservation_Zscores_colored.txt", 'rU') as file_open :
	data = file_open.readlines()
	for line in data[1:] :
		line = line.split("\t")
		name = line[0]
		top_elms_in_ank_proteins.append(name)

with open("top_elms_counts_in_binding_partners_with_conservation_Zscores_colored.txt", 'rU') as file_open :
	data = file_open.readlines()
	for line in data[1:] :
		line = line.split("\t")
		name = line[0]
		top_elms_in_binding_partners.append(name)

with open("top_elms_counts_in_interacting_pairs_Zscores_colored.txt", 'rU') as file_open :
	data = file_open.readlines()
	for line in data[1:] :
		line = line.split("\t")
		name = line[0]
		top_elms_in_interacting_pairs.append(name)

with open("top_Pfam_domains_in_all-ank-20130926.fasta_MaxHomologs_1000_Zscores_color.txt", 'rU') as file_open :
	data = file_open.readlines()
	for line in data[1:] :
		line = line.split("\t")
		name = line[0]
		top_families_in_ank_proteins.append(name)

with open("top_Pfam_domains_in_binding_partners_2038.fasta_MaxHomologs_1000_Zscores_color.txt", 'rU') as file_open :
	data = file_open.readlines()
	for line in data[1:] :
		line = line.split("\t")
		name = line[0]
		top_families_in_binding_partners.append(name)

with open("top_Pfam_domains_in_interacting_pairs_Zscores_colored_with_pfam_clan.txt", 'rU') as file_open :
	data = file_open.readlines()
	for line in data[1:] :
		line = line.split("\t")
		name = line[0]
		top_families_in_interacting_pairs.append(name)


elm_family_mapping, family_elm_mapping = {}, {}

with open("elm_interaction_domains_modified.txt", 'rU') as file_open :
	data = file_open.readlines()
	for line in data[1:] :
		line = line.split("\t")
		elm_name = line[0]
		pfam_family = line[1]
		
		if elm_name not in elm_family_mapping :
			elm_family_mapping[elm_name] = [pfam_family]
		else :
			elm_family_mapping[elm_name].append(pfam_family)

		if pfam_family not in family_elm_mapping :
			family_elm_mapping[pfam_family] = [elm_name]
		else :
			family_elm_mapping[pfam_family].append(elm_name)



counterC1, counterC2, counterC3, counterC4, counterC5, counterC6 = 0,0,0,0,0,0
counterC1bis, counterC2bis, counterC3bis, counterC4bis, counterC5bis, counterC6bis = 0,0,0,0,0,0

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

# ########################################################################        IN ANK AND BD    ######################################################################		

# for protein in ank_and_BD_list :
# 	print protein

# 	domains_in_protein = []
# 	elms_in_protein = []

# 	if os.path.isfile("pfam_in_all_ank1234_and_BD/%s_pfam.txt" % protein) == True :
# 		domain_search_protein = (getoutput("awk '{ print $6 }' pfam_in_all_ank1234_and_BD/%s_pfam.txt" % protein).split("\n"))

# 		#print domain_search_A

# 		for domain in domain_search_protein :
# 			domain_in_protein = (domain.split("."))[0]
# 			domains_in_protein.append(domain_in_protein)
		
# 		domains_in_protein = list(set(domains_in_protein))
# 		#print domains_in_partnerA

# 	else :
# 		#print "%s n'a pas de domains pfam" % partnerA
# 		pass



# 	elms_in_protein = list(set((getoutput("grep \"%s\" elms_search_in_interacting_pairs.txt | awk '{ print $2 }'" % protein).split("\n"))))

# 	#print elms_in_partnerB


# 	for pfam_family in domains_in_protein :
# 		#print pfam_family
# 		counter1bis = counter1bis + 1
# 		if pfam_family in family_elm_mapping :
# 			binding_elms = family_elm_mapping[pfam_family]
# 			for elm in binding_elms :
# 				if elm in elms_in_protein :
# 					counter1 = counter1 + 1
# 				#	print elm
# 		else :
# 			#print "This Pfam family %s was not shown to bind to any elm" % pfam_family
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

# 	domains_in_protein = []
# 	elms_in_protein = []

# 	if os.path.isfile("pfam_in_all_ank1234_and_BD/%s_pfam.txt" % protein) == True :
# 		domain_search_protein = (getoutput("awk '{ print $6 }' pfam_in_all_ank1234_and_BD/%s_pfam.txt" % protein).split("\n"))

# 		#print domain_search_A

# 		for domain in domain_search_protein :
# 			domain_in_protein = (domain.split("."))[0]
# 			domains_in_protein.append(domain_in_protein)
		
# 		domains_in_protein = list(set(domains_in_protein))
# 		#print domains_in_partnerA

# 	else :
# 		#print "%s n'a pas de domains pfam" % partnerA
# 		pass



# 	elms_in_protein = list(set((getoutput("grep \"%s\" elms_search_in_interacting_pairs.txt | awk '{ print $2 }'" % protein).split("\n"))))

# 	#print elms_in_partnerB


# 	for pfam_family in domains_in_protein :
# 		#print pfam_family
# 		counter1bis = counter1bis + 1
# 		if pfam_family in family_elm_mapping :
# 			binding_elms = family_elm_mapping[pfam_family]
# 			for elm in binding_elms :
# 				if elm in elms_in_protein :
# 					counter1 = counter1 + 1
# 				#	print elm
# 		else :
# 			#print "This Pfam family %s was not shown to bind to any elm" % pfam_family
# 			pass


# print counterC1
# print counterC1bis
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

# 	domains_in_protein = []
# 	elms_in_protein = []

# 	if os.path.isfile("pfam_in_all_ank1234_and_BD/%s_pfam.txt" % protein) == True :
# 		domain_search_protein = (getoutput("awk '{ print $6 }' pfam_in_all_ank1234_and_BD/%s_pfam.txt" % protein).split("\n"))

# 		#print domain_search_A

# 		for domain in domain_search_protein :
# 			domain_in_protein = (domain.split("."))[0]
# 			domains_in_protein.append(domain_in_protein)
		
# 		domains_in_protein = list(set(domains_in_protein))
# 		#print domains_in_partnerA

# 	else :
# 		#print "%s n'a pas de domains pfam" % partnerA
# 		pass



# 	elms_in_protein = list(set((getoutput("grep \"%s\" elms_search_in_interacting_pairs.txt | awk '{ print $2 }'" % protein).split("\n"))))

# 	#print elms_in_partnerB


# 	for pfam_family in domains_in_protein :
# 		#print pfam_family
# 		counter1bis = counter1bis + 1
# 		if pfam_family in family_elm_mapping :
# 			binding_elms = family_elm_mapping[pfam_family]
# 			for elm in binding_elms :
# 				if elm in elms_in_protein :
# 					counter1 = counter1 + 1
# 				#	print elm
# 		else :
# 			#print "This Pfam family %s was not shown to bind to any elm" % pfam_family
# 			pass


# print counter1
# print counter1bis
# print "\n"


# ########################################################################        IN ANK AND BD    ######################################################################		
# for protein in ank_and_BD_list :
# 	print protein

# 	domains_in_protein = []
# 	elms_in_protein = []

# 	elms_in_protein = list(set((getoutput("grep \"%s\" elms_search_in_interacting_pairs.txt | awk '{ print $2 }'" % protein).split("\n"))))

	

# 	#print partnerB

# 	if os.path.isfile("pfam_in_all_ank1234_and_BD/%s_pfam.txt" % protein) == True :
# 		domain_search_protein = (getoutput("awk '{ print $6 }' pfam_in_all_ank1234_and_BD/%s_pfam.txt" % protein).split("\n"))

# 		#print domain_search_protein

# 		for domain in domain_search_protein :
# 			domain_in_protein = (domain.split("."))[0]
# 			domains_in_protein.append(domain_in_protein)
		
# 		domains_in_protein = list(set(domains_in_protein))
# 		#print domains_in_protein

# 	else :
# 		#print "%s n'a pas de domains pfam" % partnerB
# 		pass



# 	for elm in elms_in_protein :
# 		#print elm
# 		counter2bis = counter2bis + 1
# 		if elm in elm_family_mapping :
# 			binding_domains = elm_family_mapping[elm]
# 			for domain in binding_domains :
# 				if domain in domains_in_protein :
# 					counter2 = counter2 + 1
# 			#		print domain
# 		else :
# 			#print "This elm %s was not shown to bind to any domain" % elm
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

# 	domains_in_protein = []
# 	elms_in_protein = []

# 	elms_in_protein = list(set((getoutput("grep \"%s\" elms_search_in_interacting_pairs.txt | awk '{ print $2 }'" % protein).split("\n"))))

	

# 	#print partnerB

# 	if os.path.isfile("pfam_in_all_ank1234_and_BD/%s_pfam.txt" % protein) == True :
# 		domain_search_protein = (getoutput("awk '{ print $6 }' pfam_in_all_ank1234_and_BD/%s_pfam.txt" % protein).split("\n"))

# 		#print domain_search_protein

# 		for domain in domain_search_protein :
# 			domain_in_protein = (domain.split("."))[0]
# 			domains_in_protein.append(domain_in_protein)
		
# 		domains_in_protein = list(set(domains_in_protein))
# 		#print domains_in_protein

# 	else :
# 		#print "%s n'a pas de domains pfam" % partnerB
# 		pass



# 	for elm in elms_in_protein :
# 		#print elm
# 		counter2bis = counter2bis + 1
# 		if elm in elm_family_mapping :
# 			binding_domains = elm_family_mapping[elm]
# 			for domain in binding_domains :
# 				if domain in domains_in_protein :
# 					counter2 = counter2 + 1
# 			#		print domain
# 		else :
# 			#print "This elm %s was not shown to bind to any domain" % elm
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

	domains_in_protein = []
	elms_in_protein = []

	elms_in_protein = list(set((getoutput("grep \"%s\" elms_search_in_interacting_pairs.txt | awk '{ print $2 }'" % protein).split("\n"))))

	

	#print partnerB

	if os.path.isfile("pfam_in_all_ank1234_and_BD/%s_pfam.txt" % protein) == True :
		domain_search_protein = (getoutput("awk '{ print $6 }' pfam_in_all_ank1234_and_BD/%s_pfam.txt" % protein).split("\n"))

		#print domain_search_protein

		for domain in domain_search_protein :
			domain_in_protein = (domain.split("."))[0]
			domains_in_protein.append(domain_in_protein)
		
		domains_in_protein = list(set(domains_in_protein))
		#print domains_in_protein

	else :
		#print "%s n'a pas de domains pfam" % partnerB
		pass



	for elm in elms_in_protein :
		#print elm
		counter2bis = counter2bis + 1
		if elm in elm_family_mapping :
			binding_domains = elm_family_mapping[elm]
			for domain in binding_domains :
				if domain in domains_in_protein :
					counter2 = counter2 + 1
			#		print domain
		else :
			#print "This elm %s was not shown to bind to any domain" % elm
			pass

print counter2
print counter2bis
print "\n"



stop = timeit.default_timer()
print stop - start 