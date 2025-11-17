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



# for i in clan_and_family_elm_mapping :
# 	print i, clan_and_family_elm_mapping[i]
# print "\n"
# for i in elm_clan_and_family_mapping :
# 	print i, elm_clan_and_family_mapping[i]



	


counter1, counter2, counter3, counter4 = 0,0,0,0
counter1bis, counter2bis, counter3bis, counter4bis = 0,0,0,0


with open("interacting_pairs_list.txt", 'rU') as file_open :
	data = file_open.readlines()
	for line in data[:10] :

		interacting_pair = line.replace("\n","")
		partnerA = line[0:6]
		partnerB = line[7:].replace("\n","")
		print interacting_pair

####################   SEARCH FOR PFAM STRUCTURE IN A | SEARCH FOR ELM IN B | PFAM IN A --> CORRESPONDING ELM IN B ? ##############
		pfam_clan_search_A = []
		domain_search_A = []
		domains_in_partnerA = []
		domains_in_partnerB = []
		pfam_structures_in_partnerA = []
		pfam_structures_in_partnerB = []
		elms_in_partnerA = []
		elms_in_partnerB = []

		if os.path.isfile("pfam_in_all_ank1234_and_BD/%s_pfam.txt" % partnerA) == True :
			pfam_clan_search_A = (getoutput("awk '{ print $15 }' pfam_in_all_ank1234_and_BD/%s_pfam.txt" % partnerA)).split("\n")
			domain_search_A = (getoutput("awk '{ print $6 }' pfam_in_all_ank1234_and_BD/%s_pfam.txt" % partnerA)).split("\n")
			
			for domain in domain_search_A :
				domain_in_A = (domain.split("."))[0]
				domains_in_partnerA.append(domain_in_A)

			# print pfam_clan_search_A, domains_in_partnerA
			# print len(pfam_clan_search_A), len(domains_in_partnerA)

			if len(domains_in_partnerA) > 0 :
				for i in range(0,len(domains_in_partnerA)) :
					if pfam_clan_search_A[i][:2] == 'CL' :
						pfam_structures_in_partnerA.append(pfam_clan_search_A[i])
					else :
						pfam_structures_in_partnerA.append(domains_in_partnerA[i])
			else :
				pass		

			pfam_structures_in_partnerA = list(set(pfam_structures_in_partnerA))
			# print pfam_structures_in_partnerA


		else :
			#print "%s n'a pas de pfam_structures" % partnerA
			pass
		



		elms_in_partnerB = list(set((getoutput("grep \"%s\" elms_search_in_interacting_pairs.txt | awk '{ print $2 }'" % partnerB).split("\n"))))

		print len(elms_in_partnerB), len(pfam_structures_in_partnerA)


		for pfam_structure in pfam_structures_in_partnerA :
			#print pfam_structure
			counter1bis = counter1bis + 1
			if pfam_structure in clan_and_family_elm_mapping :
				binding_elms = clan_and_family_elm_mapping[pfam_structure]
				for elm in binding_elms :
					if elm in elms_in_partnerB :
						counter1 = counter1 + 1
					#	print elm
			else :
				#print "This Pfam clan or family %s was not shown to bind to any elm" % pfam_structure
				pass


	# print counter1
	# print counter1bis
	# print "\n"


####################   SEARCH FOR ELM IN A | SEARCH FOR PFAM STRUCTURE IN B | ELM IN A --> CORRESPONDING PFAM IN B ? ##############

		pfam_clan_search_A = []
		domain_search_A = []
		domains_in_partnerA = []
		domains_in_partnerB = []
		pfam_structures_in_partnerA = []
		pfam_structures_in_partnerB = []
		elms_in_partnerA = []
		elms_in_partnerB = []

		elms_in_partnerA = list(set((getoutput("grep \"%s\" elms_search_in_interacting_pairs.txt | awk '{ print $2 }'" % partnerA).split("\n"))))



		


		if os.path.isfile("pfam_in_all_ank1234_and_BD/%s_pfam.txt" % partnerB) == True :
			pfam_clan_search_B = (getoutput("awk '{ print $15 }' pfam_in_all_ank1234_and_BD/%s_pfam.txt" % partnerB)).split("\n")
			domain_search_B = (getoutput("awk '{ print $6 }' pfam_in_all_ank1234_and_BD/%s_pfam.txt" % partnerB)).split("\n")
			
			for domain in domain_search_B :
				domain_in_B = (domain.split("."))[0]
				domains_in_partnerB.append(domain_in_B)

			# print pfam_clan_search_B, domains_in_partnerB
			# print len(pfam_clan_search_B), len(domains_in_partnerB)

			if len(domains_in_partnerB) > 0 :
				for i in range(0,len(domains_in_partnerB)) :
					if pfam_clan_search_B[i][:2] == 'CL' :
						pfam_structures_in_partnerB.append(pfam_clan_search_B[i])
					else :
						pfam_structures_in_partnerB.append(domains_in_partnerB[i])
			else :
				pass		

			pfam_structures_in_partnerB = list(set(pfam_structures_in_partnerB))
			# print pfam_structures_in_partnerB


		else :
			#print "%s n'a pas de pfam_structures" % partnerB
			pass

		print len(elms_in_partnerA), len(pfam_structures_in_partnerB)

		for elm in elms_in_partnerA :
			#print elm
			counter2bis = counter2bis + 1
			if elm in elm_clan_and_family_mapping :
				binding_pfam_structures = elm_clan_and_family_mapping[elm]
				for pfam_structure in binding_pfam_structures :
					if pfam_structure in pfam_structures_in_partnerB :
						counter2 = counter2 + 1
				#		print pfam_structure
			else :
				#print "This elm %s was not shown to bind to any pfam_structure" % elm
				pass

	# print counter2
	# print counter2bis
	# print "\n"


####################   SEARCH FOR PFAM STRUCTURE IN B | SEARCH FOR ELM IN A | PFAM IN B --> CORRESPONDING ELM IN A ? ##############

		pfam_clan_search_A = []
		domain_search_A = []
		domains_in_partnerA = []
		domains_in_partnerB = []
		pfam_structures_in_partnerA = []
		pfam_structures_in_partnerB = []
		elms_in_partnerA = []
		elms_in_partnerB = []

		if os.path.isfile("pfam_in_all_ank1234_and_BD/%s_pfam.txt" % partnerB) == True :
			pfam_clan_search_B = (getoutput("awk '{ print $15 }' pfam_in_all_ank1234_and_BD/%s_pfam.txt" % partnerB)).split("\n")
			domain_search_B = (getoutput("awk '{ print $6 }' pfam_in_all_ank1234_and_BD/%s_pfam.txt" % partnerB)).split("\n")
			
			for domain in domain_search_B :
				domain_in_B = (domain.split("."))[0]
				domains_in_partnerB.append(domain_in_B)

			# print pfam_clan_search_B, domains_in_partnerB
			# print len(pfam_clan_search_B), len(domains_in_partnerB)

			if len(domains_in_partnerB) > 0 :
				for i in range(0,len(domains_in_partnerB)) :
					if pfam_clan_search_B[i][:2] == 'CL' :
						pfam_structures_in_partnerB.append(pfam_clan_search_B[i])
					else :
						pfam_structures_in_partnerB.append(domains_in_partnerB[i])
			else :
				pass		

			pfam_structures_in_partnerB = list(set(pfam_structures_in_partnerB))
			# print pfam_structures_in_partnerB


		else :
			#print "%s n'a pas de pfam_structures" % partnerB
			pass
		

		elms_in_partnerA = list(set((getoutput("grep \"%s\" elms_search_in_interacting_pairs.txt | awk '{ print $2 }'" % partnerA).split("\n"))))

		#print elms_in_partnerA


		for pfam_structure in pfam_structures_in_partnerB :
			#print pfam_structure
			counter3bis = counter3bis + 1
			if pfam_structure in clan_and_family_elm_mapping :
				binding_elms = clan_and_family_elm_mapping[pfam_structure]
				for elm in binding_elms :
					if elm in elms_in_partnerA :
						counter3 = counter3+ 1
					#	print elm
			else :
				#print "This Pfam family %s was not shown to bind to any elm" % pfam_structure
				pass


	# print counter3
	# print counter3bis
	# print "\n"

####################   SEARCH FOR ELM IN B | SEARCH FOR PFAM STRUCTURE IN A | ELM IN B --> CORRESPONDING PFAM IN A ? ##############

		pfam_clan_search_A = []
		domain_search_A = []
		domains_in_partnerA = []
		domains_in_partnerB = []
		pfam_structures_in_partnerA = []
		pfam_structures_in_partnerB = []
		elms_in_partnerA = []
		elms_in_partnerB = []

		elms_in_partnerB = list(set((getoutput("grep \"%s\" elms_search_in_interacting_pairs.txt | awk '{ print $2 }'" % partnerB).split("\n"))))



		if os.path.isfile("pfam_in_all_ank1234_and_BD/%s_pfam.txt" % partnerA) == True :
			pfam_clan_search_A = (getoutput("awk '{ print $15 }' pfam_in_all_ank1234_and_BD/%s_pfam.txt" % partnerA)).split("\n")
			domain_search_A = (getoutput("awk '{ print $6 }' pfam_in_all_ank1234_and_BD/%s_pfam.txt" % partnerA)).split("\n")
			
			for domain in domain_search_A :
				domain_in_A = (domain.split("."))[0]
				domains_in_partnerA.append(domain_in_A)

			# print pfam_clan_search_A, domains_in_partnerA
			# print len(pfam_clan_search_A), len(domains_in_partnerA)

			if len(domains_in_partnerA) > 0 :
				for i in range(0,len(domains_in_partnerA)) :
					if pfam_clan_search_A[i][:2] == 'CL' :
						pfam_structures_in_partnerA.append(pfam_clan_search_A[i])
					else :
						pfam_structures_in_partnerA.append(domains_in_partnerA[i])
			else :
				pass		

			pfam_structures_in_partnerA = list(set(pfam_structures_in_partnerA))
			# print pfam_structures_in_partnerA


		else :
			#print "%s n'a pas de pfam_structures" % partnerA
			pass
		



		for elm in elms_in_partnerB :
			#print elm
			counter4bis = counter4bis + 1
			if elm in elm_clan_and_family_mapping :
				binding_pfam_structures = elm_clan_and_family_mapping[elm]
				for pfam_structure in binding_pfam_structures :
					if pfam_structure in pfam_structures_in_partnerA :
						counter4 = counter4 + 1
						#print pfam_structure
			else :
				#print "This elm %s was not shown to bind to any pfam_structure" % elm
				pass

	# print counter4
	# print counter4bis
	# print "\n"


print counter1
print counter1bis
print "\n"
print counter2
print counter2bis
print "\n"
print counter3
print counter3bis
print "\n"
print counter4
print counter4bis
print "\n"

stop = timeit.default_timer()
print stop - start 