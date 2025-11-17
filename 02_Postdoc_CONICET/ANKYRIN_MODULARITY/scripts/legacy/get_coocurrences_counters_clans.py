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

elm_clan_mapping, clan_elm_mapping = {}, {}

with open("elm_interaction_clans_update.txt", 'rU') as file_open :
	data = file_open.readlines()
	for line in data[1:] :
		line = line.split("\t")
		elm_name = line[0]
		pfam_clan = line[-1].replace("\n","")
		
		if pfam_clan != "\N" :

			if elm_name not in elm_clan_mapping :
				elm_clan_mapping[elm_name] = [pfam_clan]
			else :
				elm_clan_mapping[elm_name].append(pfam_clan)

			if pfam_clan not in clan_elm_mapping :
				clan_elm_mapping[pfam_clan] = [elm_name]
			else :
				clan_elm_mapping[pfam_clan].append(elm_name)

print len(clan_elm_mapping), len(elm_clan_mapping)


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



# for i in clan_elm_mapping :
# 	print i, clan_elm_mapping[i]
# print "\n"
# for i in elm_clan_mapping :
# 	print i, elm_clan_mapping[i]


counter1, counter2, counter3, counter4 = 0,0,0,0
counter1bis, counter2bis, counter3bis, counter4bis = 0,0,0,0


with open("interacting_pairs_list.txt", 'rU') as file_open :
	data = file_open.readlines()
	for line in data :

		interacting_pair = line.replace("\n","")
		partnerA = line[0:6]
		partnerB = line[7:].replace("\n","")
		print interacting_pair

################################################################################################################################################
		clans_in_partnerA = []
		clans_in_partnerB = []
		elms_in_partnerA = []
		elms_in_partnerB = []
		#print partnerA

		if os.path.isfile("pfam_in_all_ank1234_and_BD/%s_pfam.txt" % partnerA) == True :
			clan_search_A = (getoutput("awk '{ print $15 }' pfam_in_all_ank1234_and_BD/%s_pfam.txt" % partnerA)).split("\n")

			#print clan_search_A

			for clan_in_A in clan_search_A :
				clans_in_partnerA.append(clan_in_A)
			
			clans_in_partnerA = list(set(clans_in_partnerA))
			#print clans_in_partnerA

		else :
			#print "%s n'a pas de clans pfam" % partnerA
			pass
		

		#print partnerB


		elms_in_partnerB = list(set((getoutput("grep \"%s\" elms_search_in_interacting_pairs.txt | awk '{ print $2 }'" % partnerB).split("\n"))))

		#print elms_in_partnerB


		for pfam_clan in clans_in_partnerA :
			#print pfam_clan
			counter1bis = counter1bis + 1
			if pfam_clan in clan_elm_mapping :
				binding_elms = clan_elm_mapping[pfam_clan]
				for elm in binding_elms :
					if elm in elms_in_partnerB :
						counter1 = counter1 + 1
					#	print elm
			else :
				#print "This Pfam clan %s was not shown to bind to any elm" % pfam_clan
				pass


	# print counter1
	# print counter1bis
	# print "\n"


############################################################################################################################################################################		
		clans_in_partnerA = []
		clans_in_partnerB = []
		elms_in_partnerA = []
		elms_in_partnerB = []
		#print partnerA

		elms_in_partnerA = list(set((getoutput("grep \"%s\" elms_search_in_interacting_pairs.txt | awk '{ print $2 }'" % partnerA).split("\n"))))

		

		#print partnerB

		if os.path.isfile("pfam_in_all_ank1234_and_BD/%s_pfam.txt" % partnerB) == True :
			clan_search_B = (getoutput("awk '{ print $15 }' pfam_in_all_ank1234_and_BD/%s_pfam.txt" % partnerB).split("\n"))

			#print clan_search_A

			for clan_in_B in clan_search_B :
				clans_in_partnerB.append(clan_in_B)
			
			clans_in_partnerB = list(set(clans_in_partnerB))
			#print clans_in_partnerA

		else :
			#print "%s n'a pas de clans pfam" % partnerB
			pass



		for elm in elms_in_partnerA :
			#print elm
			counter2bis = counter2bis + 1
			if elm in elm_clan_mapping :
				binding_clans = elm_clan_mapping[elm]
				for clan in binding_clans :
					if clan in clans_in_partnerB :
						counter2 = counter2 + 1
				#		print clan
			else :
				#print "This elm %s was not shown to bind to any clan" % elm
				pass

	# print counter2
	# print counter2bis
	# print "\n"
	################################################################################################################################################
		clans_in_partnerA = []
		clans_in_partnerB = []
		elms_in_partnerA = []
		elms_in_partnerB = []
		#print partnerB

		if os.path.isfile("pfam_in_all_ank1234_and_BD/%s_pfam.txt" % partnerB) == True :
			clan_search_B = (getoutput("awk '{ print $15 }' pfam_in_all_ank1234_and_BD/%s_pfam.txt" % partnerB).split("\n"))

			#print clan_search_B

			for clan_in_B in clan_search_B :
				clans_in_partnerB.append(clan_in_B)
			
			clans_in_partnerB = list(set(clans_in_partnerB))
			#print clans_in_partnerB

		else :
			#print "%s n'a pas de clans pfam" % partnerA
			pass
		

		#print partnerB


		elms_in_partnerA = list(set((getoutput("grep \"%s\" elms_search_in_interacting_pairs.txt | awk '{ print $2 }'" % partnerA).split("\n"))))

		#print elms_in_partnerA


		for pfam_clan in clans_in_partnerB :
			#print pfam_clan
			counter3bis = counter3bis + 1
			if pfam_clan in clan_elm_mapping :
				binding_elms = clan_elm_mapping[pfam_clan]
				for elm in binding_elms :
					if elm in elms_in_partnerA :
						counter3 = counter3+ 1
					#	print elm
			else :
				#print "This Pfam family %s was not shown to bind to any elm" % pfam_clan
				pass


	# print counter3
	# print counter3bis
	# print "\n"

#############################################################################################################################################################################		
		clans_in_partnerA = []
		clans_in_partnerB = []

		elms_in_partnerA = []
		elms_in_partnerB = []
		#print partnerB

		elms_in_partnerB = list(set((getoutput("grep \"%s\" elms_search_in_interacting_pairs.txt | awk '{ print $2 }'" % partnerB).split("\n"))))

		

		#print partnerA

		if os.path.isfile("pfam_in_all_ank1234_and_BD/%s_pfam.txt" % partnerA) == True :
			clan_search_A = (getoutput("awk '{ print $15 }' pfam_in_all_ank1234_and_BD/%s_pfam.txt" % partnerA).split("\n"))

			#print clan_search_A

			for clan_in_A in clan_search_A :
				clans_in_partnerA.append(clan_in_A)
			
			clans_in_partnerA = list(set(clans_in_partnerA))
			#print clans_in_partnerA

		else :
			#print "%s n'a pas de clans pfam" % partnerB
			pass



		for elm in elms_in_partnerB :
			#print elm
			counter4bis = counter4bis + 1
			if elm in elm_clan_mapping :
				binding_clans = elm_clan_mapping[elm]
				for clan in binding_clans :
					if clan in clans_in_partnerA :
						counter4 = counter4 + 1
						#print clan
			else :
				#print "This elm %s was not shown to bind to any clan" % elm
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