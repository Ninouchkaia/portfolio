#!/usr/bin/python
# -*- coding: utf-8 -*-


# script will check for duplicates in tables


import timeit
start = timeit.default_timer()
import os
import sys
import re
from commands import getoutput # permet d'obtenir l'output d'une commande bash.


# the list of the 17 datasets to be analyzed
species_list = ["10090-M.musculus", "160490-S.pyogenes", "224308-Spectral_counting_B.subtili", "267671-L.interrogans", "3702-A.thaliana", "449447-M.aeruginosa", "4896-S.pombe", "4932-S.cerevisiae", "511145-E.coli", "593117-Spectral_counting_T.gammatolerans", "6239-C.elegans", "7227-D.melanogaster", "83332-M.tuberculosis", "9031-Spectral_counting_G.gallus", "9606-H.sapiens", "9913-B.taurus", "99287-Spectral_counting_S.typhimurium"]


for species in species_list:
	print species
	table = open("/home/nina/scripts/paxdb/abundance_datasets/tables/tables_2/%s_table2.txt" % species)
	table = table.readlines()
	mydict = {}
	for line in table[6:] : 
		line = line.split("\t")
		uniprot_id = line[0]
		string_id = line[1]
		if uniprot_id in mydict :
			mydict[uniprot_id].append(string_id)
		else :
			mydict[uniprot_id] = [string_id]
	#print mydict
		
	for key in mydict :
		if len(list(set(mydict[key]))) != len(mydict[key]) : # for each key, if the dictionary values dont contain any duplicates, len(list) and len(list(set(list))) should be the same.
				print "ERROR"













stop = timeit.default_timer()
print stop - start 