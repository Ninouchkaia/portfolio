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

ank_list = []
with open("mapping_table_unp_string_uniref50_UPPERCASE.txt", 'rU') as file_open :
	data = file_open.readlines()
	for line in data :
		line = line.split("\t")
		ank_id = line[1].replace("\n", "")
		ank_list.append(ank_id)

print len(ank_list),  len(set(ank_list))



with open("interacting_pairs_list.txt", 'rU') as file_open :
	data = file_open.readlines()
	with open("interacting_pairs_list_REORDED.txt", 'a+') as file_write :
		for line in data :
			partnerA = line[0:6]
			partnerB = line[7:].replace("\n","")
			if partnerA in ank_list : #couvre le cas ou A=ank/B=bd et le cas ou A=ank/B=ank
				interacting_pair = "%s_%s" % (partnerA,partnerB)
			elif partnerA not in ank_list : #couvre les cas ou A=bd/B=ank
				interacting_pair = "%s_%s" % (partnerB,partnerA)
			else :
				print "wat"
			
			file_write.write("%s\n" % interacting_pair)				
