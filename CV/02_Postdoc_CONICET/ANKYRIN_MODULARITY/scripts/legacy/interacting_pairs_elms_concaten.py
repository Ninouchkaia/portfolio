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


cmd = "mkdir elms_in_interacting_pairs"
os.system(cmd)

with open("interacting_pairs_list.txt", 'rU') as file_open :
	data = file_open.readlines()
	for line in data :
		interacting_pair = line.replace("\n","")
		partnerA = line[0:6]
		partnerB = line[7:].replace("\n","")

		elm_search_A = getoutput("grep \"%s\" elms_search_in_interacting_pairs.txt" % partnerA)
		elm_search_B = getoutput("grep \"%s\" elms_search_in_interacting_pairs.txt" % partnerB)

		with open("elms_in_interacting_pairs/elms_in_%s.txt" % interacting_pair , 'a+') as file_write :
			file_write.write("%s\n" % (elm_search_A))
			file_write.write("%s\n" % (elm_search_B))






stop = timeit.default_timer()
print stop - start 