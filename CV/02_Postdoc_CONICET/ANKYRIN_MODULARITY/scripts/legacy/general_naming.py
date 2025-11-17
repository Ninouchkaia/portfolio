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


filename = sys.argv[1]
ma_liste = []
# # general naming
with open(filename, 'rU') as file_open :
	data = file_open.readlines()
	for line in data[1:] :
		line = line.split("\t")
		name = line[0]
		ma_liste.append(name)
ma_liste = list(set(ma_liste))


filename2 = sys.argv[2]
ma_liste2 = []
# # general naming
with open(filename2, 'rU') as file_open :
	data = file_open.readlines()
	for line in data[1:] :
		line = line.split("\t")
		name = line[0]
		ma_liste2.append(name)
ma_liste2 = list(set(ma_liste2))

overlap = []
for name in ma_liste :
	if name in ma_liste2 :
		print name
		overlap.append(name)

print ", ".join(overlap)

stop = timeit.default_timer()
print stop - start 