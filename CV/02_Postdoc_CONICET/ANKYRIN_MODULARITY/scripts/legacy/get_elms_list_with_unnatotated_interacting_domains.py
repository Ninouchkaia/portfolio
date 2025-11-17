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

elm_name_list = []
with open("elm_patterns_20140701.txt", 'rU') as file_open :
	data = file_open.readlines()
	for line in data :
		line = line.split("\t")
		elm_name = line[0]
		elm_name_list.append(elm_name)
print len(elm_name_list)


mapped_elms = []
with open("elm_interaction_domains_modified.txt", 'rU') as file_open :
	data = file_open.readlines()
	for line in data[1:] :
		line = line.split("\t")
		elm_name = line[0]
		mapped_elms.append(elm_name)
mapped_elms = list(set(mapped_elms))

print len(mapped_elms)


for elm_name in elm_name_list :
	if elm_name not in mapped_elms :
		print elm_name










stop = timeit.default_timer()
print stop - start 