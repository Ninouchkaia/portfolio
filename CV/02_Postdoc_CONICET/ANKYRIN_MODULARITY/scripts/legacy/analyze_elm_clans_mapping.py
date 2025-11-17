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

elm_name_list, pfam_family_list, pfam_clan_list = [], [], []

with open("elm_interaction_clans.txt", 'rU') as file_open :
	data = file_open.readlines()
	for line in data[1:] :
		line = line.split("\t")
		elm_name = line[0]
		elm_name_list.append(elm_name)
		pfam_family = line[1]
		pfam_family_list.append(pfam_family)
		pfam_clan = line[-1].replace("\n", "")
		if pfam_clan != "\N" :
			pfam_clan_list.append(pfam_clan)


print len(elm_name_list), len(set(elm_name_list))
print len(pfam_family_list), len(set(pfam_family_list))
print len(pfam_clan_list), len(set(pfam_clan_list))






stop = timeit.default_timer()
print stop - start 