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

family_clan_mapping_dict = {}

with open("Pfam-A.clans.txt", 'rU') as file_open :
	data = file_open.readlines()
	for line in data[1:] :
		line = line.split("\t")
		family = line[3]
		clan = line[1]
		if clan != "\N" :
			family_clan_mapping_dict[family] = clan
# print family_clan_mapping_dict


filename = sys.argv[1]
with open("%s_with_clans.txt" % filename[:-4], 'a+') as file_write :
	with open(filename, 'rU') as file_open :
		data = file_open.readlines()
		first_line_splitted = data[0].split("\t")
		file_write.write("%s\tCLAN\t%s" % (first_line_splitted[0], '\t'.join(first_line_splitted[1:])))		
		for line in data[1:] :
			line = line.split("\t")
			domain_name = line[0]
			if domain_name in family_clan_mapping_dict :
				clan_name = family_clan_mapping_dict[domain_name]
			else :
				clan_name = "NULL"
			file_write.write("%s\t%s\t%s" % (domain_name, clan_name, '\t'.join(line[1:])))






stop = timeit.default_timer()
print stop - start 