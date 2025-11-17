#!/usr/bin/python
# -*- coding: utf-8 -*-

import timeit
start = timeit.default_timer()

import sys
import os
from collections import Counter
import re


family_clan_mapping_dict = {}

with open("Pfam-A.clans.txt", 'rU') as file_open :
	data = file_open.readlines()
	for line in data[1:] :
		line = line.split("\t")
		family = line[-2]
		clan = line[1]
		if clan != "\N" :
			family_clan_mapping_dict[family] = clan

#build dict of total swissprot counts for all PFAM domains
clan_count_dict = {}
with open("FrequenciesPfam_mapped.txt", 'rU') as file_read :
	data = file_read.readlines()
	for line in data[1:] :
		line = line.split("\t")
		domain = line[1]
		if domain in family_clan_mapping_dict :
			clan = family_clan_mapping_dict[domain]
		else :
			clan = domain

		domain_count_in_swissprot = int(line[-1].replace("\n",""))
		if clan not in clan_count_dict :
			clan_count_dict[clan] = domain_count_in_swissprot
		else :
			clan_count_dict[clan] = clan_count_dict[clan] + domain_count_in_swissprot





print len(clan_count_dict) #10'784

with open("FrequenciesPfam_mapped_CLANS.txt", 'a+') as file_write :
	for clan in clan_count_dict :
		file_write.write("%s\t%s\n" % (clan, clan_count_dict[clan])) 




stop = timeit.default_timer()
print stop - start