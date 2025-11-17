#!/usr/bin/python
# -*- coding: utf-8 -*-

import timeit
start = timeit.default_timer()
import sys
import os
import glob
import re

elm_list = []

top_elms_counts_in_ank_proteins, top_elms_counts_in_ankregions, top_elms_counts_in_non_ankregions, top_elms_counts_in_ank_repeats, top_elms_counts_in_ank_linkers = [], [], [], [], []

with open("top_elms_counts_in_ank_proteins_with_conservation_Zscores_colored.txt", 'rU') as file_open :
	data = file_open.readlines()
	for line in data :
		line = line.split("\t")
		elm_name = line[0]
		elm_list.append(elm_name)
		top_elms_counts_in_ank_proteins.append(elm_name)


with open("top_elms_counts_in_ankregions_with_conservation_Zscores_colored.txt", 'rU') as file_open :
	data = file_open.readlines()
	for line in data :
		line = line.split("\t")
		elm_name = line[0]
		elm_list.append(elm_name)
		top_elms_counts_in_ankregions.append(elm_name)

with open("top_elms_counts_in_non_ankregions_with_conservation_Zscores_colored.txt", 'rU') as file_open :
	data = file_open.readlines()
	for line in data :
		line = line.split("\t")
		elm_name = line[0]
		elm_list.append(elm_name)
		top_elms_counts_in_non_ankregions.append(elm_name)

with open("top_elms_counts_in_ank_repeats_with_conservation_Zscores_colored.txt", 'rU') as file_open :
	data = file_open.readlines()
	for line in data :
		line = line.split("\t")
		elm_name = line[0]
		elm_list.append(elm_name)
		top_elms_counts_in_ank_repeats.append(elm_name)


with open("top_elms_counts_in_ank_linkers_with_conservation_Zscores_colored.txt", 'rU') as file_open :
	data = file_open.readlines()
	for line in data :
		line = line.split("\t")
		elm_name = line[0]
		elm_list.append(elm_name)
		top_elms_counts_in_ank_linkers.append(elm_name)


elm_list = list(set(elm_list))

print len(elm_list)
for elm_name in elm_list :
	print elm_name
	if elm_name in top_elms_counts_in_ank_proteins :
		print "ank proteins"
	if elm_name in top_elms_counts_in_ankregions :
		print "ank regions"
	if elm_name in top_elms_counts_in_non_ankregions :
		print "non ank regions"
	if elm_name in top_elms_counts_in_ank_repeats :
		print "ank repeats"
	if elm_name in top_elms_counts_in_ank_linkers :
		print "ank linkers"
	print "\n"



stop = timeit.default_timer()
print stop - start 