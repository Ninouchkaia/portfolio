#!/usr/bin/python
# -*- coding: utf-8 -*-

import timeit
start = timeit.default_timer()
import os
import sys
import re

decay_dict = {}
with open("abundance_decay_yeast.txt", 'rU') as yeast_decay :
	decay_data = yeast_decay.readlines()
	for line in decay_data :
		line = line.split("\t")
		string_id = line[1]
		print string_id
		decay = line[3].replace("\n", "")
		decay_dict[string_id] = decay
print decay_dict


with open("4932-S.cerevisiae_table.txt", 'rU') as yeast_table :
	yeast_data = yeast_table.readlines()
	with open("yeast_table_with_decay.txt", 'a+') as file_write :	
		file_write.write("%s\tcost_protein (5*n(aa))\tdecay\tprotein term (cost protein * decay)\n" % yeast_data[0].replace("\n", ""))	
		for line in yeast_data[1:] :
			line = line.replace("\n", "")
			line = line.split("\t")
			n_aa = float(line[1])
			cost_protein = 5 * n_aa
			string_id = line[0]
			if string_id in decay_dict :
				decay = decay_dict[string_id]
				protein_term = float(cost_protein) * float(decay)
				file_write.write("%s\t%s\t%s\t%s\n" % ("\t".join(line), cost_protein, decay, protein_term))


stop = timeit.default_timer()
print stop - start 