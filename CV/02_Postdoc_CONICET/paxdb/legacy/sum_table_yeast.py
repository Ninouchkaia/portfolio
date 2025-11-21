#!/usr/bin/python
# -*- coding: utf-8 -*-

import timeit
start = timeit.default_timer()
import os
import sys
import re
from collections import Counter

# with open("yeast_table_with_decay.txt", 'rU') as decay_table :
# 	decay_data = decay_table.readlines()
# 	with open("yeast_new_table.txt", 'a+') as file_write :
# 		for line in decay_data :
# 			line = line.split("\t")
# 			file_write.write("%s\t%s\t%s\t%s\t%s\n" % (line[0], line[2], line[3], line[7].replace("\n",""), line[4]))

liste = []
with open("yeast_new_table.txt", 'rU') as table :
	data = table.readlines()
	for line in data[1:] :
		line = line.split("\t")
		liste.append(line[0])

print len(liste)
print liste

duplicates = [k for k,v in Counter(liste).items() if v>1]

print duplicates





stop = timeit.default_timer()
print stop - start 