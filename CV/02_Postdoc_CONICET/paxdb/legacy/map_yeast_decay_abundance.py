#!/usr/bin/python
# -*- coding: utf-8 -*-

import timeit
start = timeit.default_timer()
import os
import sys
import re
with open("SuppDataSet.txt", 'rU') as DataSet_file :
	DataSet = DataSet_file.readlines()
	with open("abundance_datasets/datasets/4932-S.cerevisiae_dataset.txt", 'rU') as yeast_abundanceFile :
		yeast_abundanceData = yeast_abundanceFile.readlines()
		with open('abundance_decay_yeast.txt', 'a+') as file_write :
			for line in DataSet[15:]:
				line = line.split("\t")
				if len(line)>1:
					yeast_id = line[0]
					decay = line[2]
					for line in yeast_abundanceData :
						line = line.replace("\n", "")	
						line = line.split("\t")
						string_id = line[1]
						if yeast_id == string_id[5:] :	
							file_write.write("%s\t%s\n" % ("\t".join(line), decay))
							print yeast_id


stop = timeit.default_timer()
print stop - start 