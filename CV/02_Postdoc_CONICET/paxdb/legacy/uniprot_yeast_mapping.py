#!/usr/bin/python
# -*- coding: utf-8 -*-

import timeit
start = timeit.default_timer()
import os
import sys
import re


with open("SuppDataSet.txt", 'rU') as SuppDataSet :
	dataset = SuppDataSet.readlines()
	with open("yeast.txt", 'rU') as uniprot_file :
		uniprot_data = uniprot_file.readlines()
		with open("uniprot_yeast_mapping.txt", 'a+') as file_write:
			for line_dataset in dataset[12:] : 
				line_dataset = line_dataset.split("\t")
				if len(line_dataset)>1:
					yeast_ORF = line_dataset[0]
					for line_uniprot in uniprot_data[60:] :
						if yeast_ORF in line_uniprot :
							uniprot = line_uniprot[95:101]
							print yeast_ORF, uniprot
							file_write.write("%s\t%s\n" % (yeast_ORF, uniprot))



stop = timeit.default_timer()
print stop - start 