#!/usr/bin/python
# -*- coding: utf-8 -*-
# jupyter notebook --browser='C:/Program Files (x86)/Google/Chrome/Application/chrome.exe %s'


import timeit
start = timeit.default_timer()

import pandas as pd
import glob

with open("A:\\Downloads\\Projects\\workFromHome\\Projects\\drug_screening\\2022_Barcodes\\data\\combined_runs_filtered_avoid_zero_reads_scaled.csv") as file_read:
	with open("A:\\Downloads\\Projects\\workFromHome\\Projects\\drug_screening\\2022_Barcodes\\data\\true_barcode_list.txt", "a+") as file_write:
		file_write.write("lentibarcode_id\tlentibarcode_sequence\n")
		data = file_read.readlines()
		i = 1
		for line in data[1:] :
			line = line.replace("\n","").split(";")
			barcode = line[0]
			file_write.write("%s\t%s\n" % (i, barcode))
			i = i+1

stop = timeit.default_timer()
print(stop - start)   