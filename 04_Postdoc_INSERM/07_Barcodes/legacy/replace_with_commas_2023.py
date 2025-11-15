#!/usr/bin/python
# -*- coding: utf-8 -*-
# jupyter notebook --browser='C:/Program Files (x86)/Google/Chrome/Application/chrome.exe %s'


import timeit
start = timeit.default_timer()

import pandas as pd
import glob

with open("A:\\Downloads\\Projects\\workFromHome\\Projects\\drug_screening\\2022_Barcodes\\data\\combined_runs_filtered_avoid_zero_reads_scaled.csv") as file_read:
	with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\drug_screening\\2022_Barcodes\\data\\combined_runs_filtered_avoid_zero_reads_scaled_commas.csv', 'a+') as file_write:
		data = file_read.readlines()
		file_write.write(data[0])
		for line in data[1:] :
			line = line.replace(".", ",")
			file_write.write(line)

stop = timeit.default_timer()
print(stop - start)   