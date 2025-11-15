#!/usr/bin/python
# -*- coding: utf-8 -*-
# jupyter notebook --browser='C:/Program Files (x86)/Google/Chrome/Application/chrome.exe %s'

import timeit
start = timeit.default_timer()

from scipy.cluster.hierarchy import linkage, dendrogram, fcluster
from scipy.spatial.distance import squareform
import glob
import statistics


import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
 

experiments = ['exp200921','exp130921','exp300821','exp040821','exp181021','exp151121','exp271221','exp010821']
with open("fold_change_result6_fillna_control_renamed_filtered4.csv", 'a+') as file_write :
	with open("result6_fillna_control_renamed_filtered4.csv") as file_read:
		data = file_read.readlines()
		file_write.write(data[0])
		for experiment in experiments :
			print(experiment)
			barcodes_dict={}
			controls_indexes = []
			for line in data[:1] :
				line = line.replace("\n","").split(";")
				for condition in line :
					if ("Contro" in condition) and (experiment in condition) :
						# print(condition, line.index(condition))
						controls_indexes.append(int(line.index(condition)))
			for line in data[1:] :
				line = line.replace("\n","").split(";")
				controls_reads = []
				for index in controls_indexes :
					controls_reads.append(float(line[index]))
				controls_avg = statistics.mean(controls_reads)
				# print(controls_reads, controls_avg)
				
				barcode = line[0]
				for read in line[1:] :
					if barcode not in barcodes_dict :
						if controls_avg != 0 :
							barcodes_dict[barcode] = [(float(read)/controls_avg)*100]
						else :
							barcodes_dict[barcode] = ['NaN']
					else :
						if controls_avg != 0 :
							barcodes_dict[barcode].append((float(read)/controls_avg)*100)
						else :
							barcodes_dict[barcode].append('NaN')
			# print(barcodes_dict)	
	for barcode in barcodes_dict :
		file_write.write("%s;%s\n" % (barcode,';'.join([str(i) for i in barcodes_dict[barcode]])))			


stop = timeit.default_timer()
print(stop - start) 