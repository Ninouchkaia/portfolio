#!/usr/bin/python
# -*- coding: utf-8 -*-
# jupyter notebook --browser='C:/Program Files (x86)/Google/Chrome/Application/chrome.exe %s'

import timeit
start = timeit.default_timer()

from scipy.cluster.hierarchy import linkage, dendrogram, fcluster
from scipy.spatial.distance import squareform
import glob
import statistics

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
 

# experiments = ['exp200921','exp130921','exp300821','exp040821','exp181021','exp151121','exp271221','exp010821']
# experiments = ['exp200921']

### we use a list of experiments that go in the same order as in the original file for simplicity ie : 
experiments = ['exp010821', 'exp040821', 'exp300821', 'exp130921', 'exp200921', 'exp181021', 'exp151121','exp271221',]

# with open("fold_change_result6_fillna_control_renamed_filtered4.csv", 'a+') as file_write :
# with open("fold_change_result6_fillna_control_renamed_filtered4_averaged_over_repliactes.csv", 'a+') as file_write :


### we build as many files (fc) as there are experiments, and we will agreggate them afterwards
with open("result6_fillna_control_renamed_filtered4.csv") as file_read:
	data = file_read.readlines()
	conditions_list = []
	barcodes_list = []
	## get conditions list in the same order as in the original data
	for line in data[:1] :
		line = line.replace("\n","").split(";")
		for condition in line:
			conditions_list.append(condition)
	print(conditions_list)
	
	# for experiment in experiments : ## we are now looping the experiments, but not necessarily in the same order as in the original file
	for experiment in experiments : ## we are now looping the experiments, in the same order as in the original file
		print(experiment)
		with open("fold_change_result6_fillna_control_renamed_filtered4_%s.csv" % experiment, 'a+') as file_write :
			barcodes_dict={}
			controls_indexes = []
			conditions_indexes = []
			conditions_list = []
			### for each barcode, we want to get the number of reads in the controls of that experiment
			### we first keep in memory the positions of the corresponding columns
			for line in data[:1] :
				line = line.replace("\n","").split(";")
				for condition in line :
					if ("Contro" not in condition) and (experiment in condition) :
						conditions_indexes.append(int(line.index(condition)))
						conditions_list.append(condition)
					if ("Contro" in condition) and (experiment in condition) :
						# print(condition, line.index(condition))
						controls_indexes.append(int(line.index(condition)))
			print(controls_indexes)
			print(conditions_indexes, len(conditions_indexes))
			file_write.write(";%s\n" % (';'.join(conditions_list)))

			### then, we go over each barcode, and calculate the average of the number of reads within the controls
			for line in data[1:] :
				line = line.replace("\n","").split(";")
				barcode = line[0]
				controls_reads = []
				
				for index in controls_indexes :
					controls_reads.append(float(line[index]))
				controls_avg = statistics.mean(controls_reads)
				# print(controls_reads, controls_avg)
				
				for index in conditions_indexes :
					read = line[index]
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