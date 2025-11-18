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
 


with open("result6_fillna_control_renamed_filtered4.csv") as file_read:
	data = file_read.readlines()
	conditions, experiments, runs = [],[],[]
	for line in data[:1] :
		line = line.replace("\n","").split(";")
		for condition in line[1:] :
			condition = condition.split("_")
			# print(('\t').join(condition))


			conditions.append('_'.join(condition[0:2]))
			experiments.append(condition[2])
			runs.append(condition[3])
conditions = list(set(conditions))
experiments = list(set(experiments))
runs = list(set(runs))

# print("conditions")
# print(conditions)
# print("experiments")
# print(experiments)
# print("runs")
# print(runs)

# print("conditions")
# for i in sorted(conditions,key=str.lower) :
# 	print(i)
# print("experiments")
# for i in experiments :
# 	print(i)
# print("runs")
# for i in runs :
# 	print(i)

with open("barcodes_in_controls_stdev_4_from_fold_changes.csv", 'a+') as file_write:
	file_write.write("experiment\tbarcode\tcntrl1\tcntrl2\tcntrl3\tcntrl4\tavg\tstdev\t(max-min)/avg\n")
	with open("fold_change_result6_fillna_control_renamed_filtered4.csv") as file_read:
		data = file_read.readlines()
		for experiment in experiments :
			# print(experiment)
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
				controls_avg = round(statistics.mean(controls_reads), 2)
				controls_stdev = round(statistics.stdev(controls_reads), 2)
				# avg_std = round(controls_avg / controls_stdev, 2)
				if controls_avg != 0 :
					max_min = round(((max(controls_reads) - min(controls_reads)) / controls_avg), 2)
				else :
					max_min = '0'
				# print(experiment,line[0], '\t'.join([str(i) for i in controls_reads]), controls_avg, controls_stdev,avg_std, max_min)
				# print(experiment,line[0], '\t'.join([str(i) for i in controls_reads]), controls_avg, controls_stdev, max_min)
				# print(experiment, controls_avg, controls_stdev,avg_std, max_min)
				file_write.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (experiment,line[0], '\t'.join([str(i) for i in controls_reads]), controls_avg, controls_stdev, max_min))
stop = timeit.default_timer()
print(stop - start) 