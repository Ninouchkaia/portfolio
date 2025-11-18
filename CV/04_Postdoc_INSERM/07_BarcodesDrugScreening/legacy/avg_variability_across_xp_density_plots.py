#!/usr/bin/python
# -*- coding: utf-8 -*-
# jupyter notebook --browser='C:/Program Files (x86)/Google/Chrome/Application/chrome.exe %s'

import timeit
start = timeit.default_timer()

import scipy.stats as stats 
import numpy as np
from statistics import mean
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


xp_list = ["exp300821","exp010821","exp040821","exp151121","exp271221","exp130921","exp200921",
"exp181021"] 

with open("barcodes_in_controls_stdev_4_sorted_barcode.csv") as file_read :
	data = file_read.readlines()
	avg_var_dict = {}
	for line in data[1:] :
		line = line.replace("\n","").split(";")
		experiment = line[0]
		barcode = line[1]
		avged_var = float(line[8])
		if barcode not in avg_var_dict :
			avg_var_dict[barcode] = [avged_var]
		else : 
			avg_var_dict[barcode].append(avged_var)

list_averages=[]
for barcode in avg_var_dict:
	list_averages.append(mean(avg_var_dict[barcode]))
	if len(avg_var_dict[barcode]) != 8 :
		print(avg_var_dict[barcode])

# with open("barcode_variability_within_controls_aveaged_across_xps.csv", "a+") as file_write :
# 	file_write.write("barcode\tavged_var\n")
# 	for barcode in avg_var_dict :
# 		file_write.write("%s\t%s\n" % (barcode,mean(avg_var_dict[barcode])))

# plt.hist(list_averages, bins=5000)

density = stats.kde.gaussian_kde(list_averages)
x = np.arange(min(list_averages), max(list_averages), 0.01)
plt.scatter(x, density(x))
plt.title('Density Plot for averaged variability (max-min)/avg \nof barcodes in controls over experiments')
plt.savefig("averaged_variability_over_xps.pdf")

stop = timeit.default_timer()
print(stop - start) 