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


with open("barcode_variability_within_controls_aveaged_across_xps.csv") as file_read :
	data = file_read.readlines()
	avg_var_dict = {}
	for line in data[1:] :
		line = line.replace("\n","").split("\t")
		barcode = line[0]
		avged_var = float(line[1])
		if barcode in avg_var_dict :
			print("KO !!")
		avg_var_dict[barcode] = avged_var

list_avg_var = []
for barcode in avg_var_dict :
	list_avg_var.append(avg_var_dict[barcode])
plt.hist(list_avg_var, bins=500)
plt.title('Frequency Plot for barcode variability (max-min)/avg \n across controls averaged over experiments (500 bins)')
plt.savefig("freq_plot_averaged_variability_over_xps.pdf")
plt.savefig("freq_plot_averaged_variability_over_xps.png")


# density = stats.kde.gaussian_kde(list_avg_var)
# x = np.arange(min(list_avg_var), max(list_avg_var), 0.1)
# plt.scatter(x, density(x))
# plt.title('Density Plot for averaged variability (max-min)/avg \nof barcodes in controls over experiments')
# plt.savefig("averaged_variability_over_xps.pdf")
# plt.show()

stop = timeit.default_timer()
print(stop - start) 