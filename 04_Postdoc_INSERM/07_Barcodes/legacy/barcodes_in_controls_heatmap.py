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
barcodes_dict = {}
# with open("barcodes_in_controls_stdev_4.csv") as file_read:
with open("barcodes_in_controls_stdev_4_from_fold_changes.csv") as file_read:

	data = file_read.readlines()
	for experiment in experiments :
		for line in data[1:] :
			line = line.replace("\n","").split("\t")
			# print(line)
			xp_id = line[0]
			if xp_id == experiment :
				barcode = line[1]
				stdev = float(line[6])
				max_min = float(line[8])
				if barcode in barcodes_dict :
					barcodes_dict[barcode].append(max_min)
				else :
					barcodes_dict[barcode] = [max_min]

# print(barcodes_dict)

barcodes_df = pd.DataFrame.from_dict(barcodes_dict, orient='index',columns=experiments)
print(barcodes_df)
barcodes_df = barcodes_df.dropna()
print(barcodes_df)

filtered_barcodes_df = barcodes_df[barcodes_df < 2] 
filtered_barcodes_df = filtered_barcodes_df.dropna()
print(filtered_barcodes_df)

# sns.heatmap(barcodes_df, yticklabels=False)
# sns.heatmap(filtered_barcodes_df, yticklabels=False)

# sns.clustermap(barcodes_df,yticklabels=False)
sns.clustermap(filtered_barcodes_df,yticklabels=False)
plt.show()

# plt.savefig("clustermap_barcodes_in_controls_max_min.pdf", format='pdf')

		
stop = timeit.default_timer()
print(stop - start) 