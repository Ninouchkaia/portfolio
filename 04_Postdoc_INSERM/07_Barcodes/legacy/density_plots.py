#!/usr/bin/python
# -*- coding: utf-8 -*-
# jupyter notebook --browser='C:/Program Files (x86)/Google/Chrome/Application/chrome.exe %s'

import timeit
start = timeit.default_timer()


import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# tips = sns.load_dataset("tips")
# print(type(tips))
# print(tips)
 
df = pd.read_csv("barcodes_in_controls_stdev_4.csv", sep='\t', header=0, index_col=0)
print(df)

# df["(max-min)/avg"].plot.kde()
# plt.savefig("variability_per_barcode_within_controls_across_XPs.pdf")

# df.groupby('experiment')['(max-min)/avg'].plot.density(legend=True)
# plt.title('Density Plots for # of reads variability (max-min)/avg of barcodes in controls')
# plt.savefig("variability_per_barcode_within_controls_across_XPs_density_histograms.pdf")

# bw_methodstr, scalar or callable, optional
# The method used to calculate the estimator bandwidth. This can be ‘scott’, ‘silverman’, a scalar constant or a callable. 
# If None (default), ‘scott’ is used. 
# See scipy.stats.gaussian_kde for more information.
# https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.gaussian_kde.html#scipy.stats.gaussian_kde

#### each experiment #####
# ax = sns.violinplot(x = 'experiment', y = '(max-min)/avg', data = df.reset_index(),
# 	# scale="count", inner="stick",
# 	bw=0.05)
# ax.set_xticklabels(ax.get_xticklabels(),rotation = 45,Fontsize=8)
# ax.figure.tight_layout()
# ax.set_title('Violin Plots for barcode variability (max-min)/avg across controls\n(8 experiments, bw=0.05)', pad=5)

#### all experiments #####
ax = sns.violinplot(x = '(max-min)/avg', data = df.reset_index(),
	# scale="count", inner="stick",
	bw=0.05)
ax.set_xticklabels(ax.get_xticklabels(),rotation = 45,Fontsize=8)
ax.figure.tight_layout()
# ax.set_title('Violin Plots for barcode variability (max-min)/avg across controls\n(all 8 experiments included, bw=0.05)', pad=5)


# plt.show()
plt.savefig("violin_plot_barcode_variability_all_xp_bw0.05.pdf")


stop = timeit.default_timer()
print(stop - start) 