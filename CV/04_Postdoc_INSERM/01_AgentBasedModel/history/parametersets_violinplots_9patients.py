#!/usr/bin/python
# -*- coding: utf-8 -*-
# jupyter notebook --browser='C:/Program Files (x86)/Google/Chrome/Application/chrome.exe %s'


import timeit
start = timeit.default_timer()

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

parameters_sets = pd.read_csv('patient_param_sets_9patients_param_renamed.csv', sep=';', index_col=0) 

# sns.violinplot(data=parameters_sets["parameter"])

# df = sns.load_dataset("titanic")
# print(df)
# sns.violinplot(data=df, x="deck", y="age", hue="alive", split=True)

parameters_sets = parameters_sets.T
parameters_sets["class"] = parameters_sets["class"].astype("int")
# parameters_sets["class"] = parameters_sets["class"].astype("category")

print(parameters_sets)



# sns.set(context="paper", palette="colorblind", style="ticks")
# g = sns.FacetGrid(parameters_sets, col=" alpha", sharey=False, size=4, aspect=.5)
# g = g.map(sns.violinplot, data=parameters_sets, y=" alpha", x="class",  bw=.15, hue="class", split=False, cut=0, scale="width")
# # g.fig.get_axes()[0].legend(False)

# sns.violinplot(data=parameters_sets, y=" alpha", x="class",  bw=.15, hue="class", split=True, cut=0, scale="width")

# fig, axes = plt.subplots(nrows=20)
# i=0
# for column in parameters_sets:
#     print(column)
#     sns.violinplot(data=parameters_sets, y=column, x="class",  bw=.15, hue="class", split=False, cut=0, scale="width",ax=axes[i])
#     i=i+1


# sns.violinplot(data=parameters_sets, y=" alpha", x="class",  bw=.15, split=False, cut=0, scale="width")


for column in parameters_sets:
	print(column)
	fig, ax = plt.subplots()
	# sns.violinplot(data=parameters_sets, y=column, x="class",  bw=.15, split=False, cut=0, scale="width") ## 1
	ax = sns.violinplot(data=parameters_sets, y=column, x="class",  bw=.30, split=False, cut=0, scale="width") ## 2
	ax.set_xticklabels(["Class A", "Class B"])
	ax.set(xlabel=None)
	sns.set(font_scale=1.5)
	# sns.violinplot(data=parameters_sets, y=column, x="class",  bw=.60, split=False, cut=0, scale="width") ## 3
	fig.savefig('%s_distrib.pdf' % column.replace(" ","") ) 
	fig.savefig('%s_distrib.png' % column.replace(" ","") ) 



# plt.show()

stop = timeit.default_timer()
print(stop - start)  

