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

# ax = sns.swarmplot(x="barcode", y="(max-min)/avg", hue="experiment", data=df)
# ax = sns.swarmplot(x=df["(max-min)/avg"])
# plt.show()

sns.set_theme(style="whitegrid")
ax = sns.swarmplot(x=df["(max-min)/avg"])
plt.savefig("test.pdf")

stop = timeit.default_timer()
print(stop - start) 