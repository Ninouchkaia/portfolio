#!/usr/bin/python
# -*- coding: utf-8 -*-
# jupyter notebook --browser='C:/Program Files (x86)/Google/Chrome/Application/chrome.exe %s'

import timeit
start = timeit.default_timer()


import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import pylab as pl

 
df = pd.read_csv("barcodes_in_controls_stdev_4.csv", sep='\t', header=0, index_col=0)
print(df)
df = df.reset_index()

# axes = df['(max-min)/avg'].hist(by=df['experiment'],bins=100)
axes = df['(max-min)/avg'].hist(bins=100)

# plt.title('Frequency Plots for barcode variability (max-min)/avg \n across controls in each experiment\n(100 bins)')
plt.savefig('histo_plots_all_xps.pdf')

stop = timeit.default_timer()
print(stop - start) 



