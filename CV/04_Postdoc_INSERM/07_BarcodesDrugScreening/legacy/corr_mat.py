#!/usr/bin/python
# -*- coding: utf-8 -*-
# jupyter notebook --browser='C:/Program Files (x86)/Google/Chrome/Application/chrome.exe %s'


import timeit
start = timeit.default_timer()

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
 

df = pd.read_csv("result6_fillna_control_renamed_filtered6.csv", sep=';', header=0, index_col=0)
print (df) #

matrix = df.corr(
    method = 'pearson',  # The method of correlation
    min_periods = 1).round(2)     # Min number of observations required)

print(matrix)

sns.heatmap(matrix, annot=True)
plt.show()


stop = timeit.default_timer()
print(stop - start)  