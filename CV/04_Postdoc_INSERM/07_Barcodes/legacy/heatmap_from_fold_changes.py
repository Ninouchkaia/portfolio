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


df = pd.read_csv("fold_changes_agreggated.csv", sep=';', header=0, index_col=0)
print(df)


### average the fold changes over the 4 replicates, for each condition...

sns.heatmap(df, annot=True)
plt.savefig("fold_changes_agreggated_heatmap.pdf")


stop = timeit.default_timer()
print(stop - start) 