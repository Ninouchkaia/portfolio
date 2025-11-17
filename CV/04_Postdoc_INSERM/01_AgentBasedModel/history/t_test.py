#!/usr/bin/python
# -*- coding: utf-8 -*-
# jupyter notebook --browser='C:/Program Files (x86)/Google/Chrome/Application/chrome.exe %s'


import timeit
start = timeit.default_timer()

import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import ttest_ind

parameters_sets = pd.read_csv('patient_param_sets_9patients_param_renamed.csv', sep=';', index_col=0) 


parameters_sets = parameters_sets.T
# parameters_sets["class"] = parameters_sets["class"].astype("category")
print(parameters_sets)

group1 = parameters_sets[parameters_sets['class']==1.0]
group2 = parameters_sets[parameters_sets['class']==2.0]
# print(group1)

#perform independent two sample t-test
# for col in parameters_sets.columns :
# 	print("%s;%s" % (col, str(ttest_ind(group1[col], group2[col])).split("=")[2][:-1]))
	# print("%s;%s" % (col, str(ttest_ind(group1[col], group2[col], equal_var=False)).split("=")[2][:-1]))

# plt.show()


from scipy import stats
for col in parameters_sets.columns :
	print("%s;%s" % (col, str(stats.ks_2samp(group1[col], group2[col])).split("=")[2][:-1]))



stop = timeit.default_timer()
print(stop - start)  

