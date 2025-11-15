#!/usr/bin/python
# -*- coding: utf-8 -*-
# jupyter notebook --browser='C:/Program Files (x86)/Google/Chrome/Application/chrome.exe %s'

import timeit
start = timeit.default_timer()

from scipy.cluster.hierarchy import linkage, dendrogram, fcluster
from scipy.spatial.distance import squareform
import glob


import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
 
# df = pd.read_csv("result6_fillna_control_renamed_filtered6.csv", sep=';', header=0, index_col=0)

controls, time_zeros, compounds = [],[],[]

for fname in glob.glob("result6_fillna_control_renamed_filtered6.csv"):
    print(fname)
    
    df = pd.read_csv(fname, sep=';', header=0, index_col=0)
    print (df.columns) #

    for col in df.columns:
        print(col)
#     for col in df.columns :
#         if 'Contro' in col :
#             print (col)
#             controls.append(col)
#         if 'Temps' in col :
#         	time_zeros.append(col)
#         else :
#         	compounds.append(col)

# print(len(controls))
# print(len(time_zeros))
# print(len(compounds))
# print(set(controls).intersection(set(time_zeros)))







stop = timeit.default_timer()
print(stop - start) 