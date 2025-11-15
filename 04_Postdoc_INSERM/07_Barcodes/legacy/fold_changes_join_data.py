#!/usr/bin/python
# -*- coding: utf-8 -*-
# jupyter notebook --browser='C:/Program Files (x86)/Google/Chrome/Application/chrome.exe %s'


import timeit
start = timeit.default_timer()

import pandas as pd
import glob
df_list = []

for fname in glob.glob("fold_change_result6_fillna_control_renamed_filtered4_*.csv"):
    print(fname)
    df = pd.read_csv(fname, sep=';', header=0, index_col=0)
    print (df) #
    df_list.append(df)

print(len(df_list))

# for i in range(0,len(df_list)-1) :
# 	print(i) 
# 	result = pd.concat([df_list[i], df_list[i+1]], axis=1)

result = pd.concat(df_list, axis=1)

result.to_csv('fold_changes_agreggated.csv', sep=';', index=True)










stop = timeit.default_timer()
print(stop - start)  