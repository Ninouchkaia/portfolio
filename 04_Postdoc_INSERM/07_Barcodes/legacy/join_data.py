#!/usr/bin/python
# -*- coding: utf-8 -*-
# jupyter notebook --browser='C:/Program Files (x86)/Google/Chrome/Application/chrome.exe %s'


import timeit
start = timeit.default_timer()

import pandas as pd
import glob
path = "A:\\Downloads\\Projects\\workFromHome\\Projects\\drug_screening\\data\\*.csv"
df_list = []
run = 1
for fname in glob.glob("result6.csv"):
    print(fname)
    df = pd.read_csv(fname, sep=';', header=0, index_col=0)
    #print (df) #
    df_list.append(df)
    run = run + 1
    if run == 3 :
        break

result = pd.concat(df_list, axis=1)
result.to_csv('result7.csv', sep=';', index=True)










stop = timeit.default_timer()
print(stop - start)  