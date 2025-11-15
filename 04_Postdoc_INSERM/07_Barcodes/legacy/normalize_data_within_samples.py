#!/usr/bin/python
# -*- coding: utf-8 -*-
# jupyter notebook --browser='C:/Program Files (x86)/Google/Chrome/Application/chrome.exe %s'


import timeit
start = timeit.default_timer()

import pandas as pd
import glob
# for fname in glob.glob("*.csv"):
for fname in glob.glob("result6_fillna_control_renamed_filtered4.csv"):
    print(fname)
    df = pd.read_csv(fname, sep=';', header=0, index_col=0)
    print ("df") #
    print (df) #

    percents_df = df.apply(lambda x: 100 * x / float(x.sum()))
    print("percents_df")
    print(percents_df)

    percents_df.to_csv('result6_fillna_control_renamed_filtered4_percents.csv', sep=';', index=True)

    



stop = timeit.default_timer()
print(stop - start)  