#!/usr/bin/python
# -*- coding: utf-8 -*-
# jupyter notebook --browser='C:/Program Files (x86)/Google/Chrome/Application/chrome.exe %s'


import timeit
start = timeit.default_timer()

import pandas as pd
import glob
# path = "A:\\Downloads\\Projects\\workFromHome\\Projects\\drug_screening\\data\\joined\\*.csv"
# for fname in glob.glob("*.csv"):
for fname in glob.glob("joined\\result7.csv"):
    print(fname)
    df = pd.read_csv(fname, sep=';', header=0, index_col=0)
    # print (df)
    df_filled = df.fillna(0)
    df_filled.to_csv('result7_fillna.csv', sep=';', index=True)






stop = timeit.default_timer()
print(stop - start)  