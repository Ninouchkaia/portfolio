#!/usr/bin/python
# -*- coding: utf-8 -*-
# jupyter notebook --browser='C:/Program Files (x86)/Google/Chrome/Application/chrome.exe %s'


import timeit
start = timeit.default_timer()

import pandas as pd
import glob

print('HELLO')
for fname in glob.glob("/home/nina/Projects/2022_Barcodes/NetBioMed2022/data/combined_runs.csv"):
    print(fname)
    df = pd.read_csv(fname, sep=';', header=0, index_col=0)
    print(df)
    df = df.fillna(0).apply(lambda x: 1000000 * x / float(x.sum()))

    df.to_csv('/home/nina/Projects/2022_Barcodes/NetBioMed2022/data/combined_runs_normalized.csv', sep=';', index=True, float_format='%.2f')

stop = timeit.default_timer()
print(stop - start)   