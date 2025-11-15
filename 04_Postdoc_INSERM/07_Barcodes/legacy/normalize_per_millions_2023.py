#!/usr/bin/python
# -*- coding: utf-8 -*-
# jupyter notebook --browser='C:/Program Files (x86)/Google/Chrome/Application/chrome.exe %s'


import timeit
start = timeit.default_timer()

import pandas as pd
import glob

with open("A:\\Downloads\\Projects\\workFromHome\\Projects\\drug_screening\\2022_Barcodes\\data\\combined_runs_filtered_avoid_zero_reads.csv") as file_read:
    df = pd.read_csv(file_read, sep=';', header=0, index_col=0)
    print(df)
    df = df.fillna(0).apply(lambda x: 1000000 * x / float(x.sum()))

    df.to_csv('A:\\Downloads\\Projects\\workFromHome\\Projects\\drug_screening\\2022_Barcodes\\data\\combined_runs_filtered_avoid_zero_reads_scaled.csv', sep=';', index=True, float_format='%.6f')

stop = timeit.default_timer()
print(stop - start)   