#!/usr/bin/python
# -*- coding: utf-8 -*-
# jupyter notebook --browser='C:/Program Files (x86)/Google/Chrome/Application/chrome.exe %s'


import timeit
start = timeit.default_timer()

import pandas as pd
import glob
import matplotlib.pyplot as plt

for fname in glob.glob("A:\\Downloads\\Projects\\workFromHome\\Projects\\drug_screening\\2022_Barcodes\\data\\combined_runs.csv"):
    print(fname)
    df = pd.read_csv(fname, sep=';', header=0, index_col=0)
    for col in df.columns:
    	print(col)


stop = timeit.default_timer()
print(stop - start)  