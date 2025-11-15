#!/usr/bin/python
# -*- coding: utf-8 -*-
# jupyter notebook --browser='C:/Program Files (x86)/Google/Chrome/Application/chrome.exe %s'


import timeit
start = timeit.default_timer()

import pandas as pd
import glob
import matplotlib inline
from matplotlib_venn import venn2, venn2_circles
from matplotlib_venn import venn3, venn3_circles
from matplotlib import pyplot as plt


rownames_list = []
path = "A:\\Downloads\\Projects\\workFromHome\\Projects\\drug_screening\\data\\*.csv"
for fname in glob.glob("*.csv"):
# for fname in glob.glob("result6.csv"):
    print(fname)
    df = pd.read_csv(fname, sep=';', header=0, index_col=0)
    # get the row names
    row_names = df.index.values.tolist()
    rownames_list.append(row_names)

for row_names in rownames_list :

stop = timeit.default_timer()
print(stop - start)  