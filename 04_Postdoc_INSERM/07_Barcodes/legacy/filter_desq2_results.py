#!/usr/bin/python
# -*- coding: utf-8 -*-
# jupyter notebook --browser='C:/Program Files (x86)/Google/Chrome/Application/chrome.exe %s'


import timeit
start = timeit.default_timer()

import pandas as pd
import glob

# for fname in glob.glob("condition_*.csv"):
for fname in glob.glob("condition_*.csv"):
    print(fname)
    # df = pd.read_csv(fname, sep=',', header=0, index_col=0)
    df = pd.read_csv(fname, sep=';', header=0, index_col=0)
    filtered_df = df[df["pvalue"] < 0.05]
    # filtered_df = filtered_df[filtered_df["log2FoldChange"] > 1]
    filtered_df = filtered_df[(filtered_df["log2FoldChange"] < -1) | (filtered_df["log2FoldChange"] > 1)]
    filtered_df.to_csv('pval_log2fc_filtered_%s' % fname[14:], sep=';', index=True)
    # filtered_df.to_csv('pval_filtered\\filtered_%s' % fname, sep=';', index=True)
    # filtered_df.to_csv('pval_log2fc_filtered_pos\\%s' % fname[14:], sep=';', index=True)





stop = timeit.default_timer()
print(stop - start)  