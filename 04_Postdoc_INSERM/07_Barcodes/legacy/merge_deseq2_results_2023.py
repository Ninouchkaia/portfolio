#!/usr/bin/python
# -*- coding: utf-8 -*-
# jupyter notebook --browser='C:/Program Files (x86)/Google/Chrome/Application/chrome.exe %s'


import timeit
start = timeit.default_timer()

import pandas as pd
import glob

colnames = []
# for fname in glob.glob("condition_*.csv"):
# for fname in glob.glob("pval_log2fc_filtered_*.csv"):
for fname in glob.glob("unfiltered_*.csv"):

    print(fname)
    colname = '_'.join(fname.split('_')[1:3])
    colnames.append(colname)
    print(colname)

merged_df = pd.DataFrame(columns = colnames)
# print(merged_df)


# # for fname in glob.glob("condition_*.csv"):
# for fname in glob.glob("pval_log2fc_filtered_*.csv"):
for fname in glob.glob("unfiltered_*.csv"):
    print(fname)  
    colname = '_'.join(fname.split('_')[1:3])
    df = pd.read_csv(fname, sep=';', header=0, index_col=0)
    logfc = df['log2FoldChange']
    merged_df[colname] = logfc

print(merged_df)

#df2 = merged_df.T.drop_duplicates().T

# merged_df = merged_df.dropna()

# print(merged_df)

#print(df2)


# df2.to_csv('merged_logfc_pval_filtered_deseq2_2023.csv', sep=';', index=True)
# merged_df.to_csv('merged_logfc_pval_filtered_deseq2_2023.csv', sep=';', index=True)
merged_df.to_csv('merged_unfiltered_deseq2_2023.csv', sep=';', index=True)





stop = timeit.default_timer()
print(stop - start)  