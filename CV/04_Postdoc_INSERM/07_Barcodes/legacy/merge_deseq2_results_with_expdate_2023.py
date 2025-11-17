#!/usr/bin/python
# -*- coding: utf-8 -*-
# jupyter notebook --browser='C:/Program Files (x86)/Google/Chrome/Application/chrome.exe %s'


import timeit
start = timeit.default_timer()

import pandas as pd
import glob

colnames = []
# for fname in glob.glob("condition_*.csv"):
for fname in glob.glob("pval_log2fc_filtered_*.csv"):
# for fname in glob.glob("unfiltered_*.csv"):

    # print(fname)
    colname = fname.split('_')
    # print(colname)
    # colname_sublist = [colname[1], colname[2], colname[5][0:-4]] # for unfiltered
    colname_sublist = [colname[3], colname[4], colname[7][0:-4]] # for filtered
    # print(colname_sublist)
    colname_joined = '_'.join(colname_sublist)
    colnames.append(colname_joined)
    print(colname_joined)

merged_df = pd.DataFrame(columns = colnames)
# print(merged_df)


# # for fname in glob.glob("condition_*.csv"):
for fname in glob.glob("pval_log2fc_filtered_*.csv"):
# for fname in glob.glob("unfiltered_*.csv"):
    # print(fname)  
    colname = fname.split('_')
    # colname_sublist = [colname[1], colname[2], colname[5][0:-4]] # for unfiltered
    colname_sublist = [colname[3], colname[4], colname[7][0:-4]] # for filtered
    colname_joined = '_'.join(colname_sublist)
    df = pd.read_csv(fname, sep=';', header=0, index_col=0)
    logfc = df['log2FoldChange']
    merged_df[colname_joined] = logfc

print(merged_df)

#df2 = merged_df.T.drop_duplicates().T

# merged_df = merged_df.dropna()

# print(merged_df)

#print(df2)


#merged_df.to_csv('merged_logfc_pval_filtered_deseq2_with_expdate_2023.csv', sep=';', index=True)
# merged_df.to_csv('merged_unfiltered_deseq2_with_expdate_2023.csv', sep=';', index=True)





stop = timeit.default_timer()
print(stop - start)  