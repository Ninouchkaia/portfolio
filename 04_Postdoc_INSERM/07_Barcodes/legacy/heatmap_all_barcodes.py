#!/usr/bin/python
# -*- coding: utf-8 -*-
# jupyter notebook --browser='C:/Program Files (x86)/Google/Chrome/Application/chrome.exe %s'

import timeit
start = timeit.default_timer()

import scipy.stats as stats 
import numpy as np
from statistics import mean
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# df = pd.read_csv("result6_fillna_control_renamed_filtered4_millions.csv", sep=';', header=0, index_col=0)



experiments = ['exp010821', 'exp040821', 'exp300821', 'exp130921', 'exp200921', 'exp181021', 'exp151121','exp271221',]

for exp in experiments :
    plt.figure()
    df = pd.read_csv("result6_fillna_control_renamed_filtered4_millions.csv", sep=';', header=0, index_col=0)

    # df.drop([col for col in df.columns if 'exp010821' not in col],axis=1,inplace=True)
    df.drop([col for col in df.columns if exp not in col],axis=1,inplace=True)
    df = np.log2(df)

    df.sort_values(by=list(df.columns), axis=0, ascending=True, inplace=False, kind='quicksort', na_position='last')
    print(df)

    xticks = df.columns

    # sns.color_palette("Blues", as_cmap=True)
    ax = plt.axes()

    res = sns.heatmap(df, annot=False,yticklabels=False, xticklabels=xticks, cmap="Blues", ax=ax)
    ax.set_title(exp)
    res.set_xticklabels(res.get_xticklabels(), fontsize = 4, rotation=90)
    res.figure.tight_layout()


    # plt.savefig("result6_fillna_control_renamed_filtered4_heatmap_exp010821.pdf")
    plt.savefig("result6_fillna_control_renamed_filtered4_millions_heatmap_log2_%s.png" % exp)



stop = timeit.default_timer()
print(stop - start) 
