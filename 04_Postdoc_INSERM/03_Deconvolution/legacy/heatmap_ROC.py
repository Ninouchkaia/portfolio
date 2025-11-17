#!/usr/bin/python
# -*- coding: utf-8 -*-

import timeit
start = timeit.default_timer()

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

df = pd.read_csv (
	'A:\\Downloads\\Projects\\workFromHome\\Projects\\deconvolution\\revision\\Figure 6\\scores_signatures_deconv_sorted-ROC AUC.csv', 
	sep=';',
	index_col=0)
print(df)

print(df.dtypes)

cols = df.columns
df[cols] = df[cols].apply(pd.to_numeric, errors='coerce')

print(df.dtypes)

df = df.sort_values(by=['mean lodo'], ascending=False)

mask = df.isnull()

plt.figure(figsize=(6, 8))

sns.heatmap(df, mask=mask,cmap="Greens",square=True)
plt.tight_layout()

# plt.show()

plt.savefig('heatmap_ROC_sorted_mean_lodo.pdf')  


stop = timeit.default_timer()
print(stop - start)  