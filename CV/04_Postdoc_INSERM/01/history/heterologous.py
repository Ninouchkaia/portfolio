#!/usr/bin/python
# -*- coding: utf-8 -*-
# jupyter notebook --browser='C:/Program Files (x86)/Google/Chrome/Application/chrome.exe %s'


import timeit
start = timeit.default_timer()

import pandas as pd
import seaborn as sns
import io
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure

# sns.set(rc={'figure.figsize':(10,10)})
# sns.set(style="ticks")
# sns.axes_style("whitegrid")


# df = pd.read_csv('allPatients-heterologous.tsv', sep='\t')
# print(df)

# sns.catplot(x='mono', y='viability', hue='patient', data=df, kind='bar')
# plt.grid()  #just add this

# plt.show()


df_via_exp = pd.read_csv('allPatients_prediction-via-exp.tsv', sep='\t')
df_via_exp['viability_mean'] = df_via_exp.mean(axis=1)
df_via_exp['viability_std'] = df_via_exp.std(axis=1)
df_via_exp = df_via_exp.iloc[::-1]

df_conc_exp = pd.read_csv('allPatients_prediction-conc-exp.tsv', sep='\t')
df_conc_exp['concentration_mean'] = df_conc_exp.mean(axis=1)
df_conc_exp['concentration_std'] = df_conc_exp.std(axis=1)
df_conc_exp = df_conc_exp.iloc[::-1]

fig, axs = plt.subplots(2,1)
fig.suptitle('B-CLL cells survival at Day 9 in heterologous co-cultures\nwith varying initial monocytes proportions', fontsize=14, y=1.02)

df_via_exp.plot(kind = "bar", y = "viability_mean", legend = False, 
	yerr = "viability_std", rot=0, color = "black",
	capsize=3, ax = axs[0])

axs[0].set_ylabel('Viability (%)')
# axs[0].grid(True)
# axs[0].set_axisbelow(True)

df_conc_exp.plot(kind = "bar", y = "concentration_mean", legend = False, 
	yerr = "concentration_std", rot=0, color = "black",
	capsize=3, ax = axs[1])

axs[1].set_xlabel('Monocytes Initial Proportion (%)')
axs[1].set_ylabel('Concentration (%)')
# axs[1].grid(True)
# axs[1].set_axisbelow(True)

plt.xticks(range(len([str(x) for x in df_conc_exp['mono']])), [str(x) for x in df_conc_exp['mono']])

# plt.show()
plt.savefig('heterologous_allPatients.svg')  



stop = timeit.default_timer()
print(stop - start)  
