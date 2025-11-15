#!/usr/bin/python
# -*- coding: utf-8 -*-
# jupyter notebook --browser='C:/Program Files (x86)/Google/Chrome/Application/chrome.exe %s'


import timeit
start = timeit.default_timer()

## usage  $ python plot_sensitivity_analysis.py perturb-gui-alpha-distrib 


import sys
import pandas as pd
import matplotlib.pyplot as plt

# simu_file_name = sys.argv[1] 

## plot simulations vs experimental data
simu_file_name = "ABM_2D_9patients_1_class1_%s" % sys.argv[1]

param = ("-").join(sys.argv[1].split("-")[1:])

exp_file_path = "A:\\Downloads\\Projects\\workFromHome\\Projects\\ABM2021\\paper\\revision\\classes_ABM_2D_9patients_1"

## get patients experimental data
viability_exp = pd.read_csv('%s\\class1_viability.tsv' % exp_file_path, sep='\t') 
concentration_exp = pd.read_csv('%s\\class1_concentration.tsv' % exp_file_path, sep='\t') 

viability_exp['Day'] = viability_exp['Day'].apply(lambda x: x * 24)
concentration_exp['Day'] = concentration_exp['Day'].apply(lambda x: x * 24)

viability_exp.set_index('Day', inplace=True)
concentration_exp.set_index('Day', inplace=True)

print(viability_exp)

viability_exp['mean'] = viability_exp.mean(axis=1)
viability_exp['std'] = viability_exp.std(axis=1)
print(viability_exp)

concentration_exp['mean'] = concentration_exp.mean(axis=1)
concentration_exp['std'] = concentration_exp.std(axis=1)


### get the simulated data
df = pd.read_csv("%s.csv" % simu_file_name, skiprows = 6)

df.sort_values(by=['[step]'])
df = df[df["[step]"]%24 == 0]


### plot
fig, axes = plt.subplots(nrows=2, ncols=1)
fig.suptitle('ABM_2D_Class1_%s' % sys.argv[1], fontsize=14)

groups = df.groupby(['%s' % param, '[step]']).agg(['mean', 'std'])

viability = groups["getViability"].unstack(0)
viability.plot(None, 'mean', yerr='std', ax=axes[0], legend=False, ylim=(0,105), xlim =(-0.5,14))
axes[0].set_xticks(range(0,24*15, 24))
axes[0].set_xticklabels(range(0,15))
x_axis = axes[0].axes.get_xaxis()
x_label = x_axis.get_label()
x_label.set_visible(False)
axes[0].set_ylabel('Viability (%)')


concentration = groups["getRemainingCellRatio"].unstack(0)
concentration.plot(None, 'mean', yerr='std', ax=axes[1], legend=True, ylim= (0,105), xlim =(-0.5,14))
axes[1].set_xticks(range(0,24*15, 24))
axes[1].set_xticklabels(range(0,15))
axes[1].set_xlabel('Day')
axes[1].set_ylabel('Concentration ratio (%)')


# axes[0].legend(loc='lower right', prop={'size': 7})#, bbox_to_anchor=(1.3,1.3), )
axes[1].legend(loc='center left', bbox_to_anchor=(1, 0.5))




# ### plot experimental
viability_exp.reset_index().plot(x='Day', y='mean', ax=axes[0], yerr='std',capsize=4, legend=False, color='black')
viability_exp.reset_index().plot.scatter (x='Day', y='mean', ax=axes[0], legend=False, color='black')
concentration_exp.reset_index().plot(x='Day', y='mean', ax=axes[1], yerr='std',capsize=4, legend=False, color='black')
concentration_exp.reset_index().plot.scatter (x='Day', y='mean', ax=axes[1], legend=False, color='black')

plt.plot()

# plt.show()


plt.savefig("sensitivity_analysis_class1\\%s_class1.pdf" % sys.argv[1])
plt.savefig("sensitivity_analysis_class1\\%s_class1.png" % sys.argv[1])


stop = timeit.default_timer()
print(stop - start)