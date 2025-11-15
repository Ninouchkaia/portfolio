#!/usr/bin/python
# -*- coding: utf-8 -*-
# jupyter notebook --browser='C:/Program Files (x86)/Google/Chrome/Application/chrome.exe %s'


import timeit
start = timeit.default_timer()

import sys
import numpy
import pandas as pd
import matplotlib.pyplot as plt

## usage : $ python plot_sim_vs_exp_with_scores_class2_averaged.py class2_averaged_stocha

exp_file_path = "A:\\Downloads\\Projects\\workFromHome\\Projects\\ABM2021\\paper\\revision\\classes_ABM_2D_9patients_1"
# simu_file_path = "A:\\Downloads\\Projects\\workFromHome\\Projects\\ABM2021\\20212201\\figures\\plots\\BehaviorSpace\\NLC-CLL-revisions\\10patients"
simu_file_name = sys.argv[1]

## get patients experimental data
viability_exp = pd.read_csv('%s\\class2_viability.tsv' % exp_file_path, sep='\t') 
concentration_exp = pd.read_csv('%s\\class2_concentration.tsv' % exp_file_path, sep='\t') 

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

fig, axes = plt.subplots(nrows=2, ncols=1)
fig.suptitle('Model fitting for %s' % simu_file_name, fontsize=12)
 

# ### plot experimental
viability_exp.reset_index().plot(x='Day', y='mean', ax=axes[0], yerr='std',capsize=4, legend=False, color='black')
viability_exp.reset_index().plot.scatter (x='Day', y='mean', ax=axes[0], legend=False, color='black')
concentration_exp.reset_index().plot(x='Day', y='mean', ax=axes[1], yerr='std',capsize=4, legend=False, color='black')
concentration_exp.reset_index().plot.scatter (x='Day', y='mean', ax=axes[1], legend=False, color='black')

### get simulation data
viability_dict, remaining_dict= {}, {}
with open("%s.csv" % (simu_file_name), 'r') as file_read :
  data = file_read.readlines()
  for line in data[7:] :
    line = line.replace('\"', '').split(',')
    run_number = int(line[0])
    step = int(line[20])
    viability = float(line[22])
    remainingCellRatio = float(line[23])
    if run_number not in viability_dict :
      viability_dict[run_number] = []
    viability_dict[run_number].append(viability)
    if run_number not in remaining_dict :
      remaining_dict[run_number] = []
    remaining_dict[run_number].append(remainingCellRatio) 
  

# my_list = [24*d for d in range(0,14)]
my_list = [t for t in list(viability_exp.index.values)]

df_viability_simu = pd.DataFrame.from_dict(viability_dict)
df_remaining_simu = pd.DataFrame.from_dict(remaining_dict)

filtered_df_viability_simu = df_viability_simu[df_viability_simu.index.isin(my_list)]
filtered_df_remaining_simu = df_remaining_simu[df_remaining_simu.index.isin(my_list)]

# print(filtered_df_viability_simu)

### calculate NRMSE and squared-R for the simulations fitness to experimental values

import sklearn.metrics

# print(viability_exp['mean'])
# print(filtered_df_viability_simu.mean(axis=1))
nrms_via = sklearn.metrics.mean_squared_error(viability_exp['mean'], filtered_df_viability_simu.mean(axis=1), squared=False) / (max(viability_exp['mean']) - min(viability_exp['mean']))
nrms_conc = sklearn.metrics.mean_squared_error(concentration_exp['mean'], filtered_df_remaining_simu.mean(axis=1), squared=False) / (max(concentration_exp['mean']) - min(concentration_exp['mean']))
nrms_via = round(nrms_via, 2)
nrms_conc = round(nrms_conc, 2)

print(nrms_via, nrms_conc)


squared_r_via = sklearn.metrics.r2_score(viability_exp['mean'], filtered_df_viability_simu.mean(axis=1))
squared_r_conc = sklearn.metrics.r2_score(concentration_exp['mean'], filtered_df_remaining_simu.mean(axis=1))
squared_r_via = round(squared_r_via, 2)
squared_r_conc = round(squared_r_conc, 2)
print(squared_r_via, squared_r_conc)


filtered_df_viability_simu.plot(ax=axes[0],legend=False, color='r', zorder=0, ylim=(50,105), xlim =(-0.5,14))#, ylim=(00,100), color='orange')
x_axis = axes[0].axes.get_xaxis()
x_label = x_axis.get_label()
x_label.set_visible(False)

axes[0].set_ylabel('Viability (%)')
# axes[0].set_xticks(range(0,25*14, 24))
axes[0].set_xticks(range(0,24*15, 24))
axes[0].set_xticklabels(range(0,15))

filtered_df_remaining_simu.plot(figsize=(8, 8),ax=axes[1], legend=False, color='r', zorder=0, ylim= (30,105), xlim =(-0.5,14))#, ylim=(00,110), color='orange')
axes[1].set_xlabel('Day')
axes[1].set_ylabel('Concentration ratio (%)')
axes[1].set_xticks(range(0,24*15, 24))
axes[1].set_xticklabels(range(0,15))



leg = axes[1].legend(['Exp', 'Simulation'], loc='best')#, bbox_to_anchor=(1.1,1.1), )


LH = leg.legendHandles
LH[-1].set_linewidth(20)
LH[1].set_linewidth(2)
LH[1].set_color('black') 
LH[-1].set_color('red') 

# Show NRMSE and squared-R results
axes[0].text(0.1, 0.1, 'NRMSE_via = %s\n R²_via = %s' % (nrms_via,squared_r_via), horizontalalignment='left', verticalalignment='baseline', transform=axes[0].transAxes)
axes[1].text(0.1, 0.1, 'NRMSE_conc = %s\n R²_conc = %s' % (nrms_conc,squared_r_conc), horizontalalignment='left', verticalalignment='baseline', transform=axes[1].transAxes)


plt.savefig('%s_model_fit_with_scores.pdf' % (simu_file_name))
plt.savefig('%s_model_fit_with_scores.png' % (simu_file_name))

stop = timeit.default_timer()
print(stop - start)  


# ### get simulation data
# sim_data = pd.read_csv(simu_file_path, skiprows=6, sep=",", header=0)
# print(sim_data['[run number]'])

# runs = sim_data['[run number]'].unique()
# print(sorted(runs))

# for run in runs :
#   sim_viability = pd.DataFrame()