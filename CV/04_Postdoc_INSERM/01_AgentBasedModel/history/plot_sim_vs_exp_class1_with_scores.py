#!/usr/bin/python
# -*- coding: utf-8 -*-
# jupyter notebook --browser='C:/Program Files (x86)/Google/Chrome/Application/chrome.exe %s'


import timeit
start = timeit.default_timer()

import sys
import numpy
import pandas as pd
import matplotlib.pyplot as plt

## usage : $ python scripts_for_alpha_distrib/plot_sim_vs_exp.py bestvia16_stocha_0

## plot simulations vs experimental data

exp_file_path = "A:\\Downloads\\Projects\\workFromHome\\Projects\\ABM2021\\paper\\revision\\\classes_ABM_2D_9patients_1"
# simu_file_path = "A:\\Downloads\\Projects\\workFromHome\\Projects\\ABM2021\\20212201\\figures\\plots\\BehaviorSpace\\NLC-CLL-revisions\\10patients"
simu_file_name = sys.argv[1]

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
    step = int(line[22])
    viability = float(line[24])
    remainingCellRatio = float(line[25])
    if run_number not in viability_dict :
      viability_dict[run_number] = []
    viability_dict[run_number].append(viability)
    if run_number not in remaining_dict :
      remaining_dict[run_number] = []
    remaining_dict[run_number].append(remainingCellRatio) 
  

my_list = [24*d for d in range(0,14)]
my_list = [t for t in viability_exp.index.values.tolist() ]

# print(my_list)

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



# ### get optimized parameters
# params = pd.read_csv("%s.csv" % (simu_file_name), skiprows=6, sep=",", header=0, nrows = 1, usecols=[3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21])
# print(params)

# params['gui-alpha-distrib'] = params['gui-alpha-distrib'].round(decimals = 3)


# # params = params.T
# params_firstrow = params[["gui-apo-mov", "gui-need-sig-mov", "gui-cll-sens-dist", "gui-mono-sens-dist", "gui-macro-sens-dist", "gui-nlc-sens-dist"]]
# params_secondrow = params[["gui-life-init-gamma","gui-alpha-distrib",  "gui-diff-mean", "gui-diff-std", "gui-sig-init", "gui-sig-init-std", "gui-nlc-threshold"]]
# params_thirdrow = params[["gui-layers", "gui-alpha", "gui-mono-phago-eff", "gui-M-phago-eff", "gui-NLC-phago-eff", "gui-M-kill-eff"]]


# # hide axes
# fig.patch.set_visible(False)
# axes[2].axis('off')
# # axes[2].axis('tight')

# axes[3].axis('off')
# # axes[3].axis('tight')

# axes[4].axis('off')
# # axes[4].axis('tight')

# table1 = axes[2].table(cellText=params_firstrow.values, colLabels=[i[4:] for i in params_firstrow.columns], loc='center', fontsize=14).scale(xscale=1,yscale=2)

# table2 = axes[3].table(cellText=params_secondrow.values, colLabels=[i[4:] for i in params_secondrow.columns], loc='center', fontsize=14).scale(xscale=1,yscale=2)
# table3 = axes[4].table(cellText=params_thirdrow.values, colLabels=[i[4:] for i in params_thirdrow.columns], loc='center', fontsize=14).scale(xscale=1,yscale=2)


# # axes[2].table.auto_set_font_size(False)
# # table2.auto_set_font_size(False)
# # table3.auto_set_font_size(False)

# # plt.subplots_adjust(left=0.125,
# #                     bottom=0.1, 
# #                     right=0.9, 
# #                     top=0.9, 
# #                     wspace=0.1, 
# #                     hspace=0.0)
# # fig.tight_layout()

# plt.subplots_adjust(top=0.95)

# plt.subplots_adjust(hspace=0.0,wspace=0.1)



# plt.show()


# plt.show()

plt.savefig('%s_model_fit_with_scores.pdf' % (simu_file_name))
plt.savefig('%s_model_fit_with_scores.png' % (simu_file_name))


stop = timeit.default_timer()
print(stop - start)  


