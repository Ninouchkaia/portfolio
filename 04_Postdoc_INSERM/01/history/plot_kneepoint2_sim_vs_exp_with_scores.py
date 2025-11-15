#!/usr/bin/python
# -*- coding: utf-8 -*-
# jupyter notebook --browser='C:/Program Files (x86)/Google/Chrome/Application/chrome.exe %s'


import timeit
start = timeit.default_timer()

import sys
import numpy
import pandas as pd
import matplotlib.pyplot as plt


## plot simulations vs experimental data

patient = sys.argv[1]
param_set = sys.argv[2]

simu_file_path = "%s\\BehaviorSpace\\%s.csv" % (patient,param_set)
# print('simu_file_path', simu_file_path)

## get patient experimental data
patient_data = pd.read_csv('%s\\%s.csv' % (patient,patient)) #, index_col=0) 
patient_data['Day'] = patient_data['Day'].apply(lambda x: x * 24)
 
### get simulation data

viability_dict, remaining_dict= {}, {}
with open(simu_file_path, 'r') as file_read :
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
	

my_list = [24*d for d in range(0,15)]
my_list = [t for t in patient_data['Day']]


df_viability_simu = pd.DataFrame.from_dict(viability_dict)
df_remaining_simu = pd.DataFrame.from_dict(remaining_dict)

filtered_df_viability_simu = df_viability_simu[df_viability_simu.index.isin(my_list)]
filtered_df_remaining_simu = df_remaining_simu[df_remaining_simu.index.isin(my_list)]

# print(filtered_df_viability_simu)

### calculate NRMSE and squared-R for the simulations fitness to experimental values

import sklearn.metrics
nrms_via = sklearn.metrics.mean_squared_error(patient_data['%s_viability' % patient], filtered_df_viability_simu.mean(axis=1), squared=False) / (max(patient_data['%s_viability' % patient]) - min(patient_data['%s_viability' % patient]))
nrms_conc = sklearn.metrics.mean_squared_error(patient_data['%s_concentration' % patient], filtered_df_remaining_simu.mean(axis=1), squared=False) / (max(patient_data['%s_concentration' % patient]) - min(patient_data['%s_concentration' % patient]))
nrms_via = round(nrms_via, 2)
nrms_conc = round(nrms_conc, 2)

print(nrms_via, nrms_conc)


squared_r_via = sklearn.metrics.r2_score(patient_data['%s_viability' % patient], filtered_df_viability_simu.mean(axis=1))
squared_r_conc = sklearn.metrics.r2_score(patient_data['%s_concentration' % patient], filtered_df_remaining_simu.mean(axis=1))
squared_r_via = round(squared_r_via, 2)
squared_r_conc = round(squared_r_conc, 2)
print(squared_r_via, squared_r_conc)



### plot the figure
fig, axes = plt.subplots(nrows=2, ncols=1, subplot_kw={'ylim': (30,105), 'xlim' : (-0.5,14)})
fig.suptitle('Model fitting for %s' % param_set, fontsize=12)

filtered_df_viability_simu.plot(ax=axes[0], legend=False, color='red',zorder=0,ylim=(50,105), xlim =(-0.5,14))
filtered_df_remaining_simu.plot(ax=axes[1], legend=False, color='red',zorder=0, ylim= (30,105), xlim =(-0.5,14))

### add experimental
patient_data.plot(x='Day', y='%s_viability' % patient, ax=axes[0], legend=False, color='black')
patient_data.plot.scatter (x='Day', y='%s_viability' % patient, ax=axes[0], legend=False, color='black')
patient_data.plot(x='Day', y='%s_concentration' % patient, ax=axes[1], legend=False, color='black')
patient_data.plot.scatter (x='Day', y='%s_concentration' % patient, ax=axes[1], legend=False, color='black')

# axes[0].set_xlabel('Day')
x_axis = axes[0].axes.get_xaxis()
x_label = x_axis.get_label()
x_label.set_visible(False)
axes[0].set_ylabel('Viability (%)')
axes[0].set_xticks(range(0,24*15, 24))
axes[0].set_xticklabels(range(0,15))

axes[1].set_xlabel('Day')
axes[1].set_ylabel('Concentration Ratio (%)')
axes[1].set_xticks(range(0,24*15, 24))
axes[1].set_xticklabels(range(0,15))


leg = axes[1].legend(['Simulation', 'Exp'])#, loc='upper right', bbox_to_anchor=(1.3,1.3), )


LH = leg.legendHandles
LH[0].set_linewidth(8)
LH[1].set_color('black') 

# Show NRMSE and squared-R results
axes[0].text(0.1, 0.1, 'NRMSE_via = %s\n R²_via = %s' % (nrms_via,squared_r_via), horizontalalignment='left', verticalalignment='baseline', transform=axes[0].transAxes)
axes[1].text(0.1, 0.1, 'NRMSE_conc = %s\n R²_conc = %s' % (nrms_conc,squared_r_conc), horizontalalignment='left', verticalalignment='baseline', transform=axes[1].transAxes)
# plt.show()

# plt.savefig('%s\\BehaviorSpace\\%s_%s_model_fit_with_scores.pdf' % (patient, patient, param_set))
# plt.savefig('%s\\BehaviorSpace\\%s_%s_model_fit_with_scores.png' % (patient, patient, param_set))


plt.savefig('plots_with_scores_kneepoint2\\%s_%s_model_fit_with_scores.pdf' % (patient, param_set))
plt.savefig('plots_with_scores_kneepoint2\\%s_%s_model_fit_with_scores.png' % (patient, param_set))

stop = timeit.default_timer()
print(stop - start)  


# ### get simulation data
# sim_data = pd.read_csv(simu_file_path, skiprows=6, sep=",", header=0)
# print(sim_data['[run number]'])

# runs = sim_data['[run number]'].unique()
# print(sorted(runs))

# for run in runs :
# 	sim_viability = pd.DataFrame()