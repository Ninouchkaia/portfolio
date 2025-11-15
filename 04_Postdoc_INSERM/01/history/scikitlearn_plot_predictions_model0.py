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
import statistics
import numpy as np

df_via_exp = pd.read_csv('allPatients_prediction-via-exp.tsv', sep='\t')
df_via_exp['viability_mean'] = df_via_exp.mean(axis=1)
df_via_exp['viability_std'] = df_via_exp.std(axis=1)
df_via_exp = df_via_exp.iloc[::-1]

df_conc_exp = pd.read_csv('allPatients_prediction-conc-exp.tsv', sep='\t')
df_conc_exp['concentration_mean'] = df_conc_exp.mean(axis=1)
df_conc_exp['concentration_std'] = df_conc_exp.std(axis=1)
df_conc_exp = df_conc_exp.iloc[::-1]


sim_file = "ABM_2D_9patients_0 prediction-day9-varying-mono-init-table"

with open("%s.csv" % sim_file, 'r') as file_read :
	data = file_read.readlines()
	myViabilityDict = {}
	myremainingCellRatioDict = {}
	for line in data[7:] :
		line = line.replace('\"', '').split(',')
		monoInit = float(line[1])
		viability = float(line[23])
		remainingCellRatio = float(line[24])
		if monoInit not in myViabilityDict :
			myViabilityDict[monoInit] = [viability]
		else :
			myViabilityDict[monoInit].append(viability)
		if monoInit not in myremainingCellRatioDict :
			myremainingCellRatioDict[monoInit] = [remainingCellRatio]
		else :
			myremainingCellRatioDict[monoInit].append(remainingCellRatio)

myViabilityDict = dict((sorted(myViabilityDict.items())))
myremainingCellRatioDict = dict((sorted(myremainingCellRatioDict.items())))
# print(myremainingCellRatioDict)

### calculate the means
myViabilityMeans = {}
for monoInit in myViabilityDict :
	meanViability = statistics.mean(myViabilityDict[monoInit])
	myViabilityMeans[monoInit] = meanViability
print(myViabilityMeans)

myConcentrationMeans = {}
for monoInit in myremainingCellRatioDict :
	meanConcentration = statistics.mean(myremainingCellRatioDict[monoInit])
	myConcentrationMeans[monoInit] = meanConcentration



### calculate the stderrors
myViabilitystd = {}
for monoInit in myViabilityDict :
	stdViability = statistics.stdev(myViabilityDict[monoInit])
	myViabilitystd[monoInit] = stdViability

myConcentrationstd = {}
for monoInit in myremainingCellRatioDict :
	stdConcentration = statistics.stdev(myremainingCellRatioDict[monoInit])
	myConcentrationstd[monoInit] = stdConcentration

## create lists for the plot
monoInitList = []
viabilityMeanList, viabilityStdList = [], []
for monoInit in myViabilityDict :
	# print('monoInit',monoInit)
	monoInitList.append(monoInit)
	viabilityMeanList.append(myViabilityMeans[monoInit])
	viabilityStdList.append(myViabilitystd[monoInit])

print(viabilityMeanList)

concentrationMeanList, concentrationStdList = [], []
for monoInit in myremainingCellRatioDict :
  # print(monoInit)
  concentrationMeanList.append(myConcentrationMeans[monoInit])
  concentrationStdList.append(myConcentrationstd[monoInit])

Monocyte_initial_proportion = np.array(monoInitList)
x_pos = np.arange(len(Monocyte_initial_proportion))

viability_mean = np.array(viabilityMeanList)
viability_std = np.array(viabilityStdList)

concentration_mean = np.array(concentrationMeanList)
concentration_std = np.array(concentrationStdList)

viability_mean = list(viability_mean)
viability_std = list(viability_std)
concentration_mean = list(concentration_mean)
concentration_std = list(concentration_std)

exp_via_means_list = df_via_exp["viability_mean"]
exp_conc_means_list = df_conc_exp["concentration_mean"]

exp_via_std_list = df_via_exp["viability_std"]
exp_conc_std_list = df_conc_exp["concentration_std"]



fig, axs = plt.subplots(2,1)
# fig.suptitle('B-CLL viability and concentration\n with varying monocytes initial proportions (3 patients)\n Experimental vs. Predictions', fontsize=10, y=1.02)
# fig.suptitle('Model Predictions \n %s' % file, fontsize=10, y=1.02)
fig.suptitle('Model Predictions in heterologous co-cultures at Day 9', fontsize=14, y=0.95)
x_labels = [str(x) for x in df_via_exp['mono']]
plt.xticks(range(len(x_labels)), x_labels, fontsize=8)

x = np.arange(len(x_labels))  # the label locations
width = 0.4  # the width of the bars

axs[0].bar (x - width/2, exp_via_means_list, width, color=['black']*len(x_labels), label='Exp')
axs[0].errorbar(x - width/2, exp_via_means_list, exp_via_std_list, ecolor = 'grey', fmt='none', capsize = 2)
axs[0].bar (x + width/2, viability_mean, width, color=['red']*len(x_labels), label='Prediction')
axs[0].errorbar(x + width/2, viability_mean, viability_std, ecolor = 'grey', fmt='none', capsize = 2)
# axs[0].set_xlabel('Monocytes Initial Proportion')
axs[0].set_ylabel('Viability (%)')


plt.sca(axs[0])
plt.xticks(range(len(x_labels)), x_labels, fontsize=8)
axs[1].bar (x - width/2, exp_conc_means_list, width, color=['black']*len(x_labels), label='Exp')
axs[1].errorbar(x - width/2, exp_conc_means_list, exp_conc_std_list, ecolor = 'grey', fmt='none', capsize = 2)
axs[1].bar (x + width/2, concentration_mean, width, color=['red']*len(x_labels), label='Prediction')
axs[1].errorbar(x + width/2, concentration_mean, concentration_std, ecolor = 'grey', fmt='none', capsize = 2)
axs[1].set_xlabel('Monocytes Initial Proportion (%)')
axs[1].set_ylabel('Concentration (%)')

axs[1].legend()#loc='upper right', bbox_to_anchor=(1.3,1.3))


import math
via_errors_pred = [via_simu_i - via_exp_i for via_simu_i, via_exp_i in zip(viability_mean, exp_via_means_list)]
print(via_errors_pred)
via_errors_pred2 = [i ** 2 for i in via_errors_pred]
print(via_errors_pred2)

conc_errors_pred = [conc_simu_i - conc_exp_i for conc_simu_i, conc_exp_i in zip(concentration_mean, exp_conc_means_list)]
print(conc_errors_pred)
conc_errors_pred2 = [i ** 2 for i in conc_errors_pred]
print(conc_errors_pred2)

## calculate the RMSE and normalized RMSE for viability
import sklearn.metrics
nrms_via = sklearn.metrics.mean_squared_error(exp_via_means_list, viability_mean, squared=False) / (max(exp_via_means_list) - min(exp_via_means_list))
nrms_conc = sklearn.metrics.mean_squared_error(exp_conc_means_list, concentration_mean, squared=False) / (max(exp_conc_means_list) - min(exp_conc_means_list))
nrms_via = round(nrms_via, 2)
nrms_conc = round(nrms_conc, 2)

print(nrms_via, nrms_conc)


squared_r_via = sklearn.metrics.r2_score(exp_via_means_list, viability_mean)
squared_r_conc = sklearn.metrics.r2_score(exp_conc_means_list, concentration_mean)
squared_r_via = round(squared_r_via, 2)
squared_r_conc = round(squared_r_conc, 2)
print(squared_r_via, squared_r_conc)

# Show NRMSE and squared-R results
# (0, 0) being the lower left of the axes and (1, 1) the upper right.
axs[0].text(0.1, 0.775, 'NRMSE_via = %s\n R²_via = %s' % (nrms_via,squared_r_via), horizontalalignment='left', verticalalignment='baseline', transform=axs[0].transAxes)
axs[1].text(0.4, 0.775, 'NRMSE_conc = %s\n R²_conc = %s' % (nrms_conc,squared_r_conc), horizontalalignment='left', verticalalignment='baseline', transform=axs[1].transAxes)

# set_matplotlib_formats('pdf')
# # # set_matplotlib_formats('svg')
# set_matplotlib_formats('png')

# plt.show()

# plt.savefig("pred_Model0.png")
plt.savefig("scikitlearn_pred_Model0.svg")
plt.savefig("scikitlearn_pred_Model0.pdf")
plt.savefig("scikitlearn_pred_Model0.png")
# files.download("pred_%s.pdf" % patient) 



stop = timeit.default_timer()
print(stop - start)  

