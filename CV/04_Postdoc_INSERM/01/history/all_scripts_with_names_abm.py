
RMSE.py 
#!/usr/bin/python
# -*- coding: utf-8 -*-
# jupyter notebook --browser='C:/Program Files (x86)/Google/Chrome/Application/chrome.exe %s'


import timeit
start = timeit.default_timer()

import sys
import numpy
import pandas as pd
import math
import statistics

## usage : $ python scripts/RMSE.py 

## calculate RMSE between simulations vs experimental data

# patient = sys.argv[1]


patients_list_with_mono = ['CRE1704-1.1%','LAU1405-2.5%','DES2105-1.25%','ORE1706-0.68%','CAS1802-1.04%', 'GER160522-0.45%','REI230522-0.95%','PUJ240522-0.34%','LAR300522-0.21%','CAZ310522-3.48%']
patients_list = []

for pat in patients_list_with_mono :
  patient_name = pat.split("-")[0]
  patients_list.append(patient_name)
# print(patients_list)

param_sets = ['stocha_best_via', 'stocha_best_conc', 'stocha_knee_point']

df_NRMSE_via_max_min = pd.DataFrame(columns = param_sets, index=patients_list)
df_NRMSE_via_mean = pd.DataFrame(columns = param_sets, index=patients_list)
df_NRMSE_via_stdev = pd.DataFrame(columns = param_sets, index=patients_list)

df_NRMSE_conc_max_min = pd.DataFrame(columns = param_sets, index=patients_list)
df_NRMSE_conc_mean = pd.DataFrame(columns = param_sets, index=patients_list)
df_NRMSE_conc_stdev = pd.DataFrame(columns = param_sets, index=patients_list)

df_NRMSE_sum_max_min = pd.DataFrame(columns = param_sets, index=patients_list)
df_NRMSE_sum_mean = pd.DataFrame(columns = param_sets, index=patients_list)
df_NRMSE_sum_stdev = pd.DataFrame(columns = param_sets, index=patients_list)

for patient in patients_list :
  print("patient: %s
" % patient)
  for param_set in param_sets :
    print("param_set: %s
" % param_set)
    simu_file_path = "%s\BehaviorSpace\%s.csv" % (patient, param_set) #(patient,param_set)
    # print('simu_file_path', simu_file_path)

    ## get patient experimental data
    patient_data = pd.read_csv('%s\%s.csv' % (patient,patient), index_col=0) #, index_col=0) 
    #patient_data['Day'] = patient_data['Day'].apply(lambda x: x * 24)
    # print(patient_data)
     
    ### get simulation data
    viability_dict, remaining_dict= {}, {}
    with open(simu_file_path, 'r') as file_read :
      data = file_read.readlines()
      for line in data[7:] :
        line = line.replace('\"', '').split(',')
        run_number = int(line[0])
        step = int(line[21])
        viability = float(line[23])
        remainingCellRatio = float(line[24])
        if run_number not in viability_dict :
          viability_dict[run_number] = []
        viability_dict[run_number].append(viability)
        if run_number not in remaining_dict :
          remaining_dict[run_number] = []
        remaining_dict[run_number].append(remainingCellRatio)	
    	

    my_list = [24*d for d in list(patient_data.index.values)]
    # print(my_list)

    df_viability_simu = pd.DataFrame.from_dict(viability_dict)
    df_remaining_simu = pd.DataFrame.from_dict(remaining_dict)

    ### keep only the time points for which we also have experimental data
    filtered_df_viability_simu = df_viability_simu[df_viability_simu.index.isin(my_list)]
    filtered_df_remaining_simu = df_remaining_simu[df_remaining_simu.index.isin(my_list)]

    # print(filtered_df_viability_simu)

    ## calculate the RMSE and normalized RMSE for viability
    # Normalized RMSE = RMSE / (max value – min value)
    mean_via_simu = filtered_df_viability_simu.mean(axis=1).transpose()
    print("mean_via_simu
", mean_via_simu)
    mean_via_exp = patient_data['%s_viability' % patient].transpose()
    print("mean_via_exp
", mean_via_exp)

    via_errors = (mean_via_simu - mean_via_exp).transpose()
    print(via_errors)
    via_errors2 = via_errors * via_errors
    sum_via_errors = via_errors2[0].sum()
    RMSE_via = math.sqrt(sum_via_errors / len(via_errors))
    print("RMSE Via for %s = %s" % (patient,RMSE_via))

    normalized_RMSE_via_max_min = RMSE_via / (max(mean_via_exp) - min(mean_via_exp))
    print("Normalized RMSE Via by (max-min) for %s = %s" % (patient,normalized_RMSE_via_max_min))

    normalized_RMSE_via_mean = RMSE_via / statistics.mean(list(mean_via_exp))
    print("Normalized RMSE Via by mean for %s = %s" % (patient,normalized_RMSE_via_mean))

    normalized_RMSE_via_stdev = RMSE_via / statistics.stdev(list(mean_via_exp))
    print("Normalized RMSE Via by stdev for %s = %s" % (patient,normalized_RMSE_via_stdev))

    df_NRMSE_via_max_min[param_set][patient] = normalized_RMSE_via_max_min
    df_NRMSE_via_mean[param_set][patient] = normalized_RMSE_via_mean
    df_NRMSE_via_stdev[param_set][patient] = normalized_RMSE_via_stdev


    ## calculate the RMSE for concentration
    mean_conc_simu = filtered_df_remaining_simu.mean(axis=1).transpose()
    print("mean_conc_simu
", mean_conc_simu)
    mean_conc_exp = patient_data['%s_concentration' % patient].transpose()
    print("mean_conc_exp
", mean_conc_exp)
    conc_errors = (mean_conc_simu - mean_conc_exp).transpose()
    conc_errors2 = conc_errors * conc_errors
    sum_conc_errors = conc_errors2[0].sum()

    RMSE_conc = math.sqrt(sum_conc_errors / len(conc_errors))
    print("RMSE Conc for %s = %s" % (patient,RMSE_conc))

    normalized_RMSE_conc_max_min = RMSE_conc / (max(mean_conc_exp) - min(mean_conc_exp))
    print("Normalized RMSE Conc by (max-min) for %s = %s" % (patient,normalized_RMSE_conc_max_min))

    normalized_RMSE_conc_mean = RMSE_conc / statistics.mean(list(mean_conc_exp))
    print("Normalized RMSE Conc by mean for %s = %s" % (patient,normalized_RMSE_conc_mean))

    normalized_RMSE_conc_stdev = RMSE_conc / statistics.stdev(list(mean_conc_exp))
    print("Normalized RMSE Conc by stdev for %s = %s" % (patient,normalized_RMSE_conc_stdev))

    df_NRMSE_conc_max_min[param_set][patient] = normalized_RMSE_conc_max_min
    df_NRMSE_conc_mean[param_set][patient] = normalized_RMSE_conc_mean
    df_NRMSE_conc_stdev[param_set][patient] = normalized_RMSE_conc_stdev

    df_NRMSE_sum_max_min[param_set][patient] = normalized_RMSE_via_max_min +  normalized_RMSE_conc_max_min
    df_NRMSE_sum_mean[param_set][patient] = normalized_RMSE_via_mean + normalized_RMSE_conc_mean
    df_NRMSE_sum_stdev[param_set][patient] = normalized_RMSE_via_stdev + normalized_RMSE_conc_stdev

df_NRMSE_via_max_min.to_csv('NRMSE_via_max_min.tsv', sep='	', index=True)  
df_NRMSE_via_mean.to_csv('NRMSE_via_mean.tsv', sep='	', index=True)  
df_NRMSE_via_stdev.to_csv('NRMSE_via_stdev.tsv', sep='	', index=True)  

df_NRMSE_conc_max_min.to_csv('NRMSE_conc_max_min.tsv', sep='	', index=True)  
df_NRMSE_conc_mean.to_csv('NRMSE_conc_mean.tsv', sep='	', index=True)  
df_NRMSE_conc_stdev.to_csv('NRMSE_conc_stdev.tsv', sep='	', index=True)  

df_NRMSE_sum_max_min.to_csv('NRMSE_sum_max_min.tsv', sep='	', index=True)  
df_NRMSE_sum_mean.to_csv('NRMSE_sum_mean.tsv', sep='	', index=True)  
df_NRMSE_sum_stdev.to_csv('NRMSE_sum_stdev.tsv', sep='	', index=True)  



stop = timeit.default_timer()
print(stop - start)  

aggregateData.py 
#!/usr/bin/python
# -*- coding: utf-8 -*-
# jupyter notebook --browser='C:/Program Files (x86)/Google/Chrome/Application/chrome.exe %s'


import timeit
start = timeit.default_timer()



import sys
import os
import csv
# with open("outputs_ABM_9_2.txt", 'a+') as file_write:
with open(sys.argv[2], 'a+') as file_write:
	# rootdir = 'A:\Downloads\Projects\workFromHome\Projects\ABM2021\explorations\ABM_9_2\ABM_9_2'
	rootdir = sys.argv[1]
	for subdir, dirs, files in os.walk(rootdir):
		for file in files :
			print(file)
			path_to_file = rootdir + "\" + file
			if file == 'population1.csv' :
				with open(path_to_file, 'r') as file_read:
					data = file_read.readlines()
					for line in data :
						file_write.write(line)
			else :
				with open(path_to_file, 'r') as file_read:
					data = file_read.readlines()
					for line in data[1:] :
						file_write.write(line)

stop = timeit.default_timer()
print(stop - start)  
		

aggregate_patients_param.py 
#!/usr/bin/python
# -*- coding: utf-8 -*-
# jupyter notebook --browser='C:/Program Files (x86)/Google/Chrome/Application/chrome.exe %s'
import timeit
start = timeit.default_timer()


import os
import pandas as pd


patients_list_with_mono = ['CRE1704-1.1%','LAU1405-2.5%','DES2105-1.25%','ORE1706-0.68%','CAS1802-1.04%', 'GER160522-0.45%','REI230522-0.95%','PUJ240522-0.34%','LAR300522-0.21%','CAZ310522-3.48%']


patients_list = []

for patient in patients_list_with_mono :
	patient_name = patient.split("-")[0]
	patients_list.append(patient_name)

print(patients_list)



stop = timeit.default_timer()
print(stop - start) 

all_scripts_with_names_abm.py 

RMSE.py 
#!/usr/bin/python
# -*- coding: utf-8 -*-
# jupyter notebook --browser='C:/Program Files (x86)/Google/Chrome/Application/chrome.exe %s'


import timeit
start = timeit.default_timer()

import sys
import numpy
import pandas as pd
import math
import statistics

## usage : $ python scripts/RMSE.py 

## calculate RMSE between simulations vs experimental data

# patient = sys.argv[1]


patients_list_with_mono = ['CRE1704-1.1%','LAU1405-2.5%','DES2105-1.25%','ORE1706-0.68%','CAS1802-1.04%', 'GER160522-0.45%','REI230522-0.95%','PUJ240522-0.34%','LAR300522-0.21%','CAZ310522-3.48%']
patients_list = []

for pat in patients_list_with_mono :
  patient_name = pat.split("-")[0]
  patients_list.append(patient_name)
# print(patients_list)

param_sets = ['stocha_best_via', 'stocha_best_conc', 'stocha_knee_point']

df_NRMSE_via_max_min = pd.DataFrame(columns = param_sets, index=patients_list)
df_NRMSE_via_mean = pd.DataFrame(columns = param_sets, index=patients_list)
df_NRMSE_via_stdev = pd.DataFrame(columns = param_sets, index=patients_list)

df_NRMSE_conc_max_min = pd.DataFrame(columns = param_sets, index=patients_list)
df_NRMSE_conc_mean = pd.DataFrame(columns = param_sets, index=patients_list)
df_NRMSE_conc_stdev = pd.DataFrame(columns = param_sets, index=patients_list)

df_NRMSE_sum_max_min = pd.DataFrame(columns = param_sets, index=patients_list)
df_NRMSE_sum_mean = pd.DataFrame(columns = param_sets, index=patients_list)
df_NRMSE_sum_stdev = pd.DataFrame(columns = param_sets, index=patients_list)

for patient in patients_list :
  print("patient: %s
" % patient)
  for param_set in param_sets :
    print("param_set: %s
" % param_set)
    simu_file_path = "%s\BehaviorSpace\%s.csv" % (patient, param_set) #(patient,param_set)
    # print('simu_file_path', simu_file_path)

    ## get patient experimental data
    patient_data = pd.read_csv('%s\%s.csv' % (patient,patient), index_col=0) #, index_col=0) 
    #patient_data['Day'] = patient_data['Day'].apply(lambda x: x * 24)
    # print(patient_data)
     
    ### get simulation data
    viability_dict, remaining_dict= {}, {}
    with open(simu_file_path, 'r') as file_read :
      data = file_read.readlines()
      for line in data[7:] :
        line = line.replace('\"', '').split(',')
        run_number = int(line[0])
        step = int(line[21])
        viability = float(line[23])
        remainingCellRatio = float(line[24])
        if run_number not in viability_dict :
          viability_dict[run_number] = []
        viability_dict[run_number].append(viability)
        if run_number not in remaining_dict :
          remaining_dict[run_number] = []
        remaining_dict[run_number].append(remainingCellRatio)	
    	

    my_list = [24*d for d in list(patient_data.index.values)]
    # print(my_list)

    df_viability_simu = pd.DataFrame.from_dict(viability_dict)
    df_remaining_simu = pd.DataFrame.from_dict(remaining_dict)

    ### keep only the time points for which we also have experimental data
    filtered_df_viability_simu = df_viability_simu[df_viability_simu.index.isin(my_list)]
    filtered_df_remaining_simu = df_remaining_simu[df_remaining_simu.index.isin(my_list)]

    # print(filtered_df_viability_simu)

    ## calculate the RMSE and normalized RMSE for viability
    # Normalized RMSE = RMSE / (max value – min value)
    mean_via_simu = filtered_df_viability_simu.mean(axis=1).transpose()
    print("mean_via_simu
", mean_via_simu)
    mean_via_exp = patient_data['%s_viability' % patient].transpose()
    print("mean_via_exp
", mean_via_exp)

    via_errors = (mean_via_simu - mean_via_exp).transpose()
    print(via_errors)
    via_errors2 = via_errors * via_errors
    sum_via_errors = via_errors2[0].sum()
    RMSE_via = math.sqrt(sum_via_errors / len(via_errors))
    print("RMSE Via for %s = %s" % (patient,RMSE_via))

    normalized_RMSE_via_max_min = RMSE_via / (max(mean_via_exp) - min(mean_via_exp))
    print("Normalized RMSE Via by (max-min) for %s = %s" % (patient,normalized_RMSE_via_max_min))

    normalized_RMSE_via_mean = RMSE_via / statistics.mean(list(mean_via_exp))
    print("Normalized RMSE Via by mean for %s = %s" % (patient,normalized_RMSE_via_mean))

    normalized_RMSE_via_stdev = RMSE_via / statistics.stdev(list(mean_via_exp))
    print("Normalized RMSE Via by stdev for %s = %s" % (patient,normalized_RMSE_via_stdev))

    df_NRMSE_via_max_min[param_set][patient] = normalized_RMSE_via_max_min
    df_NRMSE_via_mean[param_set][patient] = normalized_RMSE_via_mean
    df_NRMSE_via_stdev[param_set][patient] = normalized_RMSE_via_stdev


    ## calculate the RMSE for concentration
    mean_conc_simu = filtered_df_remaining_simu.mean(axis=1).transpose()
    print("mean_conc_simu
", mean_conc_simu)
    mean_conc_exp = patient_data['%s_concentration' % patient].transpose()
    print("mean_conc_exp
", mean_conc_exp)
    conc_errors = (mean_conc_simu - mean_conc_exp).transpose()
    conc_errors2 = conc_errors * conc_errors
    sum_conc_errors = conc_errors2[0].sum()

    RMSE_conc = math.sqrt(sum_conc_errors / len(conc_errors))
    print("RMSE Conc for %s = %s" % (patient,RMSE_conc))

    normalized_RMSE_conc_max_min = RMSE_conc / (max(mean_conc_exp) - min(mean_conc_exp))
    print("Normalized RMSE Conc by (max-min) for %s = %s" % (patient,normalized_RMSE_conc_max_min))

    normalized_RMSE_conc_mean = RMSE_conc / statistics.mean(list(mean_conc_exp))
    print("Normalized RMSE Conc by mean for %s = %s" % (patient,normalized_RMSE_conc_mean))

    normalized_RMSE_conc_stdev = RMSE_conc / statistics.stdev(list(mean_conc_exp))
    print("Normalized RMSE Conc by stdev for %s = %s" % (patient,normalized_RMSE_conc_stdev))

    df_NRMSE_conc_max_min[param_set][patient] = normalized_RMSE_conc_max_min
    df_NRMSE_conc_mean[param_set][patient] = normalized_RMSE_conc_mean
    df_NRMSE_conc_stdev[param_set][patient] = normalized_RMSE_conc_stdev

    df_NRMSE_sum_max_min[param_set][patient] = normalized_RMSE_via_max_min +  normalized_RMSE_conc_max_min
    df_NRMSE_sum_mean[param_set][patient] = normalized_RMSE_via_mean + normalized_RMSE_conc_mean
    df_NRMSE_sum_stdev[param_set][patient] = normalized_RMSE_via_stdev + normalized_RMSE_conc_stdev

df_NRMSE_via_max_min.to_csv('NRMSE_via_max_min.tsv', sep='	', index=True)  
df_NRMSE_via_mean.to_csv('NRMSE_via_mean.tsv', sep='	', index=True)  
df_NRMSE_via_stdev.to_csv('NRMSE_via_stdev.tsv', sep='	', index=True)  

df_NRMSE_conc_max_min.to_csv('NRMSE_conc_max_min.tsv', sep='	', index=True)  
df_NRMSE_conc_mean.to_csv('NRMSE_conc_mean.tsv', sep='	', index=True)  
df_NRMSE_conc_stdev.to_csv('NRMSE_conc_stdev.tsv', sep='	', index=True)  

df_NRMSE_sum_max_min.to_csv('NRMSE_sum_max_min.tsv', sep='	', index=True)  
df_NRMSE_sum_mean.to_csv('NRMSE_sum_mean.tsv', sep='	', index=True)  
df_NRMSE_sum_stdev.to_csv('NRMSE_sum_stdev.tsv', sep='	', index=True)  



stop = timeit.default_timer()
print(stop - start)  

aggregateData.py 
#!/usr/bin/python
# -*- coding: utf-8 -*-
# jupyter notebook --browser='C:/Program Files (x86)/Google/Chrome/Application/chrome.exe %s'


import timeit
start = timeit.default_timer()



import sys
import os
import csv
# with open("outputs_ABM_9_2.txt", 'a+') as file_write:
with open(sys.argv[2], 'a+') as file_write:
	# rootdir = 'A:\Downloads\Projects\workFromHome\Projects\ABM2021xplorations\ABM_9_2\ABM_9_2'
	rootdir = sys.argv[1]
	for subdir, dirs, files in os.walk(rootdir):
		for file in files :
			print(file)
			path_to_file = rootdir + "\" + file
			if file == 'population1.csv' :
				with open(path_to_file, 'r') as file_read:
					data = file_read.readlines()
					for line in data :
						file_write.write(line)
			else :
				with open(path_to_file, 'r') as file_read:
					data = file_read.readlines()
					for line in data[1:] :
						file_write.write(line)

stop = timeit.default_timer()
print(stop - start)  
		

aggregate_patients_param.py 
#!/usr/bin/python
# -*- coding: utf-8 -*-
# jupyter notebook --browser='C:/Program Files (x86)/Google/Chrome/Application/chrome.exe %s'
import timeit
start = timeit.default_timer()


import os
import pandas as pd


patients_list_with_mono = ['CRE1704-1.1%','LAU1405-2.5%','DES2105-1.25%','ORE1706-0.68%','CAS1802-1.04%', 'GER160522-0.45%','REI230522-0.95%','PUJ240522-0.34%','LAR300522-0.21%','CAZ310522-3.48%']


patients_list = []

for patient in patients_list_with_mono :
	patient_name = patient.split("-")[0]
	patients_list.append(patient_name)

print(patients_list)



stop = timeit.default_timer()
print(stop - start) 

averaged_simu_plot_shell.py 
#!/usr/bin/python
# -*- coding: utf-8 -*-
# jupyter notebook --browser='C:/Program Files (x86)/Google/Chrome/Application/chrome.exe %s'


import timeit
start = timeit.default_timer()

patients_list_with_mono = ['CRE1704-1.1%','LAU1405-2.5%','DES2105-1.25%','ORE1706-0.68%','CAS1802-1.04%', 'GER160522-0.45%','REI230522-0.95%','PUJ240522-0.34%','LAR300522-0.21%','CAZ310522-3.48%']
patients_list = []

for patient in patients_list_with_mono :
	patient_name = patient.split("-")[0]
	patients_list.append(patient_name)
print(patients_list)

with open("patient_command_averaged_class2_plot_with_scores.sh", 'w') as file_write :
	file_write.write("#!/bin/bash
")
	for patient in patients_list : 
		file_write.write("python scripts/plot_sim_vs_exp_with_scores.py %s averaged_class2_simu_%s
" % (patient,patient))



stop = timeit.default_timer()
print(stop - start)  

averaged_simu_shell.py 
#!/usr/bin/python
# -*- coding: utf-8 -*-
# jupyter notebook --browser='C:/Program Files (x86)/Google/Chrome/Application/chrome.exe %s'


import timeit
start = timeit.default_timer()

patients_list_with_mono = ['CRE1704-1.1%','LAU1405-2.5%','DES2105-1.25%','ORE1706-0.68%','CAS1802-1.04%', 'GER160522-0.45%','REI230522-0.95%','PUJ240522-0.34%','LAR300522-0.21%','CAZ310522-3.48%']
patients_list = []

for patient in patients_list_with_mono :
	patient_name = patient.split("-")[0]
	patients_list.append(patient_name)
print(patients_list)

with open("patient_command_averaged_class1_simu.sh", 'w') as file_write :
	file_write.write("#!/bin/bash
")
	for patient in patients_list : 
		file_write.write("/I/Program\ Files/NetLogo\ 6.1.0/netlogo-headless.bat --model %s/ABM_2D_%s.nlogo --setup-file class1_averaged.xml --experiment averaged_simu_%s --table %s/BehaviorSpace/averaged_class1_simu_%s.csv --threads 4

" % (patient,patient,patient,patient,patient))



stop = timeit.default_timer()
print(stop - start)  

commands_plot_sim_vs_exp_with_scores.py 
#!/usr/bin/python
# -*- coding: utf-8 -*-
# jupyter notebook --browser='C:/Program Files (x86)/Google/Chrome/Application/chrome.exe %s'


import timeit
start = timeit.default_timer()

patients_list_with_mono = ['CRE1704-1.1%','LAU1405-2.5%','DES2105-1.25%','ORE1706-0.68%','CAS1802-1.04%', 'GER160522-0.45%','REI230522-0.95%','PUJ240522-0.34%','LAR300522-0.21%','CAZ310522-3.48%']
patients_list = []

for patient in patients_list_with_mono :
	patient_name = patient.split("-")[0]
	patients_list.append(patient_name)
print(patients_list)

with open("commands_plots_with_scores.sh", 'w') as file_write :
	file_write.write("#!/bin/bash
")
	for patient in patients_list : 
		file_write.write("echo \"This is a shell script for patient %s\"
" % patient)
		file_write.write("python scripts/plot_sim_vs_exp_with_scores.py %s stocha_best_via
" % patient)
		file_write.write("python scripts/plot_sim_vs_exp_with_scores.py %s stocha_best_conc
" % patient)
		file_write.write("python scripts/plot_sim_vs_exp_with_scores.py %s stocha_knee_point
" % patient)



stop = timeit.default_timer()
print(stop - start)  



## usage : $ python scripts/plot_sim_vs_exp_with_scores.py CAS1802 best_via

copy_for_git_paretoFrontGenericStochastic.py 

import math
import sys
import pylab as plt
import numpy as np

inputFile = sys.argv[1]
outputFile = sys.argv[2]

with open(inputFile, 'r') as file_read :
# with open("A:\Downloads\Projects\workFromHome\Projects\ABM2021\explorations\outputs_ABM_9_2.txt", 'r') as file_read :
# with open("A:\Downloads\Projects\workFromHome\Projects\ABM2021\explorations\ABM_8_4\ABM_8_4\population20000.csv", 'r') as file_read :

	my_points = []
	data = file_read.readlines()
	for line in data[1:] :
		line = line.replace("
","").split(",")
		apo = int(line[0])
		needSig = int(line[1])
		layers = int(line[2])
		alpha = int(line[3])
		monoPhago = int(line[4])
		NLCPhago = int(line[5])
		M2Phago = int(line[6])
		M2Kill = int(line[7])
		cllDist = int(line[8])
		MonoDist = int(line[9])
		nlcDist = int(line[10])
		macroDist = int(line[11])
		nlcThreshold = int(line[12])
		signalInitMean = int(line[13])
		signalInitStd = int(line[14])
		diffTime = int(line[15])
		diffInitStd = int(line[16])
		gammaLifeInit = int(line[17])
		alphaDistrib = float(line[18])
		delta_fitness_via = float(line[19])
		delta_fitness_conc = float(line[20])
		euclid = math.sqrt(delta_fitness_via*delta_fitness_via + delta_fitness_conc * delta_fitness_conc)
		my_points.append([delta_fitness_via,delta_fitness_conc, euclid, 1, 
			apo, needSig, layers, alpha, 
			monoPhago, NLCPhago, M2Phago, M2Kill, 
			cllDist, MonoDist, nlcDist, macroDist, nlcThreshold, signalInitMean, signalInitStd,
			diffTime, diffInitStd, gammaLifeInit, alphaDistrib])
# print(len(my_points))

b_set = set(tuple(x) for x in my_points)
my_points_singles = [ list(x) for x in b_set ]
# print(len(my_points_singles))


my_points_sorted = sorted(my_points_singles, key=lambda x: x[2])

for i in range(0, len(my_points_sorted)) :
	x1 = my_points_sorted[i][0]
	y1 = my_points_sorted[i][1]

	for j in range(0, len(my_points_sorted)) :
		x2 = my_points_sorted[j][0]
		y2 = my_points_sorted[j][1]
		if not ((x1 == x2) and (y1 == y2)) :
			if (x2 <= x1) and (y2 <= y1) :
				my_points_sorted[i][3] = 0
				break

pareto_front = []
for point in my_points_sorted :
	if point[3] == 1 :
		pareto_front.append((point[0],point[1], point[4:]))


fig, ax = plt.subplots()

x_val = [x[0] for x in my_points_sorted]
y_val = [x[1] for x in my_points_sorted]
x_val_pareto = [x[0] for x in pareto_front]
y_val_pareto = [x[1] for x in pareto_front]

# plt.scatter(x_val, y_val,color='black',s=0.75)
# plt.scatter(x_val_pareto, y_val_pareto,color='red',s=0.75)

plt.scatter(x_val, y_val,color='black',s=3)
plt.scatter(x_val_pareto, y_val_pareto,color='red',s=3)
ax.set_xlabel(r'$\Delta Viability Fitness$', fontsize=15)
ax.set_ylabel(r'$\Delta ConcentrationFitness$', fontsize=15)
ax.set_title('Pareto front')
ax.grid(True)
plt.tight_layout()
plt.savefig('%s.png' % outputFile, bbox_inches='tight')

# plt.show()

print("len(pareto_front)",len(pareto_front))

pareto_front_sorted = sorted(pareto_front, key=lambda x: x[0])
with open("%s.txt" % outputFile, 'w') as file_write :
	file_write.write("delta_fitness_via,delta_fitness_conc, apo, needSig, layers, alpha, monoPhago, NLCPhago, M2Phago, M2Kill, cllDist, MonoDist, nlcDist, macroDist, nlcThreshold, signalInitMean, signalInitStd, diffTime, diffInitStd, LifeInitGamma, alphaDistrib
")
	for sets in pareto_front_sorted :
		line = (",".join(str(x) for x in sets)).replace("[","").replace("]","")
		file_write.write(line)
		file_write.write("
")
 

extract_param_sets_from_pareto.py 
#!/usr/bin/python
# -*- coding: utf-8 -*-
# jupyter notebook --browser='C:/Program Files (x86)/Google/Chrome/Application/chrome.exe %s'


import timeit
start = timeit.default_timer()

import sys
import numpy as np
import math
import pandas as pd

my_file = sys.argv[1]

pareto_front = pd.read_csv("%s" % my_file)
pareto_front['distances'] = np.sqrt(pareto_front['delta_fitness_via'] * pareto_front['delta_fitness_via'] + pareto_front['delta_fitness_conc'] * pareto_front['delta_fitness_conc'])

best_via_set = pareto_front[pareto_front.delta_fitness_via == pareto_front.delta_fitness_via.min()]
best_conc_set = pareto_front[pareto_front.delta_fitness_conc == pareto_front.delta_fitness_conc.min()]
knee_point_set = pareto_front[pareto_front.distances == pareto_front.distances.min()]

best_params = best_via_set
best_params = best_params.append(knee_point_set)
best_params = best_params.append(best_conc_set)
best_params = best_params.drop('distances', 1)

best_params.insert(loc=0, value=['best_via_set', 'knee_point_set', 'best_conc_set'], column="set")

best_params.to_csv('best_param_sets_%s.tsv' % my_file[:-4 ], sep='	', index=False)  



stop = timeit.default_timer()
print(stop - start)  

extract_param_sets_from_pareto_adapted.py 
#!/usr/bin/python
# -*- coding: utf-8 -*-
# jupyter notebook --browser='C:/Program Files (x86)/Google/Chrome/Application/chrome.exe %s'


import timeit
start = timeit.default_timer()

import sys
import numpy as np
import math
import pandas as pd

my_file = sys.argv[1]
output = sys.argv[2]

pareto_front = pd.read_csv("%s" % my_file)
pareto_front['distances'] = np.sqrt(pareto_front['delta_fitness_via'] * pareto_front['delta_fitness_via'] + pareto_front['delta_fitness_conc'] * pareto_front['delta_fitness_conc'])

best_via_set = pareto_front[pareto_front.delta_fitness_via == pareto_front.delta_fitness_via.min()]
best_conc_set = pareto_front[pareto_front.delta_fitness_conc == pareto_front.delta_fitness_conc.min()]
knee_point_set = pareto_front[pareto_front.distances == pareto_front.distances.min()]

best_params = best_via_set
best_params = best_params.append(knee_point_set)
best_params = best_params.append(best_conc_set)
best_params = best_params.drop('distances', 1)

best_params.insert(loc=0, value=['best_via_set', 'knee_point_set', 'best_conc_set'], column="set")

best_params.to_csv('%s' % output, sep='	', index=False)  



stop = timeit.default_timer()
print(stop - start)  

heterologous.py 
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


# df = pd.read_csv('allPatients-heterologous.tsv', sep='	')
# print(df)

# sns.catplot(x='mono', y='viability', hue='patient', data=df, kind='bar')
# plt.grid()  #just add this

# plt.show()


df_via_exp = pd.read_csv('allPatients_prediction-via-exp.tsv', sep='	')
df_via_exp['viability_mean'] = df_via_exp.mean(axis=1)
df_via_exp['viability_std'] = df_via_exp.std(axis=1)
df_via_exp = df_via_exp.iloc[::-1]

df_conc_exp = pd.read_csv('allPatients_prediction-conc-exp.tsv', sep='	')
df_conc_exp['concentration_mean'] = df_conc_exp.mean(axis=1)
df_conc_exp['concentration_std'] = df_conc_exp.std(axis=1)
df_conc_exp = df_conc_exp.iloc[::-1]

fig, axs = plt.subplots(2,1)
fig.suptitle('B-CLL cells survival at Day 9 in heterologous co-cultures
with varying initial monocytes proportions', fontsize=14, y=1.02)

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

kneepoint0_by_patient_simu_shell.py 
#!/usr/bin/python
# -*- coding: utf-8 -*-
# jupyter notebook --browser='C:/Program Files (x86)/Google/Chrome/Application/chrome.exe %s'


import timeit
start = timeit.default_timer()

patients_list_with_mono = ['CRE1704-1.1%','LAU1405-2.5%','DES2105-1.25%','ORE1706-0.68%','CAS1802-1.04%', 'GER160522-0.45%','REI230522-0.95%','PUJ240522-0.34%','LAR300522-0.21%','CAZ310522-3.48%']
patients_list = []

for patient in patients_list_with_mono :
	patient_name = patient.split("-")[0]
	patients_list.append(patient_name)
print(patients_list)

with open("kneepoint_by_patient_simu.sh", 'w') as file_write :
	file_write.write("#!/bin/bash
")
	for patient in patients_list : 
		file_write.write("/I/Program\ Files/NetLogo\ 6.1.0/netlogo-headless.bat --model %s/ABM_2D_%s.nlogo --setup-file kneepoint0_simu_by_patient.xml --experiment kneepoint0_simu_%s --table %s/BehaviorSpace/kneepoint0_simu_%s.csv --threads 4

" % (patient,patient,patient,patient,patient))



stop = timeit.default_timer()
print(stop - start)  

kneepoint0_simu_plot_shell.py 
#!/usr/bin/python
# -*- coding: utf-8 -*-
# jupyter notebook --browser='C:/Program Files (x86)/Google/Chrome/Application/chrome.exe %s'


import timeit
start = timeit.default_timer()

patients_list_with_mono = ['CRE1704-1.1%','LAU1405-2.5%','DES2105-1.25%','ORE1706-0.68%','CAS1802-1.04%', 'GER160522-0.45%','REI230522-0.95%','PUJ240522-0.34%','LAR300522-0.21%','CAZ310522-3.48%']
patients_list = []

for patient in patients_list_with_mono :
	patient_name = patient.split("-")[0]
	patients_list.append(patient_name)
print(patients_list)

with open("plot_kneepoint0_simu_by_patient_with_scores.sh", 'w') as file_write :
	file_write.write("#!/bin/bash
")
	for patient in patients_list : 
		file_write.write("python scripts/plot_kneepoint0_sim_vs_exp_with_scores.py %s kneepoint0_simu_%s
" % (patient,patient))



stop = timeit.default_timer()
print(stop - start)  

kneepoint1_class1_by_patient_simu_shell.py 
#!/usr/bin/python
# -*- coding: utf-8 -*-
# jupyter notebook --browser='C:/Program Files (x86)/Google/Chrome/Application/chrome.exe %s'


import timeit
start = timeit.default_timer()

patients_list_with_mono = ['CRE1704-1.1%','LAU1405-2.5%','DES2105-1.25%','ORE1706-0.68%','CAS1802-1.04%', 'GER160522-0.45%','REI230522-0.95%','PUJ240522-0.34%','LAR300522-0.21%','CAZ310522-3.48%']
patients_list = []

for patient in patients_list_with_mono :
	patient_name = patient.split("-")[0]
	patients_list.append(patient_name)
print(patients_list)

with open("kneepoint1_class1_by_patient_simu.sh", 'w') as file_write :
	file_write.write("#!/bin/bash
")
	for patient in patients_list : 
		file_write.write("/I/Program\ Files/NetLogo\ 6.1.0/netlogo-headless.bat --model %s/ABM_2D_%s.nlogo --setup-file kneepoint1_class1_simu_by_patient.xml --experiment kneepoint1_class1_simu_%s --table %s/BehaviorSpace/kneepoint1_class1_simu_%s.csv --threads 4

" % (patient,patient,patient,patient,patient))



stop = timeit.default_timer()
print(stop - start)  

kneepoint1_class2_by_patient_simu_shell.py 
#!/usr/bin/python
# -*- coding: utf-8 -*-
# jupyter notebook --browser='C:/Program Files (x86)/Google/Chrome/Application/chrome.exe %s'


import timeit
start = timeit.default_timer()

patients_list_with_mono = ['CRE1704-1.1%','LAU1405-2.5%','DES2105-1.25%','ORE1706-0.68%','CAS1802-1.04%', 'GER160522-0.45%','REI230522-0.95%','PUJ240522-0.34%','LAR300522-0.21%','CAZ310522-3.48%']
patients_list = []

for patient in patients_list_with_mono :
	patient_name = patient.split("-")[0]
	patients_list.append(patient_name)
print(patients_list)

with open("kneepoint1_class2_by_patient_simu.sh", 'w') as file_write :
	file_write.write("#!/bin/bash
")
	for patient in patients_list : 
		file_write.write("/I/Program\ Files/NetLogo\ 6.1.0/netlogo-headless.bat --model %s/ABM_2D_%s.nlogo --setup-file kneepoint1_class2_simu_by_patient.xml --experiment kneepoint1_class2_simu_%s --table %s/BehaviorSpace/kneepoint1_class2_simu_%s.csv --threads 4

" % (patient,patient,patient,patient,patient))



stop = timeit.default_timer()
print(stop - start)  

kneepoint1_simu_plot_shell.py 
#!/usr/bin/python
# -*- coding: utf-8 -*-
# jupyter notebook --browser='C:/Program Files (x86)/Google/Chrome/Application/chrome.exe %s'


import timeit
start = timeit.default_timer()

patients_list_with_mono = ['CRE1704-1.1%','LAU1405-2.5%','DES2105-1.25%','ORE1706-0.68%','CAS1802-1.04%', 'GER160522-0.45%','REI230522-0.95%','PUJ240522-0.34%','LAR300522-0.21%','CAZ310522-3.48%']
patients_list = []

for patient in patients_list_with_mono :
	patient_name = patient.split("-")[0]
	patients_list.append(patient_name)
print(patients_list)

with open("plot_kneepoint1_simu_by_patient_with_scores.sh", 'w') as file_write :
	file_write.write("#!/bin/bash
")
	for patient in patients_list : 
		file_write.write("python scripts/plot_kneepoint1_sim_vs_exp_with_scores.py %s kneepoint1_class1_simu_%s
" % (patient,patient))



stop = timeit.default_timer()
print(stop - start)  

kneepoint2_simu_plot_shell.py 
#!/usr/bin/python
# -*- coding: utf-8 -*-
# jupyter notebook --browser='C:/Program Files (x86)/Google/Chrome/Application/chrome.exe %s'


import timeit
start = timeit.default_timer()

patients_list_with_mono = ['CRE1704-1.1%','LAU1405-2.5%','DES2105-1.25%','ORE1706-0.68%','CAS1802-1.04%', 'GER160522-0.45%','REI230522-0.95%','PUJ240522-0.34%','LAR300522-0.21%','CAZ310522-3.48%']
patients_list = []

for patient in patients_list_with_mono :
	patient_name = patient.split("-")[0]
	patients_list.append(patient_name)
print(patients_list)

with open("plot_kneepoint2_simu_by_patient_with_scores.sh", 'w') as file_write :
	file_write.write("#!/bin/bash
")
	for patient in patients_list : 
		file_write.write("python scripts/plot_kneepoint2_sim_vs_exp_with_scores.py %s kneepoint1_class2_simu_%s
" % (patient,patient))



stop = timeit.default_timer()
print(stop - start)  

kneepoint_patient_specific_simu_plot_shell.py 
#!/usr/bin/python
# -*- coding: utf-8 -*-
# jupyter notebook --browser='C:/Program Files (x86)/Google/Chrome/Application/chrome.exe %s'


import timeit
start = timeit.default_timer()

patients_list_with_mono = ['CRE1704-1.1%','LAU1405-2.5%','DES2105-1.25%','ORE1706-0.68%','CAS1802-1.04%', 'GER160522-0.45%','REI230522-0.95%','PUJ240522-0.34%','LAR300522-0.21%','CAZ310522-3.48%']
patients_list = []

for patient in patients_list_with_mono :
	patient_name = patient.split("-")[0]
	patients_list.append(patient_name)
print(patients_list)

with open("kneepoint_patient_specific_plot_with_scores.sh", 'w') as file_write :
	file_write.write("#!/bin/bash
")
	for patient in patients_list : 
		file_write.write("python scripts/plot_sim_vs_exp_with_scores_only_kneepoint_patient_specific.py %s stocha_knee_point
" % patient)



stop = timeit.default_timer()
print(stop - start)  

list_of_commands.py 
#!/usr/bin/python
# -*- coding: utf-8 -*-
# jupyter notebook --browser='C:/Program Files (x86)/Google/Chrome/Application/chrome.exe %s'


import timeit
start = timeit.default_timer()

patients_list_with_mono = ['CRE1704-1.1%','LAU1405-2.5%','DES2105-1.25%','ORE1706-0.68%','CAS1802-1.04%', 'GER160522-0.45%','REI230522-0.95%','PUJ240522-0.34%','LAR300522-0.21%','CAZ310522-3.48%']
patients_list = []

for patient in patients_list_with_mono :
	patient_name = patient.split("-")[0]
	patients_list.append(patient_name)
print(patients_list)

for patient in patients_list : 
	with open("patient_command_%s.sh" % patient, 'w') as file_write :
		file_write.write("#!/bin/bash
")
		file_write.write("echo \"This is a shell script for patient %s\"
" % patient)
		file_write.write("python scripts/aggregateData.py %s/ABM_2D_%s %s/outputs_ABM_2D_%s.txt

" % (patient,patient,patient,patient))
		file_write.write("python scripts/remove_duplicates_generic_with_filtering_keeping_only_samples.py %s/outputs_ABM_2D_%s.txt fitnessVia 50

" % (patient, patient))
		file_write.write("python scripts/copy_for_git_paretoFrontGenericStochastic.py %s/outputs_ABM_2D_%s_duplicates_removed_filtered_only_samples_kept_50.0.txt %s/pareto_ABM_2D_%s

" % (patient,patient,patient,patient))
		file_write.write("python scripts/extract_param_sets_from_pareto_adapted.py %s/pareto_ABM_2D_%s.txt %s/best_param_sets_ABM_2D_%s.tsv

" % (patient,patient,patient,patient))
		file_write.write("python scripts/parse_best_param_for_behavior_space_adapted.py %s/best_param_sets_ABM_2D_%s.tsv %s/netlogo_best_param_sets_ABM_2D_%s.txt

" % (patient,patient,patient,patient))
		file_write.write("python scripts/make_behavior_space_experiment_file.py %s/best_param_sets_ABM_2D_%s.tsv %s/experiment_file.xml

" % (patient,patient,patient))
		file_write.write("/I/Program\ Files/NetLogo\ 6.1.0/netlogo-headless.bat --model %s/ABM_2D_%s.nlogo --setup-file %s/experiment_file.xml --experiment stocha_best_via --table %s/BehaviorSpace/stocha_best_via.csv --threads 4

" % (patient,patient,patient,patient))
		file_write.write("/I/Program\ Files/NetLogo\ 6.1.0/netlogo-headless.bat --model %s/ABM_2D_%s.nlogo --setup-file %s/experiment_file.xml --experiment stocha_knee_point --table %s/BehaviorSpace/stocha_knee_point.csv --threads 4

" % (patient,patient,patient,patient))
		file_write.write("/I/Program\ Files/NetLogo\ 6.1.0/netlogo-headless.bat --model %s/ABM_2D_%s.nlogo --setup-file %s/experiment_file.xml --experiment stocha_best_conc --table %s/BehaviorSpace/stocha_best_conc.csv --threads 4

" % (patient,patient,patient,patient))
		file_write.write("python scripts/plot_sim_vs_exp.py %s stocha_best_via

" % patient)
		file_write.write("python scripts/plot_sim_vs_exp.py %s stocha_knee_point

" % patient)
		file_write.write("python scripts/plot_sim_vs_exp.py %s stocha_best_conc

" % patient)



stop = timeit.default_timer()
print(stop - start)  



# os.system('ls -l')

# import subprocess


# os.system("I:\Program Files\NetLogo\ 6.1.0\netlogo-headless.bat --model CAS1802\ABM_2D_CAS1802.nlogo --setup-file CAS1802\experiment_file.xml --experiment stocha_best_conc --table CAS1802\BehaviorSpace\stocha_best_conc.csv --threads 4")
# subprocess.call("I:\Program Files\NetLogo\ 6.1.0\netlogo-headless.bat --model CAS1802\ABM_2D_CAS1802.nlogo --setup-file CAS1802\experiment_file.xml --experiment stocha_best_conc --table CAS1802\BehaviorSpace\stocha_best_conc.csv --threads 4")
# subprocess.run(["I:/Program Files/NetLogo 6.1.0/netlogo-headless.bat" , "--model CAS1802/ABM_2D_CAS1802.nlogo",  "--setup-file CAS1802/experiment_file.xml", "--experiment stocha_best_conc" , "--table CAS1802/BehaviorSpace/stocha_best_conc.csv",  "--threads 4"])

# subprocess.run(["ls", "-lha"],shell=True)

make_averaged_simu_by_patients.py 
#!/usr/bin/python
# -*- coding: utf-8 -*-
# jupyter notebook --browser='C:/Program Files (x86)/Google/Chrome/Application/chrome.exe %s'


import timeit
start = timeit.default_timer()


import sys

my_file = sys.argv[1]
output = sys.argv[2]

patient_dict = {}
with open("patient_dict.txt", 'r') as file_read :
    data = file_read.readlines()
    for line in data[1:] :
      line = line.replace("
", "").split(" ")
      patient_name = line[0]
      gui_prop_mono_init = line[1]
      gui_prop_apo_init = line[2]
      patient_dict[patient_name] = [gui_prop_mono_init,gui_prop_apo_init]

with open(my_file, 'r') as file_read :
  with open("%s" % output, 'w') as file_write :
    file_write.write("<experiments>
")
    data = file_read.readlines()
    for line in data[1:] :
      line = line.replace("
", "").split("	")
      gui_apo_mov = line[3]
      gui_need_sig_mov = line[4]
      gui_layers = line[5]
      gui_alpha = line[6]
      gui_mono_phago_eff = line[7]
      gui_NLC_phago_eff = line[8]
      gui_M_phago_eff = line[9]
      gui_M_kill_eff = line[10]
      gui_cll_sens_dist = line[11]
      gui_mono_sens_dist = line[12]
      gui_nlc_sens_dist = line[13]
      gui_macro_sens_dist = line[14]
      gui_nlc_threshold = line[15]
      gui_sig_init = line[16]
      gui_sig_init_std = line[17]
      gui_diff_mean = line[18]
      gui_diff_std = line[19]
      gui_life_init_gamma = line[20]
      gui_alpha_distrib = line[21]
      for patient_name in patient_dict :
        file_write.write("  <experiment name=\"averaged_simu_%s\" repetitions=\"12\" runMetricsEveryStep=\"true\">
" % patient_name)
        file_write.write("    <setup>setup</setup>
")
        file_write.write("    <go>go</go>
")
        file_write.write("    <timeLimit steps=\"312\"/>
")
        file_write.write("    <metric>getSeed</metric>
")
        file_write.write("    <metric>getViability</metric>
")
        file_write.write("    <metric>getRemainingCellRatio</metric>
")

        file_write.write("    <enumeratedValueSet variable=\"gui-prop-mono-init\"><value value=\"%s\"/></enumeratedValueSet>
" % patient_dict[patient_name][0])
        file_write.write("    <enumeratedValueSet variable=\"gui-prop-apo-init\"><value value=\"%s\"/></enumeratedValueSet>
" % patient_dict[patient_name][1])       
        file_write.write("    <enumeratedValueSet variable=\"gui-apo-mov\"><value value=\"%s\"/></enumeratedValueSet>
" % gui_apo_mov)
        file_write.write("    <enumeratedValueSet variable=\"gui-need-sig-mov\"><value value=\"%s\"/></enumeratedValueSet>
" % gui_need_sig_mov)
        file_write.write("    <enumeratedValueSet variable=\"gui-layers\"><value value=\"%s\"/></enumeratedValueSet>
" % gui_layers)
        file_write.write("    <enumeratedValueSet variable=\"gui-alpha\"><value value=\"%s\"/></enumeratedValueSet>
" % gui_alpha)
        file_write.write("    <enumeratedValueSet variable=\"gui-mono-phago-eff\"><value value=\"%s\"/></enumeratedValueSet>
" % gui_mono_phago_eff)
        file_write.write("    <enumeratedValueSet variable=\"gui-NLC-phago-eff\"><value value=\"%s\"/></enumeratedValueSet>
" % gui_NLC_phago_eff)
        file_write.write("    <enumeratedValueSet variable=\"gui-M-phago-eff\"><value value=\"%s\"/></enumeratedValueSet>
" % gui_M_phago_eff)
        file_write.write("    <enumeratedValueSet variable=\"gui-M-kill-eff\"><value value=\"%s\"/></enumeratedValueSet>
" % gui_M_kill_eff)
        file_write.write("    <enumeratedValueSet variable=\"gui-cll-sens-dist\"><value value=\"%s\"/></enumeratedValueSet>
" % gui_cll_sens_dist)
        file_write.write("    <enumeratedValueSet variable=\"gui-mono-sens-dist\"><value value=\"%s\"/></enumeratedValueSet>
" % gui_mono_sens_dist)
        file_write.write("    <enumeratedValueSet variable=\"gui-nlc-sens-dist\"><value value=\"%s\"/></enumeratedValueSet>
" % gui_nlc_sens_dist)
        file_write.write("    <enumeratedValueSet variable=\"gui-macro-sens-dist\"><value value=\"%s\"/></enumeratedValueSet>
" % gui_macro_sens_dist)
        file_write.write("    <enumeratedValueSet variable=\"gui-nlc-threshold\"><value value=\"%s\"/></enumeratedValueSet>
" % gui_nlc_threshold)
        file_write.write("    <enumeratedValueSet variable=\"gui-sig-init\"><value value=\"%s\"/></enumeratedValueSet>
" % gui_sig_init)
        file_write.write("    <enumeratedValueSet variable=\"gui-sig-init-std\"><value value=\"%s\"/></enumeratedValueSet>
" % gui_sig_init_std)
        file_write.write("    <enumeratedValueSet variable=\"gui-diff-mean\"><value value=\"%s\"/></enumeratedValueSet>
" % gui_diff_mean)
        file_write.write("    <enumeratedValueSet variable=\"gui-diff-std\"><value value=\"%s\"/></enumeratedValueSet>
" % gui_diff_std)
        file_write.write("    <enumeratedValueSet variable=\"gui-life-init-gamma\"><value value=\"%s\"/></enumeratedValueSet>
" % gui_life_init_gamma)
        file_write.write("    <enumeratedValueSet variable=\"gui-alpha-distrib\"><value value=\"%s\"/></enumeratedValueSet>
" % gui_alpha_distrib)
        file_write.write("  </experiment>
")
    
    file_write.write("  </experiments>
")



stop = timeit.default_timer()
print(stop - start)  

make_behavior_space_experiment_file.py 
#!/usr/bin/python
# -*- coding: utf-8 -*-
# jupyter notebook --browser='C:/Program Files (x86)/Google/Chrome/Application/chrome.exe %s'


import timeit
start = timeit.default_timer()


### usage $ python sensitivity_analysis/make_behavior_space_experiment_file.py best_param_sets_pareto_ABM_2D_9patients_0_50.tsv sensitivity_analysis_experiment_file.xml


import sys

my_file = sys.argv[1]
output = sys.argv[2]

gui_prop_mono_init = 1.28
gui_prop_apo_init = 4.55

ranges_dict = {}
### build the dict with param names and the value ranges to perturb
ranges_dict["gui-apo-mov"] = [0,2,10]
ranges_dict["gui-need-sig-mov"] = [0,2,10]
ranges_dict["gui-layers"] = [1,1,3]
ranges_dict["gui-alpha"] = [0,50,300]
ranges_dict["gui-mono-phago-eff"] = [0,10,100]
ranges_dict["gui-NLC-phago-eff"] = [0,10,100]
ranges_dict["gui-M-phago-eff"] = [0,10,100]
ranges_dict["gui-M-kill-eff"] = [0,1,5]
ranges_dict["gui-cll-sens-dist"] = [1,1,3]
ranges_dict["gui-mono-sens-dist"] = [1,1,3]
ranges_dict["gui-nlc-sens-dist"] = [1,1,3]
ranges_dict["gui-macro-sens-dist"] = [1,1,3]
ranges_dict["gui-nlc-threshold"] = [90,20,210]
ranges_dict["gui-sig-init"] = [0,12,72]
ranges_dict["gui-sig-init-std"] = [0,8,48]
ranges_dict["gui-diff-mean"] = [48,4,72]
ranges_dict["gui-diff-std"] = [0,8,48]
ranges_dict["gui-life-init-gamma"] = [50,125,2500]
ranges_dict["gui-alpha-distrib"] = [0.1,0.05,1.0]

knee_point_dict = {}
with open(my_file, 'r') as file_read :
  data = file_read.readlines()
  for line in data[1:] :
    # print(line)
    line = line.replace("
", "").split("	")
    # print(line[0])
    if line[0] == "knee_point_set" :
      knee_point_dict["gui-apo-mov"] = line[3]
      knee_point_dict["gui-need-sig-mov"] = line[4]
      knee_point_dict["gui-layers"] = line[5]
      knee_point_dict["gui-alpha"] = line[6]
      knee_point_dict["gui-mono-phago-eff"] = line[7]
      knee_point_dict["gui-NLC-phago-eff"] = line[8]
      knee_point_dict["gui-M-phago-eff"] = line[9]
      knee_point_dict["gui-M-kill-eff"] = line[10]
      knee_point_dict["gui-cll-sens-dist"] = line[11]
      knee_point_dict["gui-mono-sens-dist"] = line[12]
      knee_point_dict["gui-nlc-sens-dist"] = line[13]
      knee_point_dict["gui-macro-sens-dist"] = line[14]
      knee_point_dict["gui-nlc-threshold"] = line[15]
      knee_point_dict["gui-sig-init"] = line[16]
      knee_point_dict["gui-sig-init-std"] = line[17]
      knee_point_dict["gui-diff-mean"] = line[18]
      knee_point_dict["gui-diff-std"] = line[19]
      knee_point_dict["gui-life-init-gamma"] = line[20]
      knee_point_dict["gui-alpha-distrib"] = line[21]

# print(knee_point_dict)

with open("%s" % output, 'w') as file_write :
    file_write.write("<experiments>
")      
    for param_name in ranges_dict :
      # print (param_name)
      ## put all the other param in a list (remaining_params) and iterate over it
      remaining_params = list(ranges_dict.keys())
      remaining_params.remove(param_name)
      # print(remaining_params)
      print("perturb-%s" % param_name)
      file_write.write("  <experiment name=\"perturb-%s\" repetitions=\"3\" runMetricsEveryStep=\"true\">
" % param_name)
      file_write.write("    <setup>setup</setup>
")
      file_write.write("    <go>go</go>
")
      file_write.write("    <timeLimit steps=\"312\"/>
")
      file_write.write("    <metric>getSeed</metric>
")
      file_write.write("    <metric>getViability</metric>
")
      file_write.write("    <metric>getRemainingCellRatio</metric>
")

      file_write.write("    <enumeratedValueSet variable=\"gui-prop-mono-init\"><value value=\"%s\"/></enumeratedValueSet>
" % gui_prop_mono_init)
      file_write.write("    <enumeratedValueSet variable=\"gui-prop-apo-init\"><value value=\"%s\"/></enumeratedValueSet>
" % gui_prop_apo_init)  

      for other_param_name in remaining_params :
        file_write.write("    <enumeratedValueSet variable=\"%s\"><value value=\"%s\"/></enumeratedValueSet>
" % (other_param_name, knee_point_dict[other_param_name]))

      file_write.write("    <steppedValueSet variable=\"%s\" first=\"%s\" step=\"%s\" last=\"%s\"/>
" % (param_name, ranges_dict[param_name][0], ranges_dict[param_name][1], ranges_dict[param_name][2]))
      file_write.write("  </experiment>
")

    file_write.write("  </experiments>
")



stop = timeit.default_timer()
print(stop - start)

make_behavior_space_experiment_file_class1.py 
#!/usr/bin/python
# -*- coding: utf-8 -*-
# jupyter notebook --browser='C:/Program Files (x86)/Google/Chrome/Application/chrome.exe %s'


import timeit
start = timeit.default_timer()


### usage $ python sensitivity_analysis/make_behavior_space_experiment_file_class1.py class1_processing/best_param_sets_pareto_ABM_2D_9patients_1_class1_50.tsv sensitivity_analysis_experiment_file_class1.xml


import sys

my_file = sys.argv[1]
output = sys.argv[2]

gui_prop_mono_init = 1.038
gui_prop_apo_init = 5.95

ranges_dict = {}
### build the dict with param names and the value ranges to perturb
ranges_dict["gui-apo-mov"] = [0,2,10]
ranges_dict["gui-need-sig-mov"] = [0,2,10]
ranges_dict["gui-layers"] = [1,1,3]
ranges_dict["gui-alpha"] = [0,50,300]
ranges_dict["gui-mono-phago-eff"] = [0,10,100]
ranges_dict["gui-NLC-phago-eff"] = [0,10,100]
ranges_dict["gui-M-phago-eff"] = [0,10,100]
ranges_dict["gui-M-kill-eff"] = [0,1,5]
ranges_dict["gui-cll-sens-dist"] = [1,1,3]
ranges_dict["gui-mono-sens-dist"] = [1,1,3]
ranges_dict["gui-nlc-sens-dist"] = [1,1,3]
ranges_dict["gui-macro-sens-dist"] = [1,1,3]
ranges_dict["gui-nlc-threshold"] = [90,20,210]
ranges_dict["gui-sig-init"] = [0,12,72]
ranges_dict["gui-sig-init-std"] = [0,8,48]
ranges_dict["gui-diff-mean"] = [48,4,72]
ranges_dict["gui-diff-std"] = [0,8,48]
ranges_dict["gui-life-init-gamma"] = [50,125,2500]
ranges_dict["gui-alpha-distrib"] = [0.1,0.05,1.0]

knee_point_dict = {}
with open(my_file, 'r') as file_read :
  data = file_read.readlines()
  for line in data[1:] :
    # print(line)
    line = line.replace("
", "").split("	")
    # print(line[0])
    if line[0] == "knee_point_set" :
      knee_point_dict["gui-apo-mov"] = line[3]
      knee_point_dict["gui-need-sig-mov"] = line[4]
      knee_point_dict["gui-layers"] = line[5]
      knee_point_dict["gui-alpha"] = line[6]
      knee_point_dict["gui-mono-phago-eff"] = line[7]
      knee_point_dict["gui-NLC-phago-eff"] = line[8]
      knee_point_dict["gui-M-phago-eff"] = line[9]
      knee_point_dict["gui-M-kill-eff"] = line[10]
      knee_point_dict["gui-cll-sens-dist"] = line[11]
      knee_point_dict["gui-mono-sens-dist"] = line[12]
      knee_point_dict["gui-nlc-sens-dist"] = line[13]
      knee_point_dict["gui-macro-sens-dist"] = line[14]
      knee_point_dict["gui-nlc-threshold"] = line[15]
      knee_point_dict["gui-sig-init"] = line[16]
      knee_point_dict["gui-sig-init-std"] = line[17]
      knee_point_dict["gui-diff-mean"] = line[18]
      knee_point_dict["gui-diff-std"] = line[19]
      knee_point_dict["gui-life-init-gamma"] = line[20]
      knee_point_dict["gui-alpha-distrib"] = line[21]

# print(knee_point_dict)

with open("%s" % output, 'w') as file_write :
    file_write.write("<experiments>
")      
    for param_name in ranges_dict :
      # print (param_name)
      ## put all the other param in a list (remaining_params) and iterate over it
      remaining_params = list(ranges_dict.keys())
      remaining_params.remove(param_name)
      # print(remaining_params)
      print("perturb-%s" % param_name)
      file_write.write("  <experiment name=\"perturb-%s\" repetitions=\"3\" runMetricsEveryStep=\"true\">
" % param_name)
      file_write.write("    <setup>setup</setup>
")
      file_write.write("    <go>go</go>
")
      file_write.write("    <timeLimit steps=\"312\"/>
")
      file_write.write("    <metric>getSeed</metric>
")
      file_write.write("    <metric>getViability</metric>
")
      file_write.write("    <metric>getRemainingCellRatio</metric>
")

      file_write.write("    <enumeratedValueSet variable=\"gui-prop-mono-init\"><value value=\"%s\"/></enumeratedValueSet>
" % gui_prop_mono_init)
      file_write.write("    <enumeratedValueSet variable=\"gui-prop-apo-init\"><value value=\"%s\"/></enumeratedValueSet>
" % gui_prop_apo_init)  

      for other_param_name in remaining_params :
        file_write.write("    <enumeratedValueSet variable=\"%s\"><value value=\"%s\"/></enumeratedValueSet>
" % (other_param_name, knee_point_dict[other_param_name]))

      file_write.write("    <steppedValueSet variable=\"%s\" first=\"%s\" step=\"%s\" last=\"%s\"/>
" % (param_name, ranges_dict[param_name][0], ranges_dict[param_name][1], ranges_dict[param_name][2]))
      file_write.write("  </experiment>
")

    file_write.write("  </experiments>
")



stop = timeit.default_timer()
print(stop - start)

make_behavior_space_experiment_file_class2.py 
#!/usr/bin/python
# -*- coding: utf-8 -*-
# jupyter notebook --browser='C:/Program Files (x86)/Google/Chrome/Application/chrome.exe %s'


import timeit
start = timeit.default_timer()


### usage $ python sensitivity_analysis/make_behavior_space_experiment_file_class2.py class2_processing/best_param_sets_pareto_ABM_2D_9patients_1_class2_50.tsv sensitivity_analysis_experiment_file_class2.xml


import sys

my_file = sys.argv[1]
output = sys.argv[2]

gui_prop_mono_init = 1.59
gui_prop_apo_init = 2.80

ranges_dict = {}
### build the dict with param names and the value ranges to perturb
ranges_dict["gui-apo-mov"] = [0,2,10]
ranges_dict["gui-need-sig-mov"] = [0,2,10]
ranges_dict["gui-layers"] = [1,1,3]
ranges_dict["gui-alpha"] = [0,50,300]
ranges_dict["gui-mono-phago-eff"] = [0,10,100]
ranges_dict["gui-NLC-phago-eff"] = [0,10,100]
ranges_dict["gui-M-phago-eff"] = [0,10,100]
ranges_dict["gui-M-kill-eff"] = [0,1,5]
ranges_dict["gui-cll-sens-dist"] = [1,1,3]
ranges_dict["gui-mono-sens-dist"] = [1,1,3]
ranges_dict["gui-nlc-sens-dist"] = [1,1,3]
ranges_dict["gui-macro-sens-dist"] = [1,1,3]
ranges_dict["gui-nlc-threshold"] = [90,20,210]
ranges_dict["gui-sig-init"] = [0,12,72]
ranges_dict["gui-sig-init-std"] = [0,8,48]
ranges_dict["gui-diff-mean"] = [48,4,72]
ranges_dict["gui-diff-std"] = [0,8,48]
ranges_dict["gui-life-init-gamma"] = [50,125,2500]
ranges_dict["gui-alpha-distrib"] = [0.1,0.05,1.0]

knee_point_dict = {}
with open(my_file, 'r') as file_read :
  data = file_read.readlines()
  for line in data[1:] :
    # print(line)
    line = line.replace("
", "").split("	")
    # print(line[0])
    if line[0] == "knee_point_set" :
      knee_point_dict["gui-apo-mov"] = line[3]
      knee_point_dict["gui-need-sig-mov"] = line[4]
      knee_point_dict["gui-layers"] = line[5]
      knee_point_dict["gui-alpha"] = line[6]
      knee_point_dict["gui-mono-phago-eff"] = line[7]
      knee_point_dict["gui-NLC-phago-eff"] = line[8]
      knee_point_dict["gui-M-phago-eff"] = line[9]
      knee_point_dict["gui-M-kill-eff"] = line[10]
      knee_point_dict["gui-cll-sens-dist"] = line[11]
      knee_point_dict["gui-mono-sens-dist"] = line[12]
      knee_point_dict["gui-nlc-sens-dist"] = line[13]
      knee_point_dict["gui-macro-sens-dist"] = line[14]
      knee_point_dict["gui-nlc-threshold"] = line[15]
      knee_point_dict["gui-sig-init"] = line[16]
      knee_point_dict["gui-sig-init-std"] = line[17]
      knee_point_dict["gui-diff-mean"] = line[18]
      knee_point_dict["gui-diff-std"] = line[19]
      knee_point_dict["gui-life-init-gamma"] = line[20]
      knee_point_dict["gui-alpha-distrib"] = line[21]

# print(knee_point_dict)

with open("%s" % output, 'w') as file_write :
    file_write.write("<experiments>
")      
    for param_name in ranges_dict :
      # print (param_name)
      ## put all the other param in a list (remaining_params) and iterate over it
      remaining_params = list(ranges_dict.keys())
      remaining_params.remove(param_name)
      # print(remaining_params)
      print("perturb-%s" % param_name)
      file_write.write("  <experiment name=\"perturb-%s\" repetitions=\"3\" runMetricsEveryStep=\"true\">
" % param_name)
      file_write.write("    <setup>setup</setup>
")
      file_write.write("    <go>go</go>
")
      file_write.write("    <timeLimit steps=\"312\"/>
")
      file_write.write("    <metric>getSeed</metric>
")
      file_write.write("    <metric>getViability</metric>
")
      file_write.write("    <metric>getRemainingCellRatio</metric>
")

      file_write.write("    <enumeratedValueSet variable=\"gui-prop-mono-init\"><value value=\"%s\"/></enumeratedValueSet>
" % gui_prop_mono_init)
      file_write.write("    <enumeratedValueSet variable=\"gui-prop-apo-init\"><value value=\"%s\"/></enumeratedValueSet>
" % gui_prop_apo_init)  

      for other_param_name in remaining_params :
        file_write.write("    <enumeratedValueSet variable=\"%s\"><value value=\"%s\"/></enumeratedValueSet>
" % (other_param_name, knee_point_dict[other_param_name]))

      file_write.write("    <steppedValueSet variable=\"%s\" first=\"%s\" step=\"%s\" last=\"%s\"/>
" % (param_name, ranges_dict[param_name][0], ranges_dict[param_name][1], ranges_dict[param_name][2]))
      file_write.write("  </experiment>
")

    file_write.write("  </experiments>
")



stop = timeit.default_timer()
print(stop - start)

make_best_sets_all_patients.py 
#!/usr/bin/python
# -*- coding: utf-8 -*-
# jupyter notebook --browser='C:/Program Files (x86)/Google/Chrome/Application/chrome.exe %s'


import timeit
start = timeit.default_timer()

import sys
import os
import pandas as pd

# my_file = sys.argv[1]
# output = sys.argv[2]

### usage : python scripts/make_best_sets_all_patients.py best_via_set 
### usage : python scripts/make_best_sets_all_patients.py best_conc_set 
### usage : python scripts/make_best_sets_all_patients.py knee_point_set 

type_of_set = sys.argv[1] ## 'best_via_set', 'best_conc_set', 'knee_point_set'

patients_list_with_mono = ['CRE1704-1.1%','LAU1405-2.5%','DES2105-1.25%','ORE1706-0.68%','CAS1802-1.04%', 'GER160522-0.45%','REI230522-0.95%','PUJ240522-0.34%','LAR300522-0.21%','CAZ310522-3.48%']
patients_list = []

for patient in patients_list_with_mono :
	patient_name = patient.split("-")[0]
	patients_list.append(patient_name)
print(patients_list)

best_set_all_patients = pd.DataFrame(columns = patients_list)
for patient in patients_list :
	df = pd.read_csv('%s\best_param_sets_ABM_2D_%s.tsv' % (patient,patient), sep='	', header=0, index_col=0)
	# print(df)
	best_set = df.loc[type_of_set]
	# print(best_via_set)
	best_set_all_patients[patient] = best_set.T

# print(best_via_set_all_patients)

best_set_all_patients.to_csv('%ss_all_patients.tsv' % type_of_set, index=True, sep='	')


stop = timeit.default_timer()
print(stop - start)  

make_files_for_git.py 
#!/usr/bin/python
# -*- coding: utf-8 -*-
# jupyter notebook --browser='C:/Program Files (x86)/Google/Chrome/Application/chrome.exe %s'
import timeit
start = timeit.default_timer()


import os
import pandas as pd

patient_dict = {}

patients_list_with_mono = ['CRE1704-1.1%','LAU1405-2.5%','DES2105-1.25%','ORE1706-0.68%','CAS1802-1.04%','REI230522-0.95%','PUJ240522-0.34%','LAR300522-0.21%','CAZ310522-3.48%']


patients_list = []

for patient in patients_list_with_mono :
	patient_name = patient.split("-")[0]
	patient_mono = float(patient.split("-")[1][:-1])
	patients_list.append(patient_name)
	patient_dict[patient_name] = [patient_mono]

print(patients_list)

os.mkdir("patient_data_for_git")
for patient in patients_list :
	patient_alias = patients_list.index(patient) + 1
	os.mkdir('patient_data_for_git\patient_%s' % patient_alias)
	print(patient)
	print(patient_alias)
	
	with open('patient_data_for_git\patient_%s\pareto_front_patient_%s.txt' % (patient_alias,patient_alias), 'w') as file_write :
		with open('%s\pareto_ABM_2D_%s.txt' % (patient,patient), 'r') as file_read :
			data = file_read.readlines()
			data[0] = "delta_fitness_via,delta_fitness_conc, apoCellsMovementProba, needSigCellsMvtProba, layersAroundNLC, antiApoBoost, monoPhagoEff, nlcPhagoEff, macroPhagoEff, macroKillEff, cllSensingDistance, monocyteSensingDistance, nlcSensingDistance, macrophageSensingDistance, nlcThreshold, signalInitMean, signalInitStd, monoDiffThreshold , monoDiffTimeStd, gammaLifeInitRate, gammaLifeInitShape
"
		file_write.write(''.join(data))

	with open('patient_data_for_git\patient_%s\NSGAII_exploration_output_patient_%s.txt' % (patient_alias,patient_alias), 'w') as file_write :
		with open('%s\outputs_ABM_2D_%s_duplicates_removed_filtered_only_samples_kept_50.0.txt' % (patient,patient), 'r') as file_read :
			data = file_read.readlines()
			data[0] = "apoCellsMovementProba,needSigCellsMvtProba,layersAroundNLC,antiApoBoost,monoPhagoEff,nlcPhagoEff,macroPhagoEff,macroKillEff,cllSensingDistance,monocyteSensingDistance,nlcSensingDistance,macrophageSensingDistance,nlcThreshold,signalInitMean,signalInitStd,monoDiffThreshold,monoDiffTimeStd,gammaLifeInitRate,gammaLifeInitShape,fitnessVia,fitnessConc,evolution$samples
"
		file_write.write(''.join(data))
	

stop = timeit.default_timer()
print(stop - start) 

make_folders.py 
#!/usr/bin/python
# -*- coding: utf-8 -*-
# jupyter notebook --browser='C:/Program Files (x86)/Google/Chrome/Application/chrome.exe %s'
import timeit
start = timeit.default_timer()


import os
import pandas as pd

patient_dict = {}

patients_list_with_mono = ['CRE1704-1.1%','LAU1405-2.5%','DES2105-1.25%','ORE1706-0.68%','CAS1802-1.04%', 'GER160522-0.45%','REI230522-0.95%','PUJ240522-0.34%','LAR300522-0.21%','CAZ310522-3.48%']


patients_list = []

for patient in patients_list_with_mono :
	patient_name = patient.split("-")[0]
	patient_mono = float(patient.split("-")[1][:-1])
	patients_list.append(patient_name)
	patient_dict[patient_name] = [patient_mono]

print(patients_list)


### fetch viability and concentration data for each patient from the fused dataframes

## viability
df_viability = pd.read_table('fused_10patients_viability_renamed.tsv', index_col='Day')
print(df_viability)
print(df_viability.columns)

### return %apo-init cells from viability at Day 0
for patient in patients_list :
	apo_init = 100 - df_viability.at[0,patient]
	print(df_viability.at[0,patient], "{:.2f}".format(apo_init))
	patient_dict[patient].append(float("{:.2f}".format(apo_init)))

### and store it into the patient_dict	
for patient in patient_dict :
	print(patient, patient_dict[patient][0], patient_dict[patient][1])



## concentration
df_concentration = pd.read_table('fused_10patients_concentration_renamed.tsv', index_col='Day')
print(df_concentration)
print(df_concentration.columns)

# patients_list = ['CRE1704']

for patient in patients_list :
	os.mkdir('%s' % patient)
	print(patient)
	viability = df_viability.loc[:,patient]
	concentration = df_concentration.loc[:,patient]
	df_patient = pd.merge(viability, concentration, on='Day', suffixes=('_viability','_concentration'))
	df_patient=df_patient.dropna(axis=0)
	print(df_patient)
	df_patient.to_csv('%s\%s.csv' % (patient,patient), index=True)
	days_measured = [str(x * 24) for x in list(df_patient.index.values)]
	time_series = "%s" % ' '.join(days_measured)
	print(time_series)

	with open('%s\ABM_2D_%s.nlogo' % (patient,patient), 'w') as file_write :
		with open('ABM_2D_10patients_16.nlogo', 'r') as file_read :
			data = file_read.readlines()
			
			for line in data :
				if '  if (member? ticks [0 24 48 72 144 168 192 216 240 311])' in line :
					print(line)
					adapted_line = '  if (member? ticks [%s])' % time_series
					file_write.write('%s
' % adapted_line)
				else :
					file_write.write(line)

	with open('%s\ABM_2D_%s.oms' % (patient,patient), 'w') as file_write :
		with open("ABM_2D_10patients_18.oms", 'r') as file_read :
			data = file_read.readlines()
			for line in data :
				if '    NetLogo6Task(workDirectory / "ABM_2D_10patients_16.nlogo", launch, embedWorkspace = false, switch3d = false, seed = mySeed) set( ' in line :
					adapted_line = '    NetLogo6Task(workDirectory / \"ABM_2D_%s.nlogo\", launch, embedWorkspace = false, switch3d = false, seed = mySeed) set( ' % patient
					file_write.write('%s
' % adapted_line)
				
				elif '    dataFile1 := (workDirectory / "filtered_fused_10patients.csv"),' in line :
					adapted_line = '    dataFile1 := (workDirectory / \"%s.csv\"),' % patient
					file_write.write('%s
' % adapted_line)
				
				elif '        evaluation = modelTask(1.2, 4.32, simuDur) -- fitnessTask,' in line :
					adapted_line = '        evaluation = modelTask(%s, %s, simuDur) -- fitnessTask,' % (patient_dict[patient][0],patient_dict[patient][1])
					file_write.write('%s
' % adapted_line)

				elif '    ) on ifb hook(workDirectory / "ABM_2D_10_patients_18", 1) // 1 save each 1 pop' in line :
					adapted_line = '    ) on ifb hook(workDirectory / \"ABM_2D_%s\", 1) // 1 save each 1 pop' % patient
					file_write.write('%s
' % adapted_line)
				else :
					file_write.write(line)


stop = timeit.default_timer()
print(stop - start) 

make_knee_points_summary.py 
#!/usr/bin/python
# -*- coding: utf-8 -*-
# jupyter notebook --browser='C:/Program Files (x86)/Google/Chrome/Application/chrome.exe %s'


import timeit
start = timeit.default_timer()

import pandas as pd

patients_list_with_mono = ['CRE1704-1.1%','LAU1405-2.5%','DES2105-1.25%','ORE1706-0.68%','CAS1802-1.04%', 'GER160522-0.45%','REI230522-0.95%','PUJ240522-0.34%','LAR300522-0.21%','CAZ310522-3.48%']

patients_list = []
for patient in patients_list_with_mono :
	patient_name = patient.split("-")[0]
	patients_list.append(patient_name)
print(patients_list)

df_summary_knee_points = pd.DataFrame(columns = patients_list)

for patient in patients_list :
	print(patient)
	my_file = "%s\best_param_sets_ABM_2D_%s.tsv" % (patient,patient)
	df = pd.read_table(my_file, index_col=0)
	df_knee_point = df.loc['knee_point_set']
	print(df_knee_point)
	df_summary_knee_points[patient] = df_knee_point

print(df_summary_knee_points)

df_summary_knee_points.to_csv("summary_knee_points.csv", index=True, sep=',')

stop = timeit.default_timer()
print(stop - start)  

make_kneepoint0_simu_by_patients.py 
#!/usr/bin/python
# -*- coding: utf-8 -*-
# jupyter notebook --browser='C:/Program Files (x86)/Google/Chrome/Application/chrome.exe %s'


import timeit
start = timeit.default_timer()


#### $ python scripts/make_kneepoint0_simu_by_patients.py kneepointset_0.tsv kneepoint0_simu_by_patient.xml

import sys

my_file = sys.argv[1]
output = sys.argv[2]

patient_dict = {}
with open("patient_dict.txt", 'r') as file_read :
    data = file_read.readlines()
    for line in data[1:] :
      line = line.replace("
", "").split(" ")
      patient_name = line[0]
      gui_prop_mono_init = line[1]
      gui_prop_apo_init = line[2]
      patient_dict[patient_name] = [gui_prop_mono_init,gui_prop_apo_init]

with open(my_file, 'r') as file_read :
  with open("%s" % output, 'w') as file_write :
    file_write.write("<experiments>
")
    data = file_read.readlines()
    for line in data[1:] :
      line = line.replace("
", "").split("	")
      gui_apo_mov = line[3]
      gui_need_sig_mov = line[4]
      gui_layers = line[5]
      gui_alpha = line[6]
      gui_mono_phago_eff = line[7]
      gui_NLC_phago_eff = line[8]
      gui_M_phago_eff = line[9]
      gui_M_kill_eff = line[10]
      gui_cll_sens_dist = line[11]
      gui_mono_sens_dist = line[12]
      gui_nlc_sens_dist = line[13]
      gui_macro_sens_dist = line[14]
      gui_nlc_threshold = line[15]
      gui_sig_init = line[16]
      gui_sig_init_std = line[17]
      gui_diff_mean = line[18]
      gui_diff_std = line[19]
      gui_life_init_gamma = line[20]
      gui_alpha_distrib = line[21]
      for patient_name in patient_dict :
        file_write.write("  <experiment name=\"kneepoint0_simu_%s\" repetitions=\"12\" runMetricsEveryStep=\"true\">
" % patient_name)
        file_write.write("    <setup>setup</setup>
")
        file_write.write("    <go>go</go>
")
        file_write.write("    <timeLimit steps=\"312\"/>
")
        file_write.write("    <metric>getSeed</metric>
")
        file_write.write("    <metric>getViability</metric>
")
        file_write.write("    <metric>getRemainingCellRatio</metric>
")

        file_write.write("    <enumeratedValueSet variable=\"gui-prop-mono-init\"><value value=\"%s\"/></enumeratedValueSet>
" % patient_dict[patient_name][0])
        file_write.write("    <enumeratedValueSet variable=\"gui-prop-apo-init\"><value value=\"%s\"/></enumeratedValueSet>
" % patient_dict[patient_name][1])       
        file_write.write("    <enumeratedValueSet variable=\"gui-apo-mov\"><value value=\"%s\"/></enumeratedValueSet>
" % gui_apo_mov)
        file_write.write("    <enumeratedValueSet variable=\"gui-need-sig-mov\"><value value=\"%s\"/></enumeratedValueSet>
" % gui_need_sig_mov)
        file_write.write("    <enumeratedValueSet variable=\"gui-layers\"><value value=\"%s\"/></enumeratedValueSet>
" % gui_layers)
        file_write.write("    <enumeratedValueSet variable=\"gui-alpha\"><value value=\"%s\"/></enumeratedValueSet>
" % gui_alpha)
        file_write.write("    <enumeratedValueSet variable=\"gui-mono-phago-eff\"><value value=\"%s\"/></enumeratedValueSet>
" % gui_mono_phago_eff)
        file_write.write("    <enumeratedValueSet variable=\"gui-NLC-phago-eff\"><value value=\"%s\"/></enumeratedValueSet>
" % gui_NLC_phago_eff)
        file_write.write("    <enumeratedValueSet variable=\"gui-M-phago-eff\"><value value=\"%s\"/></enumeratedValueSet>
" % gui_M_phago_eff)
        file_write.write("    <enumeratedValueSet variable=\"gui-M-kill-eff\"><value value=\"%s\"/></enumeratedValueSet>
" % gui_M_kill_eff)
        file_write.write("    <enumeratedValueSet variable=\"gui-cll-sens-dist\"><value value=\"%s\"/></enumeratedValueSet>
" % gui_cll_sens_dist)
        file_write.write("    <enumeratedValueSet variable=\"gui-mono-sens-dist\"><value value=\"%s\"/></enumeratedValueSet>
" % gui_mono_sens_dist)
        file_write.write("    <enumeratedValueSet variable=\"gui-nlc-sens-dist\"><value value=\"%s\"/></enumeratedValueSet>
" % gui_nlc_sens_dist)
        file_write.write("    <enumeratedValueSet variable=\"gui-macro-sens-dist\"><value value=\"%s\"/></enumeratedValueSet>
" % gui_macro_sens_dist)
        file_write.write("    <enumeratedValueSet variable=\"gui-nlc-threshold\"><value value=\"%s\"/></enumeratedValueSet>
" % gui_nlc_threshold)
        file_write.write("    <enumeratedValueSet variable=\"gui-sig-init\"><value value=\"%s\"/></enumeratedValueSet>
" % gui_sig_init)
        file_write.write("    <enumeratedValueSet variable=\"gui-sig-init-std\"><value value=\"%s\"/></enumeratedValueSet>
" % gui_sig_init_std)
        file_write.write("    <enumeratedValueSet variable=\"gui-diff-mean\"><value value=\"%s\"/></enumeratedValueSet>
" % gui_diff_mean)
        file_write.write("    <enumeratedValueSet variable=\"gui-diff-std\"><value value=\"%s\"/></enumeratedValueSet>
" % gui_diff_std)
        file_write.write("    <enumeratedValueSet variable=\"gui-life-init-gamma\"><value value=\"%s\"/></enumeratedValueSet>
" % gui_life_init_gamma)
        file_write.write("    <enumeratedValueSet variable=\"gui-alpha-distrib\"><value value=\"%s\"/></enumeratedValueSet>
" % gui_alpha_distrib)
        file_write.write("  </experiment>
")
    
    file_write.write("  </experiments>
")



stop = timeit.default_timer()
print(stop - start)  

make_kneepoint1_class1_simu_by_patients.py 
#!/usr/bin/python
# -*- coding: utf-8 -*-
# jupyter notebook --browser='C:/Program Files (x86)/Google/Chrome/Application/chrome.exe %s'


import timeit
start = timeit.default_timer()


#### $ python scripts/make_kneepoint1_class1_simu_by_patients.py kneepointset1_class1.tsv kneepoint1_class1_simu_by_patient.xml

import sys

my_file = sys.argv[1]
output = sys.argv[2]

patient_dict = {}
with open("patient_dict.txt", 'r') as file_read :
    data = file_read.readlines()
    for line in data[1:] :
      line = line.replace("
", "").split(" ")
      patient_name = line[0]
      gui_prop_mono_init = line[1]
      gui_prop_apo_init = line[2]
      patient_dict[patient_name] = [gui_prop_mono_init,gui_prop_apo_init]

with open(my_file, 'r') as file_read :
  with open("%s" % output, 'w') as file_write :
    file_write.write("<experiments>
")
    data = file_read.readlines()
    for line in data[1:] :
      line = line.replace("
", "").split("	")
      gui_apo_mov = line[3]
      gui_need_sig_mov = line[4]
      gui_layers = line[5]
      gui_alpha = line[6]
      gui_mono_phago_eff = line[7]
      gui_NLC_phago_eff = line[8]
      gui_M_phago_eff = line[9]
      gui_M_kill_eff = line[10]
      gui_cll_sens_dist = line[11]
      gui_mono_sens_dist = line[12]
      gui_nlc_sens_dist = line[13]
      gui_macro_sens_dist = line[14]
      gui_nlc_threshold = line[15]
      gui_sig_init = line[16]
      gui_sig_init_std = line[17]
      gui_diff_mean = line[18]
      gui_diff_std = line[19]
      gui_life_init_gamma = line[20]
      gui_alpha_distrib = line[21]
      for patient_name in patient_dict :
        file_write.write("  <experiment name=\"kneepoint1_class1_simu_%s\" repetitions=\"12\" runMetricsEveryStep=\"true\">
" % patient_name)
        file_write.write("    <setup>setup</setup>
")
        file_write.write("    <go>go</go>
")
        file_write.write("    <timeLimit steps=\"312\"/>
")
        file_write.write("    <metric>getSeed</metric>
")
        file_write.write("    <metric>getViability</metric>
")
        file_write.write("    <metric>getRemainingCellRatio</metric>
")

        file_write.write("    <enumeratedValueSet variable=\"gui-prop-mono-init\"><value value=\"%s\"/></enumeratedValueSet>
" % patient_dict[patient_name][0])
        file_write.write("    <enumeratedValueSet variable=\"gui-prop-apo-init\"><value value=\"%s\"/></enumeratedValueSet>
" % patient_dict[patient_name][1])       
        file_write.write("    <enumeratedValueSet variable=\"gui-apo-mov\"><value value=\"%s\"/></enumeratedValueSet>
" % gui_apo_mov)
        file_write.write("    <enumeratedValueSet variable=\"gui-need-sig-mov\"><value value=\"%s\"/></enumeratedValueSet>
" % gui_need_sig_mov)
        file_write.write("    <enumeratedValueSet variable=\"gui-layers\"><value value=\"%s\"/></enumeratedValueSet>
" % gui_layers)
        file_write.write("    <enumeratedValueSet variable=\"gui-alpha\"><value value=\"%s\"/></enumeratedValueSet>
" % gui_alpha)
        file_write.write("    <enumeratedValueSet variable=\"gui-mono-phago-eff\"><value value=\"%s\"/></enumeratedValueSet>
" % gui_mono_phago_eff)
        file_write.write("    <enumeratedValueSet variable=\"gui-NLC-phago-eff\"><value value=\"%s\"/></enumeratedValueSet>
" % gui_NLC_phago_eff)
        file_write.write("    <enumeratedValueSet variable=\"gui-M-phago-eff\"><value value=\"%s\"/></enumeratedValueSet>
" % gui_M_phago_eff)
        file_write.write("    <enumeratedValueSet variable=\"gui-M-kill-eff\"><value value=\"%s\"/></enumeratedValueSet>
" % gui_M_kill_eff)
        file_write.write("    <enumeratedValueSet variable=\"gui-cll-sens-dist\"><value value=\"%s\"/></enumeratedValueSet>
" % gui_cll_sens_dist)
        file_write.write("    <enumeratedValueSet variable=\"gui-mono-sens-dist\"><value value=\"%s\"/></enumeratedValueSet>
" % gui_mono_sens_dist)
        file_write.write("    <enumeratedValueSet variable=\"gui-nlc-sens-dist\"><value value=\"%s\"/></enumeratedValueSet>
" % gui_nlc_sens_dist)
        file_write.write("    <enumeratedValueSet variable=\"gui-macro-sens-dist\"><value value=\"%s\"/></enumeratedValueSet>
" % gui_macro_sens_dist)
        file_write.write("    <enumeratedValueSet variable=\"gui-nlc-threshold\"><value value=\"%s\"/></enumeratedValueSet>
" % gui_nlc_threshold)
        file_write.write("    <enumeratedValueSet variable=\"gui-sig-init\"><value value=\"%s\"/></enumeratedValueSet>
" % gui_sig_init)
        file_write.write("    <enumeratedValueSet variable=\"gui-sig-init-std\"><value value=\"%s\"/></enumeratedValueSet>
" % gui_sig_init_std)
        file_write.write("    <enumeratedValueSet variable=\"gui-diff-mean\"><value value=\"%s\"/></enumeratedValueSet>
" % gui_diff_mean)
        file_write.write("    <enumeratedValueSet variable=\"gui-diff-std\"><value value=\"%s\"/></enumeratedValueSet>
" % gui_diff_std)
        file_write.write("    <enumeratedValueSet variable=\"gui-life-init-gamma\"><value value=\"%s\"/></enumeratedValueSet>
" % gui_life_init_gamma)
        file_write.write("    <enumeratedValueSet variable=\"gui-alpha-distrib\"><value value=\"%s\"/></enumeratedValueSet>
" % gui_alpha_distrib)
        file_write.write("  </experiment>
")
    
    file_write.write("  </experiments>
")



stop = timeit.default_timer()
print(stop - start)  

make_kneepoint1_class2_simu_by_patients.py 
#!/usr/bin/python
# -*- coding: utf-8 -*-
# jupyter notebook --browser='C:/Program Files (x86)/Google/Chrome/Application/chrome.exe %s'


import timeit
start = timeit.default_timer()


#### $ python scripts/make_kneepoint1_class2_simu_by_patients.py kneepointset1_class2.tsv kneepoint1_class2_simu_by_patient.xml

import sys

my_file = sys.argv[1]
output = sys.argv[2]

patient_dict = {}
with open("patient_dict.txt", 'r') as file_read :
    data = file_read.readlines()
    for line in data[1:] :
      line = line.replace("
", "").split(" ")
      patient_name = line[0]
      gui_prop_mono_init = line[1]
      gui_prop_apo_init = line[2]
      patient_dict[patient_name] = [gui_prop_mono_init,gui_prop_apo_init]

with open(my_file, 'r') as file_read :
  with open("%s" % output, 'w') as file_write :
    file_write.write("<experiments>
")
    data = file_read.readlines()
    for line in data[1:] :
      line = line.replace("
", "").split("	")
      gui_apo_mov = line[3]
      gui_need_sig_mov = line[4]
      gui_layers = line[5]
      gui_alpha = line[6]
      gui_mono_phago_eff = line[7]
      gui_NLC_phago_eff = line[8]
      gui_M_phago_eff = line[9]
      gui_M_kill_eff = line[10]
      gui_cll_sens_dist = line[11]
      gui_mono_sens_dist = line[12]
      gui_nlc_sens_dist = line[13]
      gui_macro_sens_dist = line[14]
      gui_nlc_threshold = line[15]
      gui_sig_init = line[16]
      gui_sig_init_std = line[17]
      gui_diff_mean = line[18]
      gui_diff_std = line[19]
      gui_life_init_gamma = line[20]
      gui_alpha_distrib = line[21]
      for patient_name in patient_dict :
        file_write.write("  <experiment name=\"kneepoint1_class2_simu_%s\" repetitions=\"12\" runMetricsEveryStep=\"true\">
" % patient_name)
        file_write.write("    <setup>setup</setup>
")
        file_write.write("    <go>go</go>
")
        file_write.write("    <timeLimit steps=\"312\"/>
")
        file_write.write("    <metric>getSeed</metric>
")
        file_write.write("    <metric>getViability</metric>
")
        file_write.write("    <metric>getRemainingCellRatio</metric>
")

        file_write.write("    <enumeratedValueSet variable=\"gui-prop-mono-init\"><value value=\"%s\"/></enumeratedValueSet>
" % patient_dict[patient_name][0])
        file_write.write("    <enumeratedValueSet variable=\"gui-prop-apo-init\"><value value=\"%s\"/></enumeratedValueSet>
" % patient_dict[patient_name][1])       
        file_write.write("    <enumeratedValueSet variable=\"gui-apo-mov\"><value value=\"%s\"/></enumeratedValueSet>
" % gui_apo_mov)
        file_write.write("    <enumeratedValueSet variable=\"gui-need-sig-mov\"><value value=\"%s\"/></enumeratedValueSet>
" % gui_need_sig_mov)
        file_write.write("    <enumeratedValueSet variable=\"gui-layers\"><value value=\"%s\"/></enumeratedValueSet>
" % gui_layers)
        file_write.write("    <enumeratedValueSet variable=\"gui-alpha\"><value value=\"%s\"/></enumeratedValueSet>
" % gui_alpha)
        file_write.write("    <enumeratedValueSet variable=\"gui-mono-phago-eff\"><value value=\"%s\"/></enumeratedValueSet>
" % gui_mono_phago_eff)
        file_write.write("    <enumeratedValueSet variable=\"gui-NLC-phago-eff\"><value value=\"%s\"/></enumeratedValueSet>
" % gui_NLC_phago_eff)
        file_write.write("    <enumeratedValueSet variable=\"gui-M-phago-eff\"><value value=\"%s\"/></enumeratedValueSet>
" % gui_M_phago_eff)
        file_write.write("    <enumeratedValueSet variable=\"gui-M-kill-eff\"><value value=\"%s\"/></enumeratedValueSet>
" % gui_M_kill_eff)
        file_write.write("    <enumeratedValueSet variable=\"gui-cll-sens-dist\"><value value=\"%s\"/></enumeratedValueSet>
" % gui_cll_sens_dist)
        file_write.write("    <enumeratedValueSet variable=\"gui-mono-sens-dist\"><value value=\"%s\"/></enumeratedValueSet>
" % gui_mono_sens_dist)
        file_write.write("    <enumeratedValueSet variable=\"gui-nlc-sens-dist\"><value value=\"%s\"/></enumeratedValueSet>
" % gui_nlc_sens_dist)
        file_write.write("    <enumeratedValueSet variable=\"gui-macro-sens-dist\"><value value=\"%s\"/></enumeratedValueSet>
" % gui_macro_sens_dist)
        file_write.write("    <enumeratedValueSet variable=\"gui-nlc-threshold\"><value value=\"%s\"/></enumeratedValueSet>
" % gui_nlc_threshold)
        file_write.write("    <enumeratedValueSet variable=\"gui-sig-init\"><value value=\"%s\"/></enumeratedValueSet>
" % gui_sig_init)
        file_write.write("    <enumeratedValueSet variable=\"gui-sig-init-std\"><value value=\"%s\"/></enumeratedValueSet>
" % gui_sig_init_std)
        file_write.write("    <enumeratedValueSet variable=\"gui-diff-mean\"><value value=\"%s\"/></enumeratedValueSet>
" % gui_diff_mean)
        file_write.write("    <enumeratedValueSet variable=\"gui-diff-std\"><value value=\"%s\"/></enumeratedValueSet>
" % gui_diff_std)
        file_write.write("    <enumeratedValueSet variable=\"gui-life-init-gamma\"><value value=\"%s\"/></enumeratedValueSet>
" % gui_life_init_gamma)
        file_write.write("    <enumeratedValueSet variable=\"gui-alpha-distrib\"><value value=\"%s\"/></enumeratedValueSet>
" % gui_alpha_distrib)
        file_write.write("  </experiment>
")
    
    file_write.write("  </experiments>
")



stop = timeit.default_timer()
print(stop - start)  

make_shell_commands.py 
#!/usr/bin/python
# -*- coding: utf-8 -*-
# jupyter notebook --browser='C:/Program Files (x86)/Google/Chrome/Application/chrome.exe %s'


import timeit
start = timeit.default_timer()

exp_list = ["perturb-gui-apo-mov","perturb-gui-need-sig-mov","perturb-gui-layers","perturb-gui-alpha","perturb-gui-mono-phago-eff","perturb-gui-NLC-phago-eff","perturb-gui-M-phago-eff","perturb-gui-M-kill-eff","perturb-gui-cll-sens-dist","perturb-gui-mono-sens-dist","perturb-gui-nlc-sens-dist","perturb-gui-macro-sens-dist","perturb-gui-nlc-threshold","perturb-gui-sig-init","perturb-gui-sig-init-std","perturb-gui-diff-mean","perturb-gui-diff-std","perturb-gui-life-init-gamma","perturb-gui-alpha-distrib"]


with open("sensitivity_experiments.sh" , 'a+') as file_write :
	file_write.write("#!/bin/bash
")
	for exp in exp_list : 

		file_write.write("echo \"This is a shell script for exp %s\"
" % exp)
		file_write.write("/I/Program\ Files/NetLogo\ 6.1.0/netlogo-headless.bat --model /C/Users/Nina/.openmole/nina-windows/webui/projects/ABM_2D_9patients_0.nlogo --setup-file sensitivity_analysis_experiment_file.xml --experiment %s --table ABM_2D_9patients_0_%s.csv --threads 4

" % (exp,exp))



stop = timeit.default_timer()
print(stop - start)

make_shell_commands_class1.py 
#!/usr/bin/python
# -*- coding: utf-8 -*-
# jupyter notebook --browser='C:/Program Files (x86)/Google/Chrome/Application/chrome.exe %s'


import timeit
start = timeit.default_timer()

exp_list = ["perturb-gui-apo-mov","perturb-gui-need-sig-mov","perturb-gui-layers","perturb-gui-alpha","perturb-gui-mono-phago-eff","perturb-gui-NLC-phago-eff","perturb-gui-M-phago-eff","perturb-gui-M-kill-eff","perturb-gui-cll-sens-dist","perturb-gui-mono-sens-dist","perturb-gui-nlc-sens-dist","perturb-gui-macro-sens-dist","perturb-gui-nlc-threshold","perturb-gui-sig-init","perturb-gui-sig-init-std","perturb-gui-diff-mean","perturb-gui-diff-std","perturb-gui-life-init-gamma","perturb-gui-alpha-distrib"]


with open("sensitivity_experiments_class1.sh" , 'a+') as file_write :
	file_write.write("#!/bin/bash
")
	for exp in exp_list : 

		file_write.write("echo \"This is a shell script for exp %s\"
" % exp)
		file_write.write("/I/Program\ Files/NetLogo\ 6.1.0/netlogo-headless.bat --model /C/Users/Nina/.openmole/nina-windows/webui/projects/ABM_2D_9patients_1_class1.nlogo --setup-file sensitivity_analysis_experiment_file_class1.xml --experiment %s --table ABM_2D_9patients_1_class1_%s.csv --threads 4

" % (exp,exp))



stop = timeit.default_timer()
print(stop - start)

make_shell_commands_class2.py 
#!/usr/bin/python
# -*- coding: utf-8 -*-
# jupyter notebook --browser='C:/Program Files (x86)/Google/Chrome/Application/chrome.exe %s'


import timeit
start = timeit.default_timer()

exp_list = ["perturb-gui-apo-mov","perturb-gui-need-sig-mov","perturb-gui-layers","perturb-gui-alpha","perturb-gui-mono-phago-eff","perturb-gui-NLC-phago-eff","perturb-gui-M-phago-eff","perturb-gui-M-kill-eff","perturb-gui-cll-sens-dist","perturb-gui-mono-sens-dist","perturb-gui-nlc-sens-dist","perturb-gui-macro-sens-dist","perturb-gui-nlc-threshold","perturb-gui-sig-init","perturb-gui-sig-init-std","perturb-gui-diff-mean","perturb-gui-diff-std","perturb-gui-life-init-gamma","perturb-gui-alpha-distrib"]


with open("sensitivity_experiments_class2.sh" , 'a+') as file_write :
	file_write.write("#!/bin/bash
")
	for exp in exp_list : 

		file_write.write("echo \"This is a shell script for exp %s\"
" % exp)
		file_write.write("/I/Program\ Files/NetLogo\ 6.1.0/netlogo-headless.bat --model /C/Users/Nina/.openmole/nina-windows/webui/projects/ABM_2D_9patients_1_class2.nlogo --setup-file sensitivity_analysis_experiment_file_class2.xml --experiment %s --table ABM_2D_9patients_1_class2_%s.csv --threads 4

" % (exp,exp))



stop = timeit.default_timer()
print(stop - start)

parametersets_violinplots.py 
#!/usr/bin/python
# -*- coding: utf-8 -*-
# jupyter notebook --browser='C:/Program Files (x86)/Google/Chrome/Application/chrome.exe %s'


import timeit
start = timeit.default_timer()

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

parameters_sets = pd.read_csv('patient_param_sets.csv', sep=';', index_col=0) 

# sns.violinplot(data=parameters_sets["parameter"])

# df = sns.load_dataset("titanic")
# print(df)
# sns.violinplot(data=df, x="deck", y="age", hue="alive", split=True)

parameters_sets = parameters_sets.T
parameters_sets["class"] = parameters_sets["class"].astype("int")
# parameters_sets["class"] = parameters_sets["class"].astype("category")

print(parameters_sets)



# sns.set(context="paper", palette="colorblind", style="ticks")
# g = sns.FacetGrid(parameters_sets, col=" alpha", sharey=False, size=4, aspect=.5)
# g = g.map(sns.violinplot, data=parameters_sets, y=" alpha", x="class",  bw=.15, hue="class", split=False, cut=0, scale="width")
# # g.fig.get_axes()[0].legend(False)

# sns.violinplot(data=parameters_sets, y=" alpha", x="class",  bw=.15, hue="class", split=True, cut=0, scale="width")

# fig, axes = plt.subplots(nrows=20)
# i=0
# for column in parameters_sets:
#     print(column)
#     sns.violinplot(data=parameters_sets, y=column, x="class",  bw=.15, hue="class", split=False, cut=0, scale="width",ax=axes[i])
#     i=i+1


# sns.violinplot(data=parameters_sets, y=" alpha", x="class",  bw=.15, split=False, cut=0, scale="width")


for column in parameters_sets:
	print(column)
	fig, ax = plt.subplots()
	# sns.violinplot(data=parameters_sets, y=column, x="class",  bw=.15, split=False, cut=0, scale="width") ## 1
	sns.violinplot(data=parameters_sets, y=column, x="class",  bw=.30, split=False, cut=0, scale="width") ## 2
	sns.violinplot(data=parameters_sets, y=column, x="class",  bw=.60, split=False, cut=0, scale="width") ## 3
	fig.savefig('%s_distrib.pdf' % column.replace(" ","") ) 
	fig.savefig('%s_distrib.png' % column.replace(" ","") ) 



# plt.show()

stop = timeit.default_timer()
print(stop - start)  

parametersets_violinplots_9patients.py 
#!/usr/bin/python
# -*- coding: utf-8 -*-
# jupyter notebook --browser='C:/Program Files (x86)/Google/Chrome/Application/chrome.exe %s'


import timeit
start = timeit.default_timer()

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

parameters_sets = pd.read_csv('patient_param_sets_9patients_param_renamed.csv', sep=';', index_col=0) 

# sns.violinplot(data=parameters_sets["parameter"])

# df = sns.load_dataset("titanic")
# print(df)
# sns.violinplot(data=df, x="deck", y="age", hue="alive", split=True)

parameters_sets = parameters_sets.T
parameters_sets["class"] = parameters_sets["class"].astype("int")
# parameters_sets["class"] = parameters_sets["class"].astype("category")

print(parameters_sets)



# sns.set(context="paper", palette="colorblind", style="ticks")
# g = sns.FacetGrid(parameters_sets, col=" alpha", sharey=False, size=4, aspect=.5)
# g = g.map(sns.violinplot, data=parameters_sets, y=" alpha", x="class",  bw=.15, hue="class", split=False, cut=0, scale="width")
# # g.fig.get_axes()[0].legend(False)

# sns.violinplot(data=parameters_sets, y=" alpha", x="class",  bw=.15, hue="class", split=True, cut=0, scale="width")

# fig, axes = plt.subplots(nrows=20)
# i=0
# for column in parameters_sets:
#     print(column)
#     sns.violinplot(data=parameters_sets, y=column, x="class",  bw=.15, hue="class", split=False, cut=0, scale="width",ax=axes[i])
#     i=i+1


# sns.violinplot(data=parameters_sets, y=" alpha", x="class",  bw=.15, split=False, cut=0, scale="width")


for column in parameters_sets:
	print(column)
	fig, ax = plt.subplots()
	# sns.violinplot(data=parameters_sets, y=column, x="class",  bw=.15, split=False, cut=0, scale="width") ## 1
	ax = sns.violinplot(data=parameters_sets, y=column, x="class",  bw=.30, split=False, cut=0, scale="width") ## 2
	ax.set_xticklabels(["Class A", "Class B"])
	ax.set(xlabel=None)
	sns.set(font_scale=1.5)
	# sns.violinplot(data=parameters_sets, y=column, x="class",  bw=.60, split=False, cut=0, scale="width") ## 3
	fig.savefig('%s_distrib.pdf' % column.replace(" ","") ) 
	fig.savefig('%s_distrib.png' % column.replace(" ","") ) 



# plt.show()

stop = timeit.default_timer()
print(stop - start)  

parse_best_param_for_behavior_space.py 
#!/usr/bin/python
# -*- coding: utf-8 -*-
# jupyter notebook --browser='C:/Program Files (x86)/Google/Chrome/Application/chrome.exe %s'


import timeit
start = timeit.default_timer()

import sys

my_file = sys.argv[1]

with open(my_file, 'r') as file_read :
	with open("netlogo_behavior_space_%s.txt" % my_file[:-4], 'w') as file_write :
		data = file_read.readlines()
		n = 1
		for line in data[1:] :
			line = line.replace("
", "").split("	")
			gui_apo_mov = line[3]
			gui_need_sig_mov = line[4]
			gui_layers = line[5]
			gui_alpha = line[6]
			gui_mono_phago_eff = line[7]
			gui_NLC_phago_eff = line[8]
			gui_M_phago_eff = line[9]
			gui_M_kill_eff = line[10]
			gui_cll_sens_dist = line[11]
			gui_mono_sens_dist = line[12]
			gui_nlc_sens_dist = line[13]
			gui_macro_sens_dist = line[14]
			gui_nlc_threshold = line[15]
			gui_sig_init = line[16]
			gui_sig_init_std = line[17]
			gui_diff_mean = line[18]
			gui_diff_std = line[19]
			gui_life_init_gamma = line[20]
			gui_alpha_distrib = line[21]
			if (n == 1) :
				file_write.write("Best_via_set
")
			elif (n == 2) :
				file_write.write("Knee_point_set
")
			elif (n == 3) :
				file_write.write("Best_conc_set
")
			file_write.write("[\"gui-apo-mov\" %s]
" % gui_apo_mov)
			file_write.write("[\"gui-need-sig-mov\" %s]
" % gui_need_sig_mov)
			file_write.write("[\"gui-layers\" %s]
" % gui_layers)
			file_write.write("[\"gui-alpha\" %s]
" % gui_alpha)
			file_write.write("[\"gui-mono-phago-eff\" %s]
" % gui_mono_phago_eff)
			file_write.write("[\"gui-NLC-phago-eff\" %s]
" % gui_NLC_phago_eff)
			file_write.write("[\"gui-M-phago-eff\" %s]
" % gui_M_phago_eff)
			file_write.write("[\"gui-M-kill-eff\" %s]
" % gui_M_kill_eff)
			file_write.write("[\"gui-cll-sens-dist\" %s]
" % gui_cll_sens_dist)
			file_write.write("[\"gui-mono-sens-dist\" %s]
" % gui_mono_sens_dist)
			file_write.write("[\"gui-nlc-sens-dist\" %s]
" % gui_nlc_sens_dist)
			file_write.write("[\"gui-macro-sens-dist\" %s]
" % gui_macro_sens_dist)
			file_write.write("[\"gui-nlc-threshold\" %s]
" % gui_nlc_threshold)
			file_write.write("[\"gui-sig-init\" %s]
" % gui_sig_init)
			file_write.write("[\"gui-sig-init-std\" %s]
" % gui_sig_init_std)
			file_write.write("[\"gui-diff-mean\" %s]
" % gui_diff_mean)
			file_write.write("[\"gui-diff-std\" %s]
" % gui_diff_std)
			file_write.write("[\"gui-life-init-gamma\" %s]
" % gui_life_init_gamma)
			file_write.write("[\"gui-alpha-distrib\" %s]

" % gui_alpha_distrib)

			n = n+1

print(n)
stop = timeit.default_timer()
print(stop - start)  

parse_best_param_for_behavior_space_adapted.py 
#!/usr/bin/python
# -*- coding: utf-8 -*-
# jupyter notebook --browser='C:/Program Files (x86)/Google/Chrome/Application/chrome.exe %s'


import timeit
start = timeit.default_timer()

import sys
import os

my_file = sys.argv[1]
output = sys.argv[2]

patient = my_file.split("/")[0]
with open("patient_dict.txt", 'r') as file_read :
	data = file_read.readlines()
	for line in data[1:] :
		line = line.replace("
", "").split(" ")
		patient_name = line[0]
		if patient_name == patient :
			mono_init = line[1]
			apo_init = line[2]

os.mkdir('%s/BehaviorSpace' % patient)

with open(my_file, 'r') as file_read :
	with open("%s" % output, 'w') as file_write :
		data = file_read.readlines()
		n = 1
		for line in data[1:] :
			line = line.replace("
", "").split("	")
			gui_apo_mov = line[3]
			gui_need_sig_mov = line[4]
			gui_layers = line[5]
			gui_alpha = line[6]
			gui_mono_phago_eff = line[7]
			gui_NLC_phago_eff = line[8]
			gui_M_phago_eff = line[9]
			gui_M_kill_eff = line[10]
			gui_cll_sens_dist = line[11]
			gui_mono_sens_dist = line[12]
			gui_nlc_sens_dist = line[13]
			gui_macro_sens_dist = line[14]
			gui_nlc_threshold = line[15]
			gui_sig_init = line[16]
			gui_sig_init_std = line[17]
			gui_diff_mean = line[18]
			gui_diff_std = line[19]
			gui_life_init_gamma = line[20]
			gui_alpha_distrib = line[21]
			if (n == 1) :
				file_write.write("Best_via_set
")
			elif (n == 2) :
				file_write.write("Knee_point_set
")
			elif (n == 3) :
				file_write.write("Best_conc_set
")
			file_write.write("[\"gui-prop-mono-init\" %s]
" % mono_init)
			file_write.write("[\"gui-prop-apo-init\" %s]
" % apo_init)				
			file_write.write("[\"gui-apo-mov\" %s]
" % gui_apo_mov)
			file_write.write("[\"gui-need-sig-mov\" %s]
" % gui_need_sig_mov)
			file_write.write("[\"gui-layers\" %s]
" % gui_layers)
			file_write.write("[\"gui-alpha\" %s]
" % gui_alpha)
			file_write.write("[\"gui-mono-phago-eff\" %s]
" % gui_mono_phago_eff)
			file_write.write("[\"gui-NLC-phago-eff\" %s]
" % gui_NLC_phago_eff)
			file_write.write("[\"gui-M-phago-eff\" %s]
" % gui_M_phago_eff)
			file_write.write("[\"gui-M-kill-eff\" %s]
" % gui_M_kill_eff)
			file_write.write("[\"gui-cll-sens-dist\" %s]
" % gui_cll_sens_dist)
			file_write.write("[\"gui-mono-sens-dist\" %s]
" % gui_mono_sens_dist)
			file_write.write("[\"gui-nlc-sens-dist\" %s]
" % gui_nlc_sens_dist)
			file_write.write("[\"gui-macro-sens-dist\" %s]
" % gui_macro_sens_dist)
			file_write.write("[\"gui-nlc-threshold\" %s]
" % gui_nlc_threshold)
			file_write.write("[\"gui-sig-init\" %s]
" % gui_sig_init)
			file_write.write("[\"gui-sig-init-std\" %s]
" % gui_sig_init_std)
			file_write.write("[\"gui-diff-mean\" %s]
" % gui_diff_mean)
			file_write.write("[\"gui-diff-std\" %s]
" % gui_diff_std)
			file_write.write("[\"gui-life-init-gamma\" %s]
" % gui_life_init_gamma)
			file_write.write("[\"gui-alpha-distrib\" %s]

" % gui_alpha_distrib)

			n = n+1

print(n)
stop = timeit.default_timer()
print(stop - start)  

pca_analysis.py 
#!/usr/bin/python
# -*- coding: utf-8 -*-
# jupyter notebook --browser='C:/Program Files (x86)/Google/Chrome/Application/chrome.exe %s'


import timeit
start = timeit.default_timer()

import sys
import matplotlib.pyplot as plt

# import os

#### usage : python scripts/pca_analysis.py best_via_sets_all_patients.tsv

my_file = sys.argv[1] # best_via_sets_all_patients.tsv
# output = sys.argv[2]

import pandas as pd

df = pd.read_csv(my_file, sep='	', header=0, index_col=0)
print(df)

my_data = df.T
my_data.reset_index(inplace=True)
my_data.rename(columns = {'index':'patient'}, inplace = True)
my_data.drop(["delta_fitness_via", "delta_fitness_conc"], axis = 1, inplace = True)

print(my_data)


from sklearn.preprocessing import StandardScaler
my_features = list(my_data)
my_features = my_features[1:]
print(my_features)

features = my_features

# features = []
# for feature in my_features :
# 	feature = feature.replace(" ", "")
# 	features.append(feature)

# Separating out the features
x = my_data.loc[:, features].values

# Separating out the target
y = my_data.loc[:,['patient']].values

# Standardizing the features
x = StandardScaler().fit_transform(x)


from sklearn.decomposition import PCA
pca = PCA(n_components=2)
principalComponents = pca.fit_transform(x)
principalDf = pd.DataFrame(data = principalComponents, columns = ['principal component 1', 'principal component 2'])


finalDf = pd.concat([principalDf, my_data[['patient']]], axis = 1)




fig = plt.figure(figsize = (8,8))
ax = fig.add_subplot(1,1,1) 
ax.set_xlabel('Principal Component 1', fontsize = 15)
ax.set_ylabel('Principal Component 2', fontsize = 15)


patients = my_data.patient.values.tolist()

colors = ['#00FFFF','#458B74','#000000','#0000FF','#8A2BE2','#FF4040', '#0FE721', '#FF0000', '#FC00FF', '#00FFB6', '#8B2F00']

for patient, color in zip(patients,colors):
    indicesToKeep = finalDf['patient'] == patient
    ax.scatter(finalDf.loc[indicesToKeep, 'principal component 1']
               , finalDf.loc[indicesToKeep, 'principal component 2']
               , c = color
               , s = 50)
ax.legend(patients)
ax.grid()

# plt.show()
ax.set_title('2-component PCA - %s
Explained Variance = %s' % (my_file[:-17],pca.explained_variance_ratio_), fontsize = 20)

plt.savefig('pca_analysis_%s.pdf' % my_file[:-4])
plt.savefig('pca_analysis_%s.png' % my_file[:-4])

stop = timeit.default_timer()
print(stop - start)  

plot_kneepoint0_sim_vs_exp_with_scores.py 
#!/usr/bin/python
# -*- coding: utf-8 -*-
# jupyter notebook --browser='C:/Program Files (x86)/Google/Chrome/Application/chrome.exe %s'


import timeit
start = timeit.default_timer()

import sys
import numpy
import pandas as pd
import matplotlib.pyplot as plt

## usage : $ python scripts/plot_sim_vs_exp_with_scores.py CAS1802 best_via

## plot simulations vs experimental data

patient = sys.argv[1]
param_set = sys.argv[2]

simu_file_path = "%s\BehaviorSpace\%s.csv" % (patient,param_set)
# print('simu_file_path', simu_file_path)

## get patient experimental data
patient_data = pd.read_csv('%s\%s.csv' % (patient,patient)) #, index_col=0) 
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
axes[0].text(0.1, 0.1, 'NRMSE_via = %s
 R²_via = %s' % (nrms_via,squared_r_via), horizontalalignment='left', verticalalignment='baseline', transform=axes[0].transAxes)
axes[1].text(0.1, 0.1, 'NRMSE_conc = %s
 R²_conc = %s' % (nrms_conc,squared_r_conc), horizontalalignment='left', verticalalignment='baseline', transform=axes[1].transAxes)
# plt.show()

# plt.savefig('%s\BehaviorSpace\%s_%s_model_fit_with_scores.pdf' % (patient, patient, param_set))
# plt.savefig('%s\BehaviorSpace\%s_%s_model_fit_with_scores.png' % (patient, patient, param_set))


plt.savefig('plots_with_scores_kneepoint0\%s_%s_model_fit_with_scores.pdf' % (patient, param_set))
plt.savefig('plots_with_scores_kneepoint0\%s_%s_model_fit_with_scores.png' % (patient, param_set))

stop = timeit.default_timer()
print(stop - start)  


# ### get simulation data
# sim_data = pd.read_csv(simu_file_path, skiprows=6, sep=",", header=0)
# print(sim_data['[run number]'])

# runs = sim_data['[run number]'].unique()
# print(sorted(runs))

# for run in runs :
# 	sim_viability = pd.DataFrame()

plot_kneepoint1_sim_vs_exp_with_scores.py 
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

simu_file_path = "%s\BehaviorSpace\%s.csv" % (patient,param_set)
# print('simu_file_path', simu_file_path)

## get patient experimental data
patient_data = pd.read_csv('%s\%s.csv' % (patient,patient)) #, index_col=0) 
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
axes[0].text(0.1, 0.1, 'NRMSE_via = %s
 R²_via = %s' % (nrms_via,squared_r_via), horizontalalignment='left', verticalalignment='baseline', transform=axes[0].transAxes)
axes[1].text(0.1, 0.1, 'NRMSE_conc = %s
 R²_conc = %s' % (nrms_conc,squared_r_conc), horizontalalignment='left', verticalalignment='baseline', transform=axes[1].transAxes)
# plt.show()

# plt.savefig('%s\BehaviorSpace\%s_%s_model_fit_with_scores.pdf' % (patient, patient, param_set))
# plt.savefig('%s\BehaviorSpace\%s_%s_model_fit_with_scores.png' % (patient, patient, param_set))


plt.savefig('plots_with_scores_kneepoint1\%s_%s_model_fit_with_scores.pdf' % (patient, param_set))
plt.savefig('plots_with_scores_kneepoint1\%s_%s_model_fit_with_scores.png' % (patient, param_set))

stop = timeit.default_timer()
print(stop - start)  


# ### get simulation data
# sim_data = pd.read_csv(simu_file_path, skiprows=6, sep=",", header=0)
# print(sim_data['[run number]'])

# runs = sim_data['[run number]'].unique()
# print(sorted(runs))

# for run in runs :
# 	sim_viability = pd.DataFrame()

plot_kneepoint2_sim_vs_exp_with_scores.py 
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

simu_file_path = "%s\BehaviorSpace\%s.csv" % (patient,param_set)
# print('simu_file_path', simu_file_path)

## get patient experimental data
patient_data = pd.read_csv('%s\%s.csv' % (patient,patient)) #, index_col=0) 
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
axes[0].text(0.1, 0.1, 'NRMSE_via = %s
 R²_via = %s' % (nrms_via,squared_r_via), horizontalalignment='left', verticalalignment='baseline', transform=axes[0].transAxes)
axes[1].text(0.1, 0.1, 'NRMSE_conc = %s
 R²_conc = %s' % (nrms_conc,squared_r_conc), horizontalalignment='left', verticalalignment='baseline', transform=axes[1].transAxes)
# plt.show()

# plt.savefig('%s\BehaviorSpace\%s_%s_model_fit_with_scores.pdf' % (patient, patient, param_set))
# plt.savefig('%s\BehaviorSpace\%s_%s_model_fit_with_scores.png' % (patient, patient, param_set))


plt.savefig('plots_with_scores_kneepoint2\%s_%s_model_fit_with_scores.pdf' % (patient, param_set))
plt.savefig('plots_with_scores_kneepoint2\%s_%s_model_fit_with_scores.png' % (patient, param_set))

stop = timeit.default_timer()
print(stop - start)  


# ### get simulation data
# sim_data = pd.read_csv(simu_file_path, skiprows=6, sep=",", header=0)
# print(sim_data['[run number]'])

# runs = sim_data['[run number]'].unique()
# print(sorted(runs))

# for run in runs :
# 	sim_viability = pd.DataFrame()

plot_predictions_model0.py 
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

df_via_exp = pd.read_csv('allPatients_prediction-via-exp.tsv', sep='	')
df_via_exp['viability_mean'] = df_via_exp.mean(axis=1)
df_via_exp['viability_std'] = df_via_exp.std(axis=1)
df_via_exp = df_via_exp.iloc[::-1]

df_conc_exp = pd.read_csv('allPatients_prediction-conc-exp.tsv', sep='	')
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
# fig.suptitle('B-CLL viability and concentration
 with varying monocytes initial proportions (3 patients)
 Experimental vs. Predictions', fontsize=10, y=1.02)
# fig.suptitle('Model Predictions 
 %s' % file, fontsize=10, y=1.02)
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
# Normalized RMSE = RMSE / (max value – min value)
sum_via_errors = sum(via_errors_pred2)
RMSE_via = math.sqrt(sum_via_errors / len(x_labels))
# print("RMSE Via for %s = %s" % (file,RMSE_via))
normalized_RMSE_via = RMSE_via / (max(exp_via_means_list) - min(exp_via_means_list))
print("Normalized RMSE Via for %s = %s" % ("Model 0",normalized_RMSE_via))

## calculate the R² for viability
from sklearn.metrics import r2_score
### Assume y is the actual value and f is the predicted values
# r2 = r2_score(y, f)
r2_via = r2_score(exp_via_means_list, exp_via_means_list)
print("R² Via = %s" % r2_via)

## calculate the RMSE for concentration
sum_conc_errors = sum(conc_errors_pred2)
RMSE_conc = math.sqrt(sum_conc_errors / len(x_labels))
# print("RMSE Conc for %s = %s" % (file,RMSE_conc))
normalized_RMSE_conc = RMSE_conc / (max(exp_conc_means_list) - min(exp_conc_means_list))
print("Normalized RMSE Conc for %s = %s" % ("Model 0",normalized_RMSE_conc))

## calculate the R² for concentration
# r2 = r2_score(y, f)
r2_conc = r2_score(exp_conc_means_list, concentration_mean)
print("R² Conc = %s" % r2_conc)

### round 2 decimals
normalized_RMSE_via = round(normalized_RMSE_via, 2)
normalized_RMSE_conc = round(normalized_RMSE_conc, 2)
r2_via = round(r2_via, 2)
r2_conc = round(r2_conc, 2)

# Show NRMSE and squared-R results
# (0, 0) being the lower left of the axes and (1, 1) the upper right.
axs[0].text(0.1, 0.775, 'NRMSE_via = %s
 R²_via = %s' % (normalized_RMSE_via,r2_via), horizontalalignment='left', verticalalignment='baseline', transform=axs[0].transAxes)
axs[1].text(0.4, 0.775, 'NRMSE_conc = %s
 R²_conc = %s' % (normalized_RMSE_conc,r2_conc), horizontalalignment='left', verticalalignment='baseline', transform=axs[1].transAxes)

# set_matplotlib_formats('pdf')
# # # set_matplotlib_formats('svg')
# set_matplotlib_formats('png')

# plt.show()

# plt.savefig("pred_Model0.png")
plt.savefig("pred_Model0.svg")
# files.download("pred_%s.pdf" % patient) 



stop = timeit.default_timer()
print(stop - start)  

plot_sensitivity_analysis.py 
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
simu_file_name = "ABM_2D_9patients_0_%s" % sys.argv[1]

param = ("-").join(sys.argv[1].split("-")[1:])

exp_file_path = "A:\Downloads\Projects\workFromHome\Projects\ABM2021\paper\revision\9patients"

## get patients experimental data
viability_exp = pd.read_csv('%s\filtered_fused_9patients_viability.tsv' % exp_file_path, sep='	') 
concentration_exp = pd.read_csv('%s\filtered_fused_9patients_concentration.tsv' % exp_file_path, sep='	') 

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
fig.suptitle('ABM_2D_9patients_0_%s' % sys.argv[1], fontsize=14)

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


plt.savefig("%s.pdf" % sys.argv[1])
plt.savefig("%s.png" % sys.argv[1])


stop = timeit.default_timer()
print(stop - start)

plot_sensitivity_analysis_class1.py 
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

exp_file_path = "A:\Downloads\Projects\workFromHome\Projects\ABM2021\paper\revision\classes_ABM_2D_9patients_1"

## get patients experimental data
viability_exp = pd.read_csv('%s\class1_viability.tsv' % exp_file_path, sep='	') 
concentration_exp = pd.read_csv('%s\class1_concentration.tsv' % exp_file_path, sep='	') 

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


plt.savefig("sensitivity_analysis_class1\%s_class1.pdf" % sys.argv[1])
plt.savefig("sensitivity_analysis_class1\%s_class1.png" % sys.argv[1])


stop = timeit.default_timer()
print(stop - start)

plot_sensitivity_analysis_class2.py 
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
simu_file_name = "ABM_2D_9patients_1_class2_%s" % sys.argv[1]

param = ("-").join(sys.argv[1].split("-")[1:])

exp_file_path = "A:\Downloads\Projects\workFromHome\Projects\ABM2021\paper\revision\classes_ABM_2D_9patients_1"

## get patients experimental data
viability_exp = pd.read_csv('%s\class2_viability.tsv' % exp_file_path, sep='	') 
concentration_exp = pd.read_csv('%s\class2_concentration.tsv' % exp_file_path, sep='	') 

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
fig.suptitle('ABM_2D_Class2_%s' % sys.argv[1], fontsize=14)

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


plt.savefig("sensitivity_analysis_class2\%s_class2.pdf" % sys.argv[1])
plt.savefig("sensitivity_analysis_class2\%s_class2.png" % sys.argv[1])


stop = timeit.default_timer()
print(stop - start)

plot_sim_vs_exp.py 
#!/usr/bin/python
# -*- coding: utf-8 -*-
# jupyter notebook --browser='C:/Program Files (x86)/Google/Chrome/Application/chrome.exe %s'


import timeit
start = timeit.default_timer()

import sys
import numpy
import pandas as pd
import matplotlib.pyplot as plt

## usage : $ python scripts/plot_sim_vs_exp.py CAS1802 best_via

## plot simulations vs experimental data

patient = sys.argv[1]
param_set = sys.argv[2]

simu_file_path = "%s\BehaviorSpace\%s.csv" % (patient,param_set)
# print('simu_file_path', simu_file_path)

## get patient experimental data
patient_data = pd.read_csv('%s\%s.csv' % (patient,patient)) #, index_col=0) 
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
fig.suptitle('Model fitting for %s_%s' % (patient, param_set), fontsize=12)

filtered_df_viability_simu.plot(ax=axes[0], legend=False, color='red',zorder=0)#, ylim=(60,100))
filtered_df_remaining_simu.plot(ax=axes[1], legend=False, color='red',zorder=0)#, ylim=(60,100))

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
axes[0].set_xticks(range(0,24*14, 24))
axes[0].set_xticklabels(range(0,14))

axes[1].set_xlabel('Day')
axes[1].set_ylabel('Concentration Ratio (%)')
axes[1].set_xticks(range(0,24*14, 24))
axes[1].set_xticklabels(range(0,14))


leg = axes[1].legend(['Simulation', 'Exp'])#, loc='upper right', bbox_to_anchor=(1.3,1.3), )


LH = leg.legendHandles
LH[0].set_linewidth(8)
LH[1].set_color('black') 

# Show NRMSE and squared-R results
axes[0].text(0.1, 0.1, 'NRMSE_via = %s
 R²_via = %s' % (nrms_via,squared_r_via), horizontalalignment='left', verticalalignment='baseline', transform=axes[0].transAxes)
axes[1].text(0.1, 0.1, 'NRMSE_conc = %s
 R²_conc = %s' % (nrms_conc,squared_r_conc), horizontalalignment='left', verticalalignment='baseline', transform=axes[1].transAxes)
plt.show()

# plt.savefig('%s\BehaviorSpace\%s_%s_model_fit.pdf' % (patient, patient, param_set))
# plt.savefig('%s\BehaviorSpace\%s_%s_model_fit.png' % (patient, patient, param_set))

stop = timeit.default_timer()
print(stop - start)  


# ### get simulation data
# sim_data = pd.read_csv(simu_file_path, skiprows=6, sep=",", header=0)
# print(sim_data['[run number]'])

# runs = sim_data['[run number]'].unique()
# print(sorted(runs))

# for run in runs :
# 	sim_viability = pd.DataFrame()

plot_sim_vs_exp_class1.py 
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

exp_file_path = "A:\Downloads\Projects\workFromHome\Projects\ABM2021\paper\revision\classes"
# simu_file_path = "A:\Downloads\Projects\workFromHome\Projects\ABM2021\20212201\figures\plots\BehaviorSpace\NLC-CLL-revisions\10patients"
simu_file_name = sys.argv[1]

## get patients experimental data
viability_exp = pd.read_csv('%s\class1_viability.tsv' % exp_file_path, sep='	') 
concentration_exp = pd.read_csv('%s\class1_concentration.tsv' % exp_file_path, sep='	') 

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

fig, axes = plt.subplots(nrows=5, ncols=1)
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
# print(my_list)

df_viability_simu = pd.DataFrame.from_dict(viability_dict)
df_remaining_simu = pd.DataFrame.from_dict(remaining_dict)

filtered_df_viability_simu = df_viability_simu[df_viability_simu.index.isin(my_list)]
filtered_df_remaining_simu = df_remaining_simu[df_remaining_simu.index.isin(my_list)]

# print(filtered_df_viability_simu)


filtered_df_viability_simu.plot(ax=axes[0],legend=False, color='r', zorder=0)#, ylim=(00,100), color='orange')
x_axis = axes[0].axes.get_xaxis()
x_label = x_axis.get_label()
x_label.set_visible(False)

axes[0].set_ylabel('Viability (%)')
# axes[0].set_xticks(range(0,25*14, 24))
axes[0].set_xticks(range(0,24*14, 24))
axes[0].set_xticklabels(range(0,14))

filtered_df_remaining_simu.plot(figsize=(8, 8),ax=axes[1], legend=False, color='r', zorder=0)#, ylim=(00,110), color='orange')
axes[1].set_xlabel('Day')
axes[1].set_ylabel('Concentration ratio (%)')
axes[1].set_xticks(range(0,24*14, 24))
axes[1].set_xticklabels(range(0,14))



# leg = axes[1].legend(['Simulation', 'Exp'], loc='upper right', bbox_to_anchor=(1.1,1.1), )


# LH = leg.legendHandles
# LH[-1].set_linewidth(20)
# LH[1].set_linewidth(2)
# LH[1].set_color('black') 
# LH[-1].set_color('red') 



### get optimized parameters
params = pd.read_csv("%s.csv" % (simu_file_name), skiprows=6, sep=",", header=0, nrows = 1, usecols=[3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21])
print(params)

params['gui-alpha-distrib'] = params['gui-alpha-distrib'].round(decimals = 3)


# params = params.T
params_firstrow = params[["gui-apo-mov", "gui-need-sig-mov", "gui-cll-sens-dist", "gui-mono-sens-dist", "gui-macro-sens-dist", "gui-nlc-sens-dist"]]
params_secondrow = params[["gui-life-init-gamma","gui-alpha-distrib",  "gui-diff-mean", "gui-diff-std", "gui-sig-init", "gui-sig-init-std", "gui-nlc-threshold"]]
params_thirdrow = params[["gui-layers", "gui-alpha", "gui-mono-phago-eff", "gui-M-phago-eff", "gui-NLC-phago-eff", "gui-M-kill-eff"]]


# hide axes
fig.patch.set_visible(False)
axes[2].axis('off')
# axes[2].axis('tight')

axes[3].axis('off')
# axes[3].axis('tight')

axes[4].axis('off')
# axes[4].axis('tight')

table1 = axes[2].table(cellText=params_firstrow.values, colLabels=[i[4:] for i in params_firstrow.columns], loc='center', fontsize=14).scale(xscale=1,yscale=2)

table2 = axes[3].table(cellText=params_secondrow.values, colLabels=[i[4:] for i in params_secondrow.columns], loc='center', fontsize=14).scale(xscale=1,yscale=2)
table3 = axes[4].table(cellText=params_thirdrow.values, colLabels=[i[4:] for i in params_thirdrow.columns], loc='center', fontsize=14).scale(xscale=1,yscale=2)


# axes[2].table.auto_set_font_size(False)
# table2.auto_set_font_size(False)
# table3.auto_set_font_size(False)

# plt.subplots_adjust(left=0.125,
#                     bottom=0.1, 
#                     right=0.9, 
#                     top=0.9, 
#                     wspace=0.1, 
#                     hspace=0.0)
# fig.tight_layout()

plt.subplots_adjust(top=0.95)

plt.subplots_adjust(hspace=0.0,wspace=0.1)



# plt.show()


# plt.show()

plt.savefig('%s_model_fit.pdf' % (simu_file_name))
plt.savefig('%s_model_fit.png' % (simu_file_name))


stop = timeit.default_timer()
print(stop - start)  

plot_sim_vs_exp_class1_with_scores.py 
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

exp_file_path = "A:\Downloads\Projects\workFromHome\Projects\ABM2021\paper\revision\
plot_sim_vs_exp_class2.py 
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

exp_file_path = "A:\Downloads\Projects\workFromHome\Projects\ABM2021\paper\revision\classes"
# simu_file_path = "A:\Downloads\Projects\workFromHome\Projects\ABM2021\20212201\figures\plots\BehaviorSpace\NLC-CLL-revisions\10patients"
simu_file_name = sys.argv[1]

## get patients experimental data
viability_exp = pd.read_csv('%s\class2_viability.tsv' % exp_file_path, sep='	') 
concentration_exp = pd.read_csv('%s\class2_concentration.tsv' % exp_file_path, sep='	') 

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

fig, axes = plt.subplots(nrows=5, ncols=1)
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
# print(my_list)

df_viability_simu = pd.DataFrame.from_dict(viability_dict)
df_remaining_simu = pd.DataFrame.from_dict(remaining_dict)

filtered_df_viability_simu = df_viability_simu[df_viability_simu.index.isin(my_list)]
filtered_df_remaining_simu = df_remaining_simu[df_remaining_simu.index.isin(my_list)]

# print(filtered_df_viability_simu)


filtered_df_viability_simu.plot(ax=axes[0],legend=False, color='r', zorder=0)#, ylim=(00,100), color='orange')
x_axis = axes[0].axes.get_xaxis()
x_label = x_axis.get_label()
x_label.set_visible(False)

axes[0].set_ylabel('Viability (%)')
# axes[0].set_xticks(range(0,25*14, 24))
axes[0].set_xticks(range(0,24*14, 24))
axes[0].set_xticklabels(range(0,14))

filtered_df_remaining_simu.plot(figsize=(8, 8),ax=axes[1], legend=False, color='r', zorder=0)#, ylim=(00,110), color='orange')
axes[1].set_xlabel('Day')
axes[1].set_ylabel('Concentration ratio (%)')
axes[1].set_xticks(range(0,24*14, 24))
axes[1].set_xticklabels(range(0,14))



# leg = axes[1].legend(['Simulation', 'Exp'], loc='upper right', bbox_to_anchor=(1.1,1.1), )


# LH = leg.legendHandles
# LH[-1].set_linewidth(20)
# LH[1].set_linewidth(2)
# LH[1].set_color('black') 
# LH[-1].set_color('red') 



### get optimized parameters
params = pd.read_csv("%s.csv" % (simu_file_name), skiprows=6, sep=",", header=0, nrows = 1, usecols=[3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21])
print(params)

params['gui-alpha-distrib'] = params['gui-alpha-distrib'].round(decimals = 3)


# params = params.T
params_firstrow = params[["gui-apo-mov", "gui-need-sig-mov", "gui-cll-sens-dist", "gui-mono-sens-dist", "gui-macro-sens-dist", "gui-nlc-sens-dist"]]
params_secondrow = params[["gui-life-init-gamma","gui-alpha-distrib",  "gui-diff-mean", "gui-diff-std", "gui-sig-init", "gui-sig-init-std", "gui-nlc-threshold"]]
params_thirdrow = params[["gui-layers", "gui-alpha", "gui-mono-phago-eff", "gui-M-phago-eff", "gui-NLC-phago-eff", "gui-M-kill-eff"]]


# hide axes
fig.patch.set_visible(False)
axes[2].axis('off')
# axes[2].axis('tight')

axes[3].axis('off')
# axes[3].axis('tight')

axes[4].axis('off')
# axes[4].axis('tight')

table1 = axes[2].table(cellText=params_firstrow.values, colLabels=[i[4:] for i in params_firstrow.columns], loc='center', fontsize=14).scale(xscale=1,yscale=2)

table2 = axes[3].table(cellText=params_secondrow.values, colLabels=[i[4:] for i in params_secondrow.columns], loc='center', fontsize=14).scale(xscale=1,yscale=2)
table3 = axes[4].table(cellText=params_thirdrow.values, colLabels=[i[4:] for i in params_thirdrow.columns], loc='center', fontsize=14).scale(xscale=1,yscale=2)


# axes[2].table.auto_set_font_size(False)
# table2.auto_set_font_size(False)
# table3.auto_set_font_size(False)

# plt.subplots_adjust(left=0.125,
#                     bottom=0.1, 
#                     right=0.9, 
#                     top=0.9, 
#                     wspace=0.1, 
#                     hspace=0.0)
# fig.tight_layout()

plt.subplots_adjust(top=0.95)

plt.subplots_adjust(hspace=0.0,wspace=0.1)



# plt.show()


# plt.show()

plt.savefig('%s_model_fit.pdf' % (simu_file_name))
plt.savefig('%s_model_fit.png' % (simu_file_name))


stop = timeit.default_timer()
print(stop - start)  

plot_sim_vs_exp_class2_with_scores.py 
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

exp_file_path = "A:\Downloads\Projects\workFromHome\Projects\ABM2021\paper\revision\
plot_sim_vs_exp_with_scores.py 
#!/usr/bin/python
# -*- coding: utf-8 -*-
# jupyter notebook --browser='C:/Program Files (x86)/Google/Chrome/Application/chrome.exe %s'


import timeit
start = timeit.default_timer()

import sys
import numpy
import pandas as pd
import matplotlib.pyplot as plt

## usage : $ python scripts/plot_sim_vs_exp_with_scores.py CAS1802 best_via

## plot simulations vs experimental data

patient = sys.argv[1]
param_set = sys.argv[2]

simu_file_path = "%s\BehaviorSpace\%s.csv" % (patient,param_set)
# print('simu_file_path', simu_file_path)

## get patient experimental data
patient_data = pd.read_csv('%s\%s.csv' % (patient,patient)) #, index_col=0) 
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
fig.suptitle('Model fitting for %s_%s' % (patient, param_set), fontsize=12)

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
axes[0].text(0.1, 0.1, 'NRMSE_via = %s
 R²_via = %s' % (nrms_via,squared_r_via), horizontalalignment='left', verticalalignment='baseline', transform=axes[0].transAxes)
axes[1].text(0.1, 0.1, 'NRMSE_conc = %s
 R²_conc = %s' % (nrms_conc,squared_r_conc), horizontalalignment='left', verticalalignment='baseline', transform=axes[1].transAxes)
# plt.show()

# plt.savefig('%s\BehaviorSpace\%s_%s_model_fit_with_scores.pdf' % (patient, patient, param_set))
# plt.savefig('%s\BehaviorSpace\%s_%s_model_fit_with_scores.png' % (patient, patient, param_set))


plt.savefig('plots_with_scores_byclass\%s_%s_model_fit_with_scores.pdf' % (patient, param_set))
plt.savefig('plots_with_scores_byclass\%s_%s_model_fit_with_scores.png' % (patient, param_set))

stop = timeit.default_timer()
print(stop - start)  


# ### get simulation data
# sim_data = pd.read_csv(simu_file_path, skiprows=6, sep=",", header=0)
# print(sim_data['[run number]'])

# runs = sim_data['[run number]'].unique()
# print(sorted(runs))

# for run in runs :
# 	sim_viability = pd.DataFrame()

plot_sim_vs_exp_with_scores_Model2.py 
#!/usr/bin/python
# -*- coding: utf-8 -*-
# jupyter notebook --browser='C:/Program Files (x86)/Google/Chrome/Application/chrome.exe %s'


import timeit
start = timeit.default_timer()

import sys
import numpy
import pandas as pd
import matplotlib.pyplot as plt

## usage : $ python plot_sim_vs_exp_with_scores_Model2.py Model2_stocha

exp_file_path = "A:\Downloads\Projects\workFromHome\Projects\ABM2021\paper\revision\classes_ABM_2D_9patients_1"
# simu_file_path = "A:\Downloads\Projects\workFromHome\Projects\ABM2021\20212201\figures\plots\BehaviorSpace\NLC-CLL-revisions\10patients"
simu_file_name = sys.argv[1]

## get patients experimental data
viability_exp = pd.read_csv('%s\class2_viability.tsv' % exp_file_path, sep='	') 
concentration_exp = pd.read_csv('%s\class2_concentration.tsv' % exp_file_path, sep='	') 

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
axes[0].text(0.1, 0.1, 'NRMSE_via = %s
 R²_via = %s' % (nrms_via,squared_r_via), horizontalalignment='left', verticalalignment='baseline', transform=axes[0].transAxes)
axes[1].text(0.1, 0.1, 'NRMSE_conc = %s
 R²_conc = %s' % (nrms_conc,squared_r_conc), horizontalalignment='left', verticalalignment='baseline', transform=axes[1].transAxes)


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

plot_sim_vs_exp_with_scores_class1_averaged.py 
#!/usr/bin/python
# -*- coding: utf-8 -*-
# jupyter notebook --browser='C:/Program Files (x86)/Google/Chrome/Application/chrome.exe %s'


import timeit
start = timeit.default_timer()

import sys
import numpy
import pandas as pd
import matplotlib.pyplot as plt

## usage : $ python plot_sim_vs_exp_with_scores_class1_averaged.py class1_averaged_stocha

exp_file_path = "A:\Downloads\Projects\workFromHome\Projects\ABM2021\paper\revision\classes_ABM_2D_9patients_1"
# simu_file_path = "A:\Downloads\Projects\workFromHome\Projects\ABM2021\20212201\figures\plots\BehaviorSpace\NLC-CLL-revisions\10patients"
simu_file_name = sys.argv[1]

## get patients experimental data
viability_exp = pd.read_csv('%s\class1_viability.tsv' % exp_file_path, sep='	') 
concentration_exp = pd.read_csv('%s\class1_concentration.tsv' % exp_file_path, sep='	') 

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
axes[0].text(0.1, 0.1, 'NRMSE_via = %s
 R²_via = %s' % (nrms_via,squared_r_via), horizontalalignment='left', verticalalignment='baseline', transform=axes[0].transAxes)
axes[1].text(0.1, 0.1, 'NRMSE_conc = %s
 R²_conc = %s' % (nrms_conc,squared_r_conc), horizontalalignment='left', verticalalignment='baseline', transform=axes[1].transAxes)


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

plot_sim_vs_exp_with_scores_class2_averaged.py 
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

exp_file_path = "A:\Downloads\Projects\workFromHome\Projects\ABM2021\paper\revision\classes_ABM_2D_9patients_1"
# simu_file_path = "A:\Downloads\Projects\workFromHome\Projects\ABM2021\20212201\figures\plots\BehaviorSpace\NLC-CLL-revisions\10patients"
simu_file_name = sys.argv[1]

## get patients experimental data
viability_exp = pd.read_csv('%s\class2_viability.tsv' % exp_file_path, sep='	') 
concentration_exp = pd.read_csv('%s\class2_concentration.tsv' % exp_file_path, sep='	') 

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
axes[0].text(0.1, 0.1, 'NRMSE_via = %s
 R²_via = %s' % (nrms_via,squared_r_via), horizontalalignment='left', verticalalignment='baseline', transform=axes[0].transAxes)
axes[1].text(0.1, 0.1, 'NRMSE_conc = %s
 R²_conc = %s' % (nrms_conc,squared_r_conc), horizontalalignment='left', verticalalignment='baseline', transform=axes[1].transAxes)


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

plot_sim_vs_exp_with_scores_model_1.py 
#!/usr/bin/python
# -*- coding: utf-8 -*-
# jupyter notebook --browser='C:/Program Files (x86)/Google/Chrome/Application/chrome.exe %s'


import timeit
start = timeit.default_timer()

import sys
import numpy
import pandas as pd
import matplotlib.pyplot as plt

## usage : $ python plot_sim_vs_exp_with_scores_model_1.py Model1_stocha

exp_file_path = "A:\Downloads\Projects\workFromHome\Projects\ABM2021\paper\revision\classes_ABM_2D_9patients_1"
# simu_file_path = "A:\Downloads\Projects\workFromHome\Projects\ABM2021\20212201\figures\plots\BehaviorSpace\NLC-CLL-revisions\10patients"
simu_file_name = sys.argv[1]

## get patients experimental data
viability_exp = pd.read_csv('%s\class1_viability.tsv' % exp_file_path, sep='	') 
concentration_exp = pd.read_csv('%s\class1_concentration.tsv' % exp_file_path, sep='	') 

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
axes[0].text(0.1, 0.1, 'NRMSE_via = %s
 R²_via = %s' % (nrms_via,squared_r_via), horizontalalignment='left', verticalalignment='baseline', transform=axes[0].transAxes)
axes[1].text(0.1, 0.1, 'NRMSE_conc = %s
 R²_conc = %s' % (nrms_conc,squared_r_conc), horizontalalignment='left', verticalalignment='baseline', transform=axes[1].transAxes)


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

plot_sim_vs_exp_with_scores_only_kneepoint_patient_specific.py 
#!/usr/bin/python
# -*- coding: utf-8 -*-
# jupyter notebook --browser='C:/Program Files (x86)/Google/Chrome/Application/chrome.exe %s'


import timeit
start = timeit.default_timer()

import sys
import numpy
import pandas as pd
import matplotlib.pyplot as plt

## usage : $ python scripts/plot_sim_vs_exp_with_scores.py CAS1802 best_via

## plot simulations vs experimental data

patient = sys.argv[1]
param_set = sys.argv[2]

simu_file_path = "%s\BehaviorSpace\%s.csv" % (patient,param_set)
# print('simu_file_path', simu_file_path)

## get patient experimental data
patient_data = pd.read_csv('%s\%s.csv' % (patient,patient)) #, index_col=0) 
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
fig.suptitle('Model fitting for %s_%s' % (patient, param_set), fontsize=12)

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
axes[0].text(0.1, 0.1, 'NRMSE_via = %s
 R²_via = %s' % (nrms_via,squared_r_via), horizontalalignment='left', verticalalignment='baseline', transform=axes[0].transAxes)
axes[1].text(0.1, 0.1, 'NRMSE_conc = %s
 R²_conc = %s' % (nrms_conc,squared_r_conc), horizontalalignment='left', verticalalignment='baseline', transform=axes[1].transAxes)
# plt.show()

# plt.savefig('%s\BehaviorSpace\%s_%s_model_fit_with_scores.pdf' % (patient, patient, param_set))
# plt.savefig('%s\BehaviorSpace\%s_%s_model_fit_with_scores.png' % (patient, patient, param_set))


plt.savefig('plots_with_scores_only_kneepoint_patient_specific\%s_%s_model_fit_with_scores.pdf' % (patient, param_set))
plt.savefig('plots_with_scores_only_kneepoint_patient_specific\%s_%s_model_fit_with_scores.png' % (patient, param_set))

stop = timeit.default_timer()
print(stop - start)  


# ### get simulation data
# sim_data = pd.read_csv(simu_file_path, skiprows=6, sep=",", header=0)
# print(sim_data['[run number]'])

# runs = sim_data['[run number]'].unique()
# print(sorted(runs))

# for run in runs :
# 	sim_viability = pd.DataFrame()

remove_duplicates_generic_with_filtering_keeping_only_samples.py 
#!/usr/bin/python
# -*- coding: utf-8 -*-
# jupyter notebook --browser='C:/Program Files (x86)/Google/Chrome/Application/chrome.exe %s'


import timeit
start = timeit.default_timer()
import math
import sys
import pylab as plt
import numpy as np
import pandas as pd



inputFile = sys.argv[1]
sortby = sys.argv[2]
samples_min = float(sys.argv[3])
outputFile = "%s_duplicates_removed_filtered_only_samples_kept_%s.txt" % (sys.argv[1][:-4], samples_min)

df = pd.read_csv(inputFile) 
print(df)
df["evolution$samples"] = pd.to_numeric(df["evolution$samples"], downcast="float")

df = df[df["evolution$samples"] >= samples_min]
print(df)


df = df.drop(columns=['evolution$generation'])
# print(df)

# dropping ALL duplicate values
df = df.drop_duplicates()
# print(df)

# column_mapdict = {
# "A": "a", 
# "B": "c"
# }

# df = df.rename(columns=column_mapdict)



df.sort_values("%s" % sortby, inplace = True)
print(df)

df.to_csv(outputFile,index = False)




stop = timeit.default_timer()
print(stop - start)  

scikitlearn_plot_predictions_model0.py 
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

df_via_exp = pd.read_csv('allPatients_prediction-via-exp.tsv', sep='	')
df_via_exp['viability_mean'] = df_via_exp.mean(axis=1)
df_via_exp['viability_std'] = df_via_exp.std(axis=1)
df_via_exp = df_via_exp.iloc[::-1]

df_conc_exp = pd.read_csv('allPatients_prediction-conc-exp.tsv', sep='	')
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
# fig.suptitle('B-CLL viability and concentration
 with varying monocytes initial proportions (3 patients)
 Experimental vs. Predictions', fontsize=10, y=1.02)
# fig.suptitle('Model Predictions 
 %s' % file, fontsize=10, y=1.02)
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
axs[0].text(0.1, 0.775, 'NRMSE_via = %s
 R²_via = %s' % (nrms_via,squared_r_via), horizontalalignment='left', verticalalignment='baseline', transform=axs[0].transAxes)
axs[1].text(0.4, 0.775, 'NRMSE_conc = %s
 R²_conc = %s' % (nrms_conc,squared_r_conc), horizontalalignment='left', verticalalignment='baseline', transform=axs[1].transAxes)

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

shell_plots.py 
#!/usr/bin/python
# -*- coding: utf-8 -*-
# jupyter notebook --browser='C:/Program Files (x86)/Google/Chrome/Application/chrome.exe %s'


import timeit
start = timeit.default_timer()

exp_list = ["perturb-gui-apo-mov","perturb-gui-need-sig-mov","perturb-gui-layers","perturb-gui-alpha","perturb-gui-mono-phago-eff","perturb-gui-NLC-phago-eff","perturb-gui-M-phago-eff","perturb-gui-M-kill-eff","perturb-gui-cll-sens-dist","perturb-gui-mono-sens-dist","perturb-gui-nlc-sens-dist","perturb-gui-macro-sens-dist","perturb-gui-nlc-threshold","perturb-gui-sig-init","perturb-gui-sig-init-std","perturb-gui-diff-mean","perturb-gui-diff-std","perturb-gui-life-init-gamma","perturb-gui-alpha-distrib"]


with open("plot_sensitivity_experiments.sh" , 'a+') as file_write :
	file_write.write("#!/bin/bash
")
	for exp in exp_list : 

		file_write.write("echo \"This is a shell script for exp %s\"
" % exp)
		file_write.write(" python sensitivity_analysis/plot_sensitivity_analysis.py %s
" % exp)



stop = timeit.default_timer()
print(stop - start)

shell_plots_class1.py 
#!/usr/bin/python
# -*- coding: utf-8 -*-
# jupyter notebook --browser='C:/Program Files (x86)/Google/Chrome/Application/chrome.exe %s'


import timeit
start = timeit.default_timer()

exp_list = ["perturb-gui-apo-mov","perturb-gui-need-sig-mov","perturb-gui-layers","perturb-gui-alpha","perturb-gui-mono-phago-eff","perturb-gui-NLC-phago-eff","perturb-gui-M-phago-eff","perturb-gui-M-kill-eff","perturb-gui-cll-sens-dist","perturb-gui-mono-sens-dist","perturb-gui-nlc-sens-dist","perturb-gui-macro-sens-dist","perturb-gui-nlc-threshold","perturb-gui-sig-init","perturb-gui-sig-init-std","perturb-gui-diff-mean","perturb-gui-diff-std","perturb-gui-life-init-gamma","perturb-gui-alpha-distrib"]


with open("plot_sensitivity_class_1_experiments.sh" , 'a+') as file_write :
	file_write.write("#!/bin/bash
")
	for exp in exp_list : 

		file_write.write("echo \"This is a shell script for exp %s\"
" % exp)
		file_write.write(" python sensitivity_analysis/plot_sensitivity_analysis_class1.py %s
" % exp)



stop = timeit.default_timer()
print(stop - start)

shell_plots_class2.py 
#!/usr/bin/python
# -*- coding: utf-8 -*-
# jupyter notebook --browser='C:/Program Files (x86)/Google/Chrome/Application/chrome.exe %s'


import timeit
start = timeit.default_timer()

exp_list = ["perturb-gui-apo-mov","perturb-gui-need-sig-mov","perturb-gui-layers","perturb-gui-alpha","perturb-gui-mono-phago-eff","perturb-gui-NLC-phago-eff","perturb-gui-M-phago-eff","perturb-gui-M-kill-eff","perturb-gui-cll-sens-dist","perturb-gui-mono-sens-dist","perturb-gui-nlc-sens-dist","perturb-gui-macro-sens-dist","perturb-gui-nlc-threshold","perturb-gui-sig-init","perturb-gui-sig-init-std","perturb-gui-diff-mean","perturb-gui-diff-std","perturb-gui-life-init-gamma","perturb-gui-alpha-distrib"]


with open("plot_sensitivity_class_2_experiments.sh" , 'a+') as file_write :
	file_write.write("#!/bin/bash
")
	for exp in exp_list : 

		file_write.write("echo \"This is a shell script for exp %s\"
" % exp)
		file_write.write(" python sensitivity_analysis/plot_sensitivity_analysis_class2.py %s
" % exp)



stop = timeit.default_timer()
print(stop - start)

t_test.py 
#!/usr/bin/python
# -*- coding: utf-8 -*-
# jupyter notebook --browser='C:/Program Files (x86)/Google/Chrome/Application/chrome.exe %s'


import timeit
start = timeit.default_timer()

import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import ttest_ind

parameters_sets = pd.read_csv('patient_param_sets_9patients_param_renamed.csv', sep=';', index_col=0) 


parameters_sets = parameters_sets.T
# parameters_sets["class"] = parameters_sets["class"].astype("category")
print(parameters_sets)

group1 = parameters_sets[parameters_sets['class']==1.0]
group2 = parameters_sets[parameters_sets['class']==2.0]
# print(group1)

#perform independent two sample t-test
# for col in parameters_sets.columns :
# 	print("%s;%s" % (col, str(ttest_ind(group1[col], group2[col])).split("=")[2][:-1]))
	# print("%s;%s" % (col, str(ttest_ind(group1[col], group2[col], equal_var=False)).split("=")[2][:-1]))

# plt.show()


from scipy import stats
for col in parameters_sets.columns :
	print("%s;%s" % (col, str(stats.ks_2samp(group1[col], group2[col])).split("=")[2][:-1]))



stop = timeit.default_timer()
print(stop - start)  
