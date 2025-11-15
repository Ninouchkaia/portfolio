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
  print("patient: %s\n" % patient)
  for param_set in param_sets :
    print("param_set: %s\n" % param_set)
    simu_file_path = "%s\\BehaviorSpace\\%s.csv" % (patient, param_set) #(patient,param_set)
    # print('simu_file_path', simu_file_path)

    ## get patient experimental data
    patient_data = pd.read_csv('%s\\%s.csv' % (patient,patient), index_col=0) #, index_col=0) 
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
    print("mean_via_simu\n", mean_via_simu)
    mean_via_exp = patient_data['%s_viability' % patient].transpose()
    print("mean_via_exp\n", mean_via_exp)

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
    print("mean_conc_simu\n", mean_conc_simu)
    mean_conc_exp = patient_data['%s_concentration' % patient].transpose()
    print("mean_conc_exp\n", mean_conc_exp)
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

df_NRMSE_via_max_min.to_csv('NRMSE_via_max_min.tsv', sep='\t', index=True)  
df_NRMSE_via_mean.to_csv('NRMSE_via_mean.tsv', sep='\t', index=True)  
df_NRMSE_via_stdev.to_csv('NRMSE_via_stdev.tsv', sep='\t', index=True)  

df_NRMSE_conc_max_min.to_csv('NRMSE_conc_max_min.tsv', sep='\t', index=True)  
df_NRMSE_conc_mean.to_csv('NRMSE_conc_mean.tsv', sep='\t', index=True)  
df_NRMSE_conc_stdev.to_csv('NRMSE_conc_stdev.tsv', sep='\t', index=True)  

df_NRMSE_sum_max_min.to_csv('NRMSE_sum_max_min.tsv', sep='\t', index=True)  
df_NRMSE_sum_mean.to_csv('NRMSE_sum_mean.tsv', sep='\t', index=True)  
df_NRMSE_sum_stdev.to_csv('NRMSE_sum_stdev.tsv', sep='\t', index=True)  



stop = timeit.default_timer()
print(stop - start)  