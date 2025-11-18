#!/usr/bin/python
# -*- coding: utf-8 -*-
# jupyter notebook --browser='C:/Program Files (x86)/Google/Chrome/Application/chrome.exe %s'

import timeit
start = timeit.default_timer()

from scipy.cluster.hierarchy import linkage, dendrogram, fcluster
from scipy.spatial.distance import squareform
import glob
import statistics

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
 

with open("result6_fillna_control_renamed_filtered4_exp130921_T0_renamed.csv") as file_read:
    data = file_read.readlines()
    ### parse headers
    drugs_list, exp_list = [], []
    drugs_index_dict, exp_index_dict = {}, {}
    barcode_dict = {}
    index_content_dict = {}
    drug_dict = {}
    exp_dict = {}
    exp_dict_full = {}
    for line in data[:1] :
        line = line.replace("\n","").split(";")
        for condition in line[1:] :
            conditions = condition.split("_")
            drug = conditions[0][:-1]
            replicate = conditions[0][-1]
            # print(drug,replicate)
            drugs_list.append(drug)
            exp = conditions[2]
            exp_list.append(exp)
            dose = conditions[1]
            run = conditions[3]
            sample = conditions[4]
            drug_dose = "%s_%s" % (drug,dose)
            
            if drug_dose not in drugs_index_dict :
                drugs_index_dict[drug_dose] = [line.index(condition)]
            else :
                drugs_index_dict[drug_dose].append(line.index(condition))
            if exp not in exp_index_dict :
                exp_index_dict[exp] = [line.index(condition)]
            else :
                exp_index_dict[exp].append(line.index(condition))

            index_content_dict[line.index(condition)] = '_'.join([drug,replicate,dose,exp,run,sample])

            ################# fill the drug_dict ############################
            if drug_dose not in drug_dict :
                drug_dict[drug_dose] = {}
                
            if exp not in drug_dict[drug_dose] :
                drug_dict[drug_dose][exp] = {}
            
            if replicate not in drug_dict[drug_dose][exp] :
                # drug_dict[drug][exp][replicate] = [(line.index(condition),'_'.join([drug,replicate,dose,exp,run,sample]))]
                # drug_dict[drug][exp][replicate] = [(line.index(condition))]
                
                drug_dict[drug_dose][exp][replicate] = (line.index(condition))
            
            # else :
            #   drug_dict[drug][exp][replicate].append(line.index(condition))
            #################################################################

            ################# fill the exp_dict #############################
            if exp not in exp_dict :
                exp_dict[exp] = {}
                exp_dict_full[exp] = {}
                
            if drug_dose not in exp_dict[exp] :
                exp_dict[exp][drug_dose] = {}
                exp_dict_full[exp][drug_dose] = {}
            
            if replicate not in exp_dict[exp][drug_dose] :
                exp_dict_full[exp][drug_dose][replicate] = [(line.index(condition),'_'.join([drug_dose,replicate,exp,run,sample]))]
                # exp_dict[exp][drug][replicate] = [line.index(condition)]

                exp_dict[exp][drug_dose][replicate] = line.index(condition)

            # else :
            #   exp_dict[exp][drug][replicate].append(line.index(condition))

                


    drugs_list = list(set(drugs_list))
    exp_list = list(set(exp_list))

print("drugs_list", drugs_list)
print("exp_list", exp_list) 

# for drug in drugs_index_dict :
#   print(drug, drugs_index_dict[drug])

# for exp in exp_index_dict :
#   print(exp, exp_index_dict[exp])

print("########  index_content_dict ########")
for index in index_content_dict :
    print(index, index_content_dict[index])
    # print(index, index_content_dict[index][0], index_content_dict[index][1], index_content_dict[index][2])

print("########  drug_dict ########")
for drug_dose in drug_dict :
    print(drug_dose, drug_dict[drug_dose])

print("########  exp_dict ########")
for exp in exp_dict :
    print(exp, exp_dict[exp])

#### build file displaying the reads numbers as fold changes compared to 

with open("result6_fillna_control_renamed_filtered4_exp130921_T0_renamed.csv") as file_read:
    data = file_read.readlines()
    barcode_dict, barcode_dict_avg, barcode_dict_avg_fc = {}, {}, {}
    # barcode_dict_avg_short_names = {}
    counter = 0
    # for exp in exp_dict_full :
    #   print(exp)
    #   col_names = []
    #   indexes_list = []
    #   drug_doses = list(exp_dict_full[exp].keys())
    #   for drug_dose in drug_doses :
    #       replicates = list(exp_dict_full[exp][drug_dose].keys())
    #       for replicate in replicates :
    #           col_name = exp_dict_full[exp][drug_dose][replicate][0][1]
    #           col_names.append(col_name)
    #           indexes_list.append(exp_dict_full[exp][drug_dose][replicate][0][0])

    # print(indexes_list)
    # print(col_names)

    for line in data[1:] :
        line = line.replace("\n","").split(";")
        barcode = line[0]
        reads = list(map(float, line[1:])) ## les index ici sont décalés de 1 par rapport aux valeurs des keys de mon index_content_dict
        barcode_dict[barcode] = {}
        barcode_dict_avg[barcode] = {}
        # barcode_dict_avg_short_names[barcode] = {}
        barcode_dict_avg_fc[barcode] = {}

        for exp in exp_dict :
            barcode_dict[barcode][exp] = {}
            contro_indexes = exp_dict[exp]['Contro_000u'].values()
            avg_controls = statistics.mean([reads[i-1] for i in contro_indexes])

            for drug in exp_dict[exp] :
                condition = "%s_%s" % (drug, exp)
                # condition_short = drug
                barcode_dict[barcode][exp] = {}
                replicate_reads = []
                for replicate in exp_dict[exp][drug] :
                    replicate_index = int(exp_dict[exp][drug][replicate])
                    replicate_reads.append(reads[replicate_index - 1])
                avg_replicates = statistics.mean(replicate_reads)
                if avg_controls == 0 :
                    # fold_change = float('NaN')
                    fold_change = float(0)
                    counter = counter + 1
                else :
                    fold_change = avg_replicates / avg_controls
                    

                barcode_dict[barcode][exp][drug] = {'avg_controls' : avg_controls, 'avg_replicates' : avg_replicates, 'fold_change' : fold_change}
                
                if condition not in barcode_dict_avg[barcode] :
                    barcode_dict_avg[barcode][condition] = avg_replicates
                else :
                    print("PROBLEM")

                # if condition not in barcode_dict_avg_short_names[barcode] :
                #     barcode_dict_avg_short_names[barcode][condition] = avg_replicates
                # else :
                #     print("PROBLEM")

                if condition not in barcode_dict_avg_fc[barcode] :
                    barcode_dict_avg_fc[barcode][condition] = fold_change
                else :
                    print("PROBLEM")


# with open("result6_fillna_control_renamed_filtered4_exp130921_T0_renamed_replicates_averages.csv", 'w') as file_write :
 
#     file_write.write("%s;%s\n" % (barcode, ';'.join("{:.4f}".format(i) for i in reads_to_fold_changes)))
# print("counter", counter)
# final_df = pd.DataFrame(barcode_dict_avg)
# final_df = final_df.T
# print (final_df)

# print("counter", counter)
# final_df_short_names = pd.DataFrame(barcode_dict_avg_short_names)
# final_df_short_names = final_df_short_names.T
# print (final_df_short_names)

# final_df_short_names.to_csv('result6_fillna_control_renamed_filtered4_exp130921_T0_renamed_replicates_averages_short_names.csv', sep=';', index=True)

final_df_fc = pd.DataFrame(barcode_dict_avg_fc)
final_df_fc = final_df_fc.T
print (final_df_fc)

final_df_fc.to_csv('result6_fillna_control_renamed_filtered4_exp130921_T0_renamed_replicates_averages_fc_zeros.csv', sep=';', index=True)


# print(barcode_dict['CAAGTAGACGATTAGCATTGACTGAAACATGGCAGACGCGA'])


stop = timeit.default_timer()
print(stop - start) 