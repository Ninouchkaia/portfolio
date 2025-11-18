
all_scripts_with_names_barcodes.py 


analyze_barcode_var_avged_over_xps.py 
#!/usr/bin/python
# -*- coding: utf-8 -*-
# jupyter notebook --browser='C:/Program Files (x86)/Google/Chrome/Application/chrome.exe %s'

import timeit
start = timeit.default_timer()

import scipy.stats as stats 
import numpy as np
from statistics import mean
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


with open("barcode_variability_within_controls_aveaged_across_xps.csv") as file_read :
	data = file_read.readlines()
	avg_var_dict = {}
	for line in data[1:] :
		line = line.replace("
","").split("	")
		barcode = line[0]
		avged_var = float(line[1])
		if barcode in avg_var_dict :
			print("KO !!")
		avg_var_dict[barcode] = avged_var

list_avg_var = []
for barcode in avg_var_dict :
	list_avg_var.append(avg_var_dict[barcode])
plt.hist(list_avg_var, bins=500)
plt.title('Frequency Plot for barcode variability (max-min)/avg 
 across controls averaged over experiments (500 bins)')
plt.savefig("freq_plot_averaged_variability_over_xps.pdf")
plt.savefig("freq_plot_averaged_variability_over_xps.png")


# density = stats.kde.gaussian_kde(list_avg_var)
# x = np.arange(min(list_avg_var), max(list_avg_var), 0.1)
# plt.scatter(x, density(x))
# plt.title('Density Plot for averaged variability (max-min)/avg 
of barcodes in controls over experiments')
# plt.savefig("averaged_variability_over_xps.pdf")
# plt.show()

stop = timeit.default_timer()
print(stop - start) 

analyze_runs.py 
#!/usr/bin/python
# -*- coding: utf-8 -*-
# jupyter notebook --browser='C:/Program Files (x86)/Google/Chrome/Application/chrome.exe %s'


import timeit
start = timeit.default_timer()

import pandas as pd
import glob
path = "A:\Downloads\Projects\workFromHome\Projects\drug_screening\data\*.csv"
for fname in glob.glob("*.csv"):
    print(fname)
    df = pd.read_csv(fname, sep=';')
    experiments = []
    for col in df.columns :
        if '_exp' in col :
            # print(col)
            # print(col[13:22])
            if col[13:22] not in experiments :
                experiments.append(col[13:22])
    print(experiments)




stop = timeit.default_timer()
print(stop - start)  

avg_replicates.py 
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
        line = line.replace("
","").split(";")
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
        line = line.replace("
","").split(";")
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
 
#     file_write.write("%s;%s
" % (barcode, ';'.join("{:.4f}".format(i) for i in reads_to_fold_changes)))
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

avg_variability_across_xp_density_plots.py 
#!/usr/bin/python
# -*- coding: utf-8 -*-
# jupyter notebook --browser='C:/Program Files (x86)/Google/Chrome/Application/chrome.exe %s'

import timeit
start = timeit.default_timer()

import scipy.stats as stats 
import numpy as np
from statistics import mean
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


xp_list = ["exp300821","exp010821","exp040821","exp151121","exp271221","exp130921","exp200921",
"exp181021"] 

with open("barcodes_in_controls_stdev_4_sorted_barcode.csv") as file_read :
	data = file_read.readlines()
	avg_var_dict = {}
	for line in data[1:] :
		line = line.replace("
","").split(";")
		experiment = line[0]
		barcode = line[1]
		avged_var = float(line[8])
		if barcode not in avg_var_dict :
			avg_var_dict[barcode] = [avged_var]
		else : 
			avg_var_dict[barcode].append(avged_var)

list_averages=[]
for barcode in avg_var_dict:
	list_averages.append(mean(avg_var_dict[barcode]))
	if len(avg_var_dict[barcode]) != 8 :
		print(avg_var_dict[barcode])

# with open("barcode_variability_within_controls_aveaged_across_xps.csv", "a+") as file_write :
# 	file_write.write("barcode	avged_var
")
# 	for barcode in avg_var_dict :
# 		file_write.write("%s	%s
" % (barcode,mean(avg_var_dict[barcode])))

# plt.hist(list_averages, bins=5000)

density = stats.kde.gaussian_kde(list_averages)
x = np.arange(min(list_averages), max(list_averages), 0.01)
plt.scatter(x, density(x))
plt.title('Density Plot for averaged variability (max-min)/avg 
of barcodes in controls over experiments')
plt.savefig("averaged_variability_over_xps.pdf")

stop = timeit.default_timer()
print(stop - start) 

avoid_zero_reads.py 


import timeit
start = timeit.default_timer()

import pandas as pd
import glob

with open("A:\Downloads\Projects\workFromHome\Projects\drug_screening\2022_Barcodes\data\combined_runs_filtered_avoid_zero_reads.csv", 'a+') as file_write:
    with open("A:\Downloads\Projects\workFromHome\Projects\drug_screening\2022_Barcodes\data\combined_runs_filtered.csv") as file_read:
        data = file_read.readlines()
        file_write.write(data[0])
        
        for line in data[:1] : # we parse only the column names, 1st line of the matrix
            line = line.replace("
","").split(";")

            # create a list of all experiment names
            exp_list = []
            for full_id in line[1:] : # we start from 1 to avoid the empty first cell
                #(print(full_id))

                experiment_name = (full_id.split("_"))[2] # we get the experiment name
                exp_list.append(experiment_name) # and append it to exp_list
            exp_list = list(set(exp_list))
            # print(len(exp_list))
            # for exp in exp_list :
            # print(exp)

            # get indexes of the treatment replicates for each drug/dose/exp combination and store in a dict
            drug_dose_exp_dict = {}
            for i in range(1,len(line)): # we start from 1 to avoid the empty first cell
                #print(line[i])
                drug = line[i].split("_")[0][:-1]
                dose = line[i].split("_")[1]
                exp = line[i].split("_")[2]
                drug_dose_exp = "_".join([drug, dose, exp])
                if drug_dose_exp not in drug_dose_exp_dict :
                    drug_dose_exp_dict[drug_dose_exp] = [i]
                else :
                    drug_dose_exp_dict[drug_dose_exp].append(i)
            # print(drug_dose_exp_dict)
            for drug_dose_exp in drug_dose_exp_dict :
                print("%s: %s" % (drug_dose_exp, drug_dose_exp_dict[drug_dose_exp]))

            # get indexes of the controls for each exp and store in a dict
            control_dict = {}
            for k in range(1,len(line)): # we start from 1 to avoid the empty first cell
                #print(line[k])
                exp = line[k].split("_")[2]
                if "Contro" in line[k] :
                    if exp not in control_dict :
                        control_dict[exp] = [k]
                    else :
                        control_dict[exp].append(k)
            # print(control_dict)
            # for exp in control_dict :
            #     print("%s: %s" % (exp, control_dict[exp]))

        # in each experiment, if sum reads for a treatment is zero and sum controls is non-zero, store the barcode and the treatment id in a dict
        # and write a new matrix converting the zero reads into 0.001
        zero_reads_dict = {}
        for line in data[1:] :
            line = line.replace("
","").split(";")
            line_fill = [i for i in line]
            for drug_dose_exp in drug_dose_exp_dict :
                #print(drug_dose_exp)
                
                # get experiment id
                experiment = drug_dose_exp.split("_")[2]
                
                # get control reads indexes from that experiment
                control_reads_indexes = control_dict[experiment]
                
                # make sum of controls
                sum_controls = sum([float(line[j]) for j in control_reads_indexes])
                
                # if sum controls is non-zero, then check sum of treatments
                if (sum_controls != 0) :

                    treatment_indexes = drug_dose_exp_dict[drug_dose_exp]
                    reads = [float(line[j]) for j in treatment_indexes]
                    sum_reads = sum(reads)
                    #print(sum_reads)
                    if (float(sum_reads) == float(0)) :
                        barcode = line[0]
                        #print(barcode, drug_dose_exp)
                        if (barcode not in zero_reads_dict) :
                            zero_reads_dict[barcode] = [drug_dose_exp]
                        else :
                            zero_reads_dict[barcode].append(drug_dose_exp)
                        # replace in line_fill the zero values by 0.001
                        for h in treatment_indexes :
                            line_fill[h] = '0.01'
            
            # fill the new line in new file
            file_write.write(';'.join(line_fill))
            file_write.write('
')


        # for barcode in zero_reads_dict :
        #     print(barcode, zero_reads_dict[barcode])

        #print(len(zero_reads_dict)) # 11777 (out of 12305... seems quite high)

stop = timeit.default_timer()
print(stop - start)  

avoid_zero_reads_by_replacing_all_zeros.py 


import timeit
start = timeit.default_timer()

import pandas as pd
import glob

with open("A:\Downloads\Projects\workFromHome\Projects\drug_screening\2022_Barcodes\data\combined_runs_filtered_avoid_zero_reads_all_replaced.csv", 'w') as file_write:
    with open("A:\Downloads\Projects\workFromHome\Projects\drug_screening\2022_Barcodes\data\combined_runs_filtered.csv") as file_read:
        data = file_read.readlines()
        file_write.write(data[0])
        
        ## write a new matrix converting the zero reads into 0.001
        for line in data[1:] :
            line = line.replace("
","").split(";")
            line_fill = [i for i in line]
            for i in range(1, len(line)) :
                if float(line[i]) == float(0) :
                    line_fill[i] = str(0.01)

            # fill the new line in new file
            file_write.write(';'.join(line_fill))
            file_write.write('
')




stop = timeit.default_timer()
print(stop - start)  

barcodes_in_controls_heatmap.py 
#!/usr/bin/python
# -*- coding: utf-8 -*-
# jupyter notebook --browser='C:/Program Files (x86)/Google/Chrome/Application/chrome.exe %s'

import timeit
start = timeit.default_timer()

from scipy.cluster.hierarchy import linkage, dendrogram, fcluster
from scipy.spatial.distance import squareform
import glob
import statistics


import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
 

experiments = ['exp200921','exp130921','exp300821','exp040821','exp181021','exp151121','exp271221','exp010821']
barcodes_dict = {}
# with open("barcodes_in_controls_stdev_4.csv") as file_read:
with open("barcodes_in_controls_stdev_4_from_fold_changes.csv") as file_read:

	data = file_read.readlines()
	for experiment in experiments :
		for line in data[1:] :
			line = line.replace("
","").split("	")
			# print(line)
			xp_id = line[0]
			if xp_id == experiment :
				barcode = line[1]
				stdev = float(line[6])
				max_min = float(line[8])
				if barcode in barcodes_dict :
					barcodes_dict[barcode].append(max_min)
				else :
					barcodes_dict[barcode] = [max_min]

# print(barcodes_dict)

barcodes_df = pd.DataFrame.from_dict(barcodes_dict, orient='index',columns=experiments)
print(barcodes_df)
barcodes_df = barcodes_df.dropna()
print(barcodes_df)

filtered_barcodes_df = barcodes_df[barcodes_df < 2] 
filtered_barcodes_df = filtered_barcodes_df.dropna()
print(filtered_barcodes_df)

# sns.heatmap(barcodes_df, yticklabels=False)
# sns.heatmap(filtered_barcodes_df, yticklabels=False)

# sns.clustermap(barcodes_df,yticklabels=False)
sns.clustermap(filtered_barcodes_df,yticklabels=False)
plt.show()

# plt.savefig("clustermap_barcodes_in_controls_max_min.pdf", format='pdf')

		
stop = timeit.default_timer()
print(stop - start) 

build_design.py 
#!/usr/bin/python
# -*- coding: utf-8 -*-
# jupyter notebook --browser='C:/Program Files (x86)/Google/Chrome/Application/chrome.exe %s'

import timeit
start = timeit.default_timer()


 
with open("design4.csv", 'w') as file_write :
	file_write.write(";exp;replicate;condition
")
	with open("A:\Downloads\Projects\workFromHome\Projects\drug_screening\2022_Barcodes\data\combined_runs_normalized_filtered.csv") as file_read:
		data = file_read.readlines()
		for line in data[:1] :
			line = line.replace("
","").split(";")
			for condition in line[1:] :
				file_write.write("%s;" % condition)
				conditions = condition.split("_")
				exp = conditions[2]
				drug = conditions[0][:-1]
				replicate = conditions[0][-1]
				dosage = conditions[1]
				file_write.write("%s;" % exp)
				file_write.write("%s;" % replicate)
				if 'Contro' in drug :
					file_write.write("control
")
				else :
					file_write.write("%s_%s
" % (drug,dosage) )


			


stop = timeit.default_timer()
print(stop - start) 

build_design_2023.py 
#!/usr/bin/python
# -*- coding: utf-8 -*-
# jupyter notebook --browser='C:/Program Files (x86)/Google/Chrome/Application/chrome.exe %s'

import timeit
start = timeit.default_timer()


 
with open("design_6_2023.csv", 'w') as file_write :
	file_write.write(";run;exp;replicate;condition
")
	with open("A:\Downloads\Projects\workFromHome\Projects\drug_screening\2022_Barcodes\data\combined_runs_filtered_avoid_zero_reads_x100.csv") as file_read:
		data = file_read.readlines()
		for line in data[:1] :
			line = line.replace("
","").split(";")
			for condition in line[1:] :
				file_write.write("%s;" % condition)
				conditions = condition.split("_")
				run = conditions[3]
				exp = conditions[2]
				drug = conditions[0][:-1]
				replicate = conditions[0][-1]
				dosage = conditions[1]
				file_write.write("%s;" % run)
				file_write.write("%s;" % exp)
				file_write.write("%s;" % replicate)
				if 'Contro' in drug :
					file_write.write("control
")
					# file_write.write("control_%s_%s
" % (exp,run))
				else :
					file_write.write("%s_%s
" % (drug,dosage) )
					# file_write.write("%s_%s_%s_%s
" % (drug,dosage,exp,run) )


			


stop = timeit.default_timer()
print(stop - start) 

check_df.py 
#!/usr/bin/python
# -*- coding: utf-8 -*-
# jupyter notebook --browser='C:/Program Files (x86)/Google/Chrome/Application/chrome.exe %s'


import timeit
start = timeit.default_timer()

import pandas as pd
import glob
import matplotlib.pyplot as plt

for fname in glob.glob("A:\Downloads\Projects\workFromHome\Projects\drug_screening\2022_Barcodes\data\combined_runs_normalized_filtered.csv"):
    print(fname)
    df = pd.read_csv(fname, sep=';', header=0, index_col=0)
    # df2 = df.sum(axis = 0, skipna = True)
    # ax =df2.plot(kind='line')
    # x_axis = ax.axes.get_xaxis()
    # x_axis.set_visible(False)
    # plt.show()
    






stop = timeit.default_timer()
print(stop - start)  

check_df_2023.py 
#!/usr/bin/python
# -*- coding: utf-8 -*-
# jupyter notebook --browser='C:/Program Files (x86)/Google/Chrome/Application/chrome.exe %s'


import timeit
start = timeit.default_timer()

import pandas as pd
import glob
import matplotlib.pyplot as plt

for fname in glob.glob("A:\Downloads\Projects\workFromHome\Projects\drug_screening\2022_Barcodes\data\combined_runs.csv"):
    print(fname)
    df = pd.read_csv(fname, sep=';', header=0, index_col=0)
    for col in df.columns:
    	print(col)


stop = timeit.default_timer()
print(stop - start)  

check_unique_barcodes.py 
#!/usr/bin/python
# -*- coding: utf-8 -*-
# jupyter notebook --browser='C:/Program Files (x86)/Google/Chrome/Application/chrome.exe %s'

import timeit
start = timeit.default_timer()

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


df = pd.read_csv("result6_fillna_control_renamed_filtered4.csv", sep=';', header=0, index_col=0)
print(df)

df_control = df
df_control.drop([col for col in df_control.columns if 'Contro' not in col],axis=1,inplace=True)
print(df_control)

control_unique_barcodes = []
for column_name in df_control.columns:
    column = df_control[column_name]
    # Get the count of Zeros in column 
    count = (column == 0).sum()
    # print(column_name, '0:', count, 'non-0:', len(df_control.index)-count)
    control_unique_barcodes.append(len(df_control.index)-count)


df = pd.read_csv("result6_fillna_control_renamed_filtered4.csv", sep=';', header=0, index_col=0)
df_drug = df
df_drug.drop([col for col in df_drug.columns if 'Contro' in col],axis=1,inplace=True)
df_drug.drop([col for col in df_drug.columns if 'Temps0' in col],axis=1,inplace=True)

print(df_drug)

drug_unique_barcodes = []
for column_name in df_drug.columns:
    column = df_drug[column_name]
    # Get the count of Zeros in column 
    count = (column == 0).sum()
    # print(column_name, '0:', count, 'non-0:', len(df_drug.index)-count)
    drug_unique_barcodes.append(len(df_drug.index)-count)


print("control_unique_barcodes", control_unique_barcodes)
print("drug_unique_barcodes", drug_unique_barcodes)

control_unique_barcodes_arr = np.array(control_unique_barcodes)
drug_unique_barcodes_arr = np.array(drug_unique_barcodes)

data = [drug_unique_barcodes_arr,control_unique_barcodes_arr]

# fig = plt.figure(figsize =(10, 7))
# ax = fig.add_axes([0,0,1, 1])

red_square = dict(markerfacecolor='r', marker='s')
fig, ax = plt.subplots()
ax.set_title('Unique Barcodes in Controls vs. Drugs')
ax.boxplot(data, vert=False, flierprops=red_square)

# fig.canvas.draw()

labels = ["Drugs","Controls"]

ax.set_yticklabels(labels)


 
# show plot
plt.savefig("unique_barcodes.png")



stop = timeit.default_timer()
print(stop - start) 

clean_combined_file.py 
#!/usr/bin/python
# -*- coding: utf-8 -*-
# jupyter notebook --browser='C:/Program Files (x86)/Google/Chrome/Application/chrome.exe %s'


import timeit
start = timeit.default_timer()

import pandas as pd
import glob

print('HELLO')
for fname in glob.glob("/home/nina/Projects/2022_Barcodes/NetBioMed2022/data/combined_runs.csv"):
    print(fname)
    df = pd.read_csv(fname, sep=';', header=0, index_col=0)
    print(df)
    df = df.fillna(0).apply(lambda x: 1000000 * x / float(x.sum()))

    df.to_csv('/home/nina/Projects/2022_Barcodes/NetBioMed2022/data/combined_runs_normalized.csv', sep=';', index=True, float_format='%.2f')

stop = timeit.default_timer()
print(stop - start)   

compare_dataframes.py 
#!/usr/bin/python
# -*- coding: utf-8 -*-
# jupyter notebook --browser='C:/Program Files (x86)/Google/Chrome/Application/chrome.exe %s'


import timeit
start = timeit.default_timer()

import pandas as pd
import glob
import matplotlib inline
from matplotlib_venn import venn2, venn2_circles
from matplotlib_venn import venn3, venn3_circles
from matplotlib import pyplot as plt


rownames_list = []
path = "A:\Downloads\Projects\workFromHome\Projects\drug_screening\data\*.csv"
for fname in glob.glob("*.csv"):
# for fname in glob.glob("result6.csv"):
    print(fname)
    df = pd.read_csv(fname, sep=';', header=0, index_col=0)
    # get the row names
    row_names = df.index.values.tolist()
    rownames_list.append(row_names)

for row_names in rownames_list :

stop = timeit.default_timer()
print(stop - start)  

convert_fold_changes_to_classes.py 
#!/usr/bin/python
# -*- coding: utf-8 -*-
# jupyter notebook --browser='C:/Program Files (x86)/Google/Chrome/Application/chrome.exe %s'

import timeit
start = timeit.default_timer()

import scipy.stats as stats 
import numpy as np
from statistics import mean
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

with open("fold_changes_classes.csv", 'a+') as file_write :
	with open("fold_changes_agreggated.csv") as file_read:
		data = file_read.readlines()
		file_write.write(data[0])
		for line in data[1:] :
			line = line.replace("
","").split(';')
			barcode = line[0]
			fc_class_list = []
			for fc in line[1:] :
				# print(fc)
				if fc == '' :
					fc = 'NaN'
				fc = float(fc)
				if (fc > 500) :
					class_fc = 'A'
				elif (500 >= fc > 200) :
					class_fc = 'B'
				elif (20 < fc <= 50) :
					class_fc = 'D'
				elif (fc <= 20) :
					class_fc = 'E'
				elif (50 < fc <= 200) :
					class_fc = 'F'
				fc_class_list.append(class_fc)
			file_write.write("%s;%s
" % (barcode, ';'.join(fc_class_list)))



stop = timeit.default_timer()
print(stop - start) 

corr_dendrogram.py 
#!/usr/bin/python
# -*- coding: utf-8 -*-
# jupyter notebook --browser='C:/Program Files (x86)/Google/Chrome/Application/chrome.exe %s'

import timeit
start = timeit.default_timer()

from scipy.cluster.hierarchy import linkage, dendrogram, fcluster
from scipy.spatial.distance import squareform


import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
 

# df = pd.read_csv("result6_fillna_control_renamed_filtered6.csv", sep=';', header=0, index_col=0)

df1 = pd.read_csv("result6_fillna_control_renamed_filtered4_exp130921_T0_renamed_replicates_averages.csv", sep=';', header=0, index_col=0)
print (df1) #
df1 = df1.astype(float)
df1 = df1.sort_index(axis=1, ascending=True, key=lambda x: x.str.lower())

matrix1 = df1.corr(
    method = 'pearson',  # The method of correlation
    min_periods = 1).round(2)     # Min number of observations required)

# m1 = sns.heatmap(matrix1,annot=False, yticklabels=True, xticklabels=True)
# m1.set_yticklabels(m1.get_ymajorticklabels(), size = 2)
# m1.set_xticklabels(m1.get_xmajorticklabels(), size = 2)
# m1.tick_params(width=0.5 )#left=False, bottom=False)
# m1.set_title('avg_reads_correlations')

# plt.savefig("avg_reads_correlations.pdf", format='pdf')

df2 = pd.read_csv("result6_fillna_control_renamed_filtered4_exp130921_T0_renamed_replicates_averages_fc_zeros.csv", sep=';', header=0, index_col=0)
print (df2) #
df2 = df2.astype(float)
df2 = df2.sort_index(axis=1, ascending=True, key=lambda x: x.str.lower())

matrix2 = df2.corr(
    method = 'pearson',  # The method of correlation
    min_periods = 1).round(2)     # Min number of observations required)

# m2 = sns.heatmap(matrix2,annot=False, yticklabels=True, xticklabels=True)
# m2.set_yticklabels(m2.get_ymajorticklabels(), size = 2)
# m2.set_xticklabels(m2.get_xmajorticklabels(), size = 2)
# m2.tick_params(width=0.5 )#left=False, bottom=False)
# m2.set_title('fold_changes_correlations')
# plt.savefig("fold_changes_correlations.pdf", format='pdf')


diff = matrix1 - matrix2
h = sns.heatmap(diff,annot=False, yticklabels=True, xticklabels=True)
# h.set_yticklabels(h.get_ymajorticklabels(), size = 2)
# h.set_xticklabels(h.get_xmajorticklabels(), size = 2)
# h.tick_params(width=0.5 )#left=False, bottom=False)
# h.set_title('difference between correlations based on avg_reads vs. fold changes')
# plt.savefig("diff_correl.pdf", format='pdf')


clustered_diff = sns.clustermap(diff, method="complete", cmap='RdBu', annot=False, yticklabels=True, xticklabels=True,
               annot_kws={"size": 7}, vmin=-1, vmax=1, figsize=(15,12));

clustered_diff.ax_heatmap.set_xticklabels(clustered_diff.ax_heatmap.get_xmajorticklabels(), fontsize = 6)
clustered_diff.ax_heatmap.set_yticklabels(clustered_diff.ax_heatmap.get_ymajorticklabels(), fontsize = 6)
clustered_diff.fig.suptitle('"clustered_diff between correlations based on avg_reads vs. fold changes"') 


plt.savefig("clustered_diff.pdf", format='pdf')

plt.show()





# plt.show()

# matrix.to_csv('correl_result6_fillna_control_renamed_filtered4_exp130921_T0_renamed_replicates_averages.csv', sep=';', index=True)


##################  DENDROGRAM ########################################
# plt.figure(figsize=(12,5))
# dissimilarity = 1 - abs(matrix)
# Z = linkage(squareform(dissimilarity), 'complete')

# with plt.rc_context({'lines.linewidth': 0.5}):
# 	dendrogram(Z, labels=df.columns, orientation='top', leaf_rotation=90, leaf_font_size = 1)

# plt.savefig("dendro.pdf", format='pdf')
# # plt.show()
########################################################################

############# CORRPLOT + DENDROGRAM ####################################
# g = sns.clustermap(matrix, method="complete", cmap='RdBu', annot=False, yticklabels=True, xticklabels=True,
#                annot_kws={"size": 7}, vmin=-1, vmax=1, figsize=(15,12));

# g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xmajorticklabels(), fontsize = 6)
# g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_ymajorticklabels(), fontsize = 6)

# plt.savefig("corr_dendro_result6_fillna_control_renamed_filtered4_exp130921_T0_renamed_replicates_averages_fc_zeros.pdf", format='pdf')

# plt.show()
########################################################################

stop = timeit.default_timer()
print(stop - start) 

corr_graph.py 
import timeit
start = timeit.default_timer()

import pandas as pd
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import seaborn as sns
 
df = pd.read_csv("result6_fillna_control_renamed_filtered6.csv", sep=';', header=0, index_col=0)
print (df) #

matrix = df.corr(
    method = 'pearson',  # The method of correlation
    min_periods = 1).round(2)     # Min number of observations required)

# # Transform it in a links data frame (3 columns only):
links = matrix.stack().reset_index()
print(links)
links.columns = ['var1', 'var2', 'value']

# Keep only correlation over a threshold and remove self correlation (cor(A,A)=1)
links_filtered=links.loc[ (links['value'] > 0.6) & (links['var1'] != links['var2']) ]
 
# Build your graph
G=nx.from_pandas_edgelist(links_filtered, 'var1', 'var2')
 
# Plot the network:
nx.draw(G, with_labels=True, node_color='orange', node_size=50, edge_color='black', linewidths=1, font_size=5)
# nx.draw(G, with_labels=False, node_color='orange', node_size=50, edge_color='black', linewidths=1, font_size=5)
plt.show() 

stop = timeit.default_timer()
print(stop - start)  

corr_mat.py 
#!/usr/bin/python
# -*- coding: utf-8 -*-
# jupyter notebook --browser='C:/Program Files (x86)/Google/Chrome/Application/chrome.exe %s'


import timeit
start = timeit.default_timer()

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
 

df = pd.read_csv("result6_fillna_control_renamed_filtered6.csv", sep=';', header=0, index_col=0)
print (df) #

matrix = df.corr(
    method = 'pearson',  # The method of correlation
    min_periods = 1).round(2)     # Min number of observations required)

print(matrix)

sns.heatmap(matrix, annot=True)
plt.show()


stop = timeit.default_timer()
print(stop - start)  

correl_deseq2.py 
#!/usr/bin/python
# -*- coding: utf-8 -*-
# jupyter notebook --browser='C:/Program Files (x86)/Google/Chrome/Application/chrome.exe %s'

import timeit
start = timeit.default_timer()

from scipy.cluster.hierarchy import linkage, dendrogram, fcluster
from scipy.spatial.distance import squareform


import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
 

# df = pd.read_csv("result6_fillna_control_renamed_filtered6.csv", sep=';', header=0, index_col=0)

df1 = pd.read_csv("A:\Downloads\Projects\workFromHome\Projects\drug_screening\2022_Barcodes\results\merged_pval_filtered_deseq2.csv", sep=';', header=0, index_col=0)
print (df1) #
df1 = df1.astype(float)
# df1 = df1.dropna(inplace=True)
df1 = df1.sort_index(axis=1, ascending=True, key=lambda x: x.str.lower())

df1 = df1.fillna(0)
print(df1)
df1.to_csv('merged_pval_filtered_deseq2_fillna.csv', sep=';', index=True)


# matrix1 = df1.corr(
#     method = 'pearson',  # The method of correlation
#     min_periods = 1).round(2)     # Min number of observations required)

# print(matrix1)

# matrix1 = matrix1.fillna(0)
# print(matrix1)

# matrix1 = matrix1.sort_index(axis=1, ascending=True, key=lambda x: x.str.lower())


# matrix1.to_csv('correl_matrix1_merged_pval_filtered_deseq2_fillna.csv', sep=';', index=True)

# m1 = sns.heatmap(matrix1,annot=False, yticklabels=True, xticklabels=True)
# m1.set_yticklabels(m1.get_ymajorticklabels(), size = 2)
# m1.set_xticklabels(m1.get_xmajorticklabels(), size = 2)
# m1.tick_params(width=0.5 )#left=False, bottom=False)
# # m1.set_title('deseq2_logfc_correlations')

# plt.savefig("deseq2_logfc_correlations_pval_filtered_ordered.pdf", format='pdf')


# ##################  DENDROGRAM ########################################
# plt.figure(figsize=(12,5))
# dissimilarity = 1 - abs(matrix1)
# Z = linkage(squareform(dissimilarity), 'complete')

# with plt.rc_context({'lines.linewidth': 0.5}):
# 	dendrogram(Z, labels=df.columns, orientation='top', leaf_rotation=90, leaf_font_size = 1)

# plt.savefig("dendro.pdf", format='pdf')
# plt.show()
# ########################################################################

########### CORRPLOT + DENDROGRAM ####################################
# g = sns.clustermap(matrix1, method="complete", cmap='RdBu', annot=False, yticklabels=True, xticklabels=True,
#                annot_kws={"size": 7}, vmin=-1, vmax=1, figsize=(15,12));

# g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xmajorticklabels(), fontsize = 6)
# g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_ymajorticklabels(), fontsize = 6)
# g.ax_row_dendrogram.set_visible(False)



# plt.savefig("deseq2_pvalfiltered_logfc_correlations_clustered.pdf", format='pdf')

# plt.show()
#######################################################################
# import scipy.spatial as sp, scipy.cluster.hierarchy as hc

# row_dism = 1 - df1.T.corr()
# row_linkage = hc.linkage(sp.distance.squareform(row_dism), method='complete')
# col_dism = 1 - df1.corr()
# col_linkage = hc.linkage(sp.distance.squareform(col_dism), method='complete')

# sns.clustermap(df1,figsize=(5, 5),row_linkage=row_linkage, col_linkage=col_linkage)

stop = timeit.default_timer()
print(stop - start) 

correl_deseq2_2023.py 
#!/usr/bin/python
# -*- coding: utf-8 -*-
# jupyter notebook --browser='C:/Program Files (x86)/Google/Chrome/Application/chrome.exe %s'

import timeit
start = timeit.default_timer()

from scipy.cluster.hierarchy import linkage, dendrogram, fcluster
from scipy.spatial.distance import squareform


import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
 

# df = pd.read_csv("result6_fillna_control_renamed_filtered6.csv", sep=';', header=0, index_col=0)

df1 = pd.read_csv("merged_logfc_pval_filtered_deseq2_2023.csv", sep=';', header=0, index_col=0)
print (df1) #
df1 = df1.astype(float)
# df1 = df1.dropna(inplace=True)
df1 = df1.sort_index(axis=1, ascending=True, key=lambda x: x.str.lower())

df1 = df1.fillna('NA')
df1.drop(columns=df1.columns[0], axis=1,  inplace=True)
print(df1)
df1.to_csv('merged_logfc_pval_filtered_deseq2_2023_fillna.csv', sep=';', index=True)


matrix1 = df1.corr(
    method = 'pearson',  # The method of correlation
    min_periods = 1).round(2)     # Min number of observations required)

print(matrix1)

# matrix1 = matrix1.fillna(0)
# print(matrix1)

# matrix1 = matrix1.sort_index(axis=1, ascending=True, key=lambda x: x.str.lower())


# matrix1.to_csv('correl_matrix1_merged_pval_filtered_deseq2_fillna.csv', sep=';', index=True)

# m1 = sns.heatmap(matrix1,annot=False, yticklabels=True, xticklabels=True)
# m1.set_yticklabels(m1.get_ymajorticklabels(), size = 2)
# m1.set_xticklabels(m1.get_xmajorticklabels(), size = 2)
# m1.tick_params(width=0.5 )#left=False, bottom=False)
# # m1.set_title('deseq2_logfc_correlations')

# plt.savefig("deseq2_logfc_correlations_pval_filtered_ordered_2023.pdf", format='pdf')


# ##################  DENDROGRAM ########################################
# plt.figure(figsize=(12,5))
# dissimilarity = 1 - abs(matrix1)
# Z = linkage(squareform(dissimilarity), 'complete')

# with plt.rc_context({'lines.linewidth': 0.5}):
# 	dendrogram(Z, labels=df.columns, orientation='top', leaf_rotation=90, leaf_font_size = 1)

# plt.savefig("dendro.pdf", format='pdf')
# plt.show()
# ########################################################################

########### CORRPLOT + DENDROGRAM ####################################
# g = sns.clustermap(matrix1, method="complete", cmap='RdBu', annot=False, yticklabels=True, xticklabels=True,
#                annot_kws={"size": 7}, vmin=-1, vmax=1, figsize=(15,12));

# g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xmajorticklabels(), fontsize = 6)
# g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_ymajorticklabels(), fontsize = 6)
# g.ax_row_dendrogram.set_visible(False)



# plt.savefig("deseq2_pvalfiltered_logfc_correlations_clustered.pdf", format='pdf')

# plt.show()
#######################################################################
# import scipy.spatial as sp, scipy.cluster.hierarchy as hc

# row_dism = 1 - df1.T.corr()
# row_linkage = hc.linkage(sp.distance.squareform(row_dism), method='complete')
# col_dism = 1 - df1.corr()
# col_linkage = hc.linkage(sp.distance.squareform(col_dism), method='complete')

# sns.clustermap(df1,figsize=(5, 5),row_linkage=row_linkage, col_linkage=col_linkage)

stop = timeit.default_timer()
print(stop - start) 

count_unique.py 
#!/usr/bin/python
# -*- coding: utf-8 -*-
# jupyter notebook --browser='C:/Program Files (x86)/Google/Chrome/Application/chrome.exe %s'


import timeit
start = timeit.default_timer()

with open("colnames_annotated_2023.csv", 'r') as file_read :
    data = file_read.readlines()
    print(data)
    annot_list = []
    for line in data[1:]:
        line = line.split(";")
        annot = line[1].replace("
","")
        annot_list.append(annot)
    print(annot_list)

    print(len(annot_list))

annot_list = list(set(annot_list))
print(len(annot_list))


stop = timeit.default_timer()
print(stop - start)  

data_analysis.py 
#!/usr/bin/python
# -*- coding: utf-8 -*-
# jupyter notebook --browser='C:/Program Files (x86)/Google/Chrome/Application/chrome.exe %s'


import timeit
start = timeit.default_timer()

import pandas as pd
import glob
path = "A:\Downloads\Projects\workFromHome\Projects\drug_screening\data\*.csv"
# for fname in glob.glob("*.csv"):
for fname in glob.glob("result6_fillna_control_renamed_filtered6.csv"):
    print(fname)
    df = pd.read_csv(fname, sep=';', header=0, index_col=0)
    print (df) #
    # get the control column names
    control_names, time_zeros = [], []
    for col in df.columns :
        if 'CtrlMs' in col :
            print (col)
            control_names.append(col)
        if 'Temps' in col :
            print(col)
            time_zeros.append(col)
    # make a new column the sum of the 4 controls
print(control_names)
print(time_zeros)
for i in control_names :
    print (i)

# df['sum_controls'] = df[control_names[0]] + df[control_names[1]] + df[control_names[2]] + df[control_names[3]]
# df['sum_controls'] = df[control_names[4]] + df[control_names[5]] + df[control_names[6]] + df[control_names[7]]
# df['sum_controls'] = df[control_names[0]] + df[control_names[1]] + df[control_names[2]] + df[control_names[3]] + df[control_names[4]] + df[control_names[5]] + df[control_names[6]] + df[control_names[7]]
# df['sum_controls'] = df[control_names[8]] + df[control_names[9]] + df[control_names[10]] + df[control_names[11]]

# print(df['sum_controls'])

# check if the barcodes are unique
# barcodes = df.iloc[:,0].tolist()
# print(len(barcodes))
# barcodes=list(set(barcodes))
# print(len(barcodes))
# print(barcodes[0:100])

# drop all barcodes (rows) that are not present in the controls
# df0 = df.drop(df[df.sum_controls == 0].index)
# print("drop barcodes that have 0 reads in the summed controls: ", len(df0.index))

# # drop all samples which do not belong to the experiment exp300821
# df0_filtered_exp300821 = df0.drop((x for x in df.columns.tolist() if 'exp300821' not in x), axis=1)
# df0_filtered_exp300821.to_csv('Run210929_filtered_exp300821.csv', sep=';', index=True)




# df1 = df.drop(df[df.sum_controls < 5].index)
# print("drop barcodes that have <5 reads in the summed controls: ", len(df1.index))

# df2 = df.drop(df[df.sum_controls < 10].index)
# print("drop barcodes that have <10 reads in the summed controls: ", len(df2.index))

# df3 = df.drop(df[df.sum_controls < 50].index)
# print("drop barcodes that have <50 reads in the summed controls: ", len(df3.index))

# df4 = df.drop(df[df.sum_controls < 100].index)
# print("drop barcodes that have <100 reads in the summed controls: ", len(df4.index))

####################################################################################
# ### store the retained/filtered seq.barcodes from filtered df0 into a set
# barcodes_filtered = set(df0.iloc[:,0].tolist())

# ## store barcodes from run_211011 into another set
# for fname in glob.glob("Run211011_raw data.csv"):
#     print(fname)
#     df_211011 = pd.read_csv(fname, sep=';')
#     print (df_211011) # [311688 rows x 49 columns]
# barcodes_211011 = df_211011.iloc[:,0].tolist()
# print(len(barcodes_211011))
# barcodes_211011=list(set(barcodes_211011))
# print(len(barcodes_211011))
# barcodes_211011 = set(barcodes_211011)

# intersection_filtered_210929_and_211011 = barcodes_210929_filtered.intersection(barcodes_211011)
# print(len(intersection_filtered_210929_and_211011))
##########################################################################################







stop = timeit.default_timer()
print(stop - start)  

data_analysis_controls1.py 
#!/usr/bin/python
# -*- coding: utf-8 -*-
# jupyter notebook --browser='C:/Program Files (x86)/Google/Chrome/Application/chrome.exe %s'


import timeit
start = timeit.default_timer()

import pandas as pd
import glob
path = "A:\Downloads\Projects\workFromHome\Projects\drug_screening\data\*.csv"
# for fname in glob.glob("*.csv"):
for fname in glob.glob("Run210929_raw data.csv"):
    print(fname)
    df = pd.read_csv(fname, sep=';', header=0, index_col=0)
    print (df) #
    # get the control column names
    control_names = []
    for col in df.columns :
        if 'Contro' in col :
            print (col)
            control_names.append(col)
    # make a new column the sum of the 4 controls
# print(control_names)

df['sum_controls'] = df[control_names[0]] + df[control_names[1]] + df[control_names[2]] + df[control_names[3]]
# df['sum_controls'] = df[control_names[4]] + df[control_names[5]] + df[control_names[6]] + df[control_names[7]]
# df['sum_controls'] = df[control_names[0]] + df[control_names[1]] + df[control_names[2]] + df[control_names[3]] + df[control_names[4]] + df[control_names[5]] + df[control_names[6]] + df[control_names[7]]
# df['sum_controls'] = df[control_names[8]] + df[control_names[9]] + df[control_names[10]] + df[control_names[11]]

# print(df['sum_controls'])

# check if the barcodes are unique
# barcodes = df.iloc[:,0].tolist()
# print(len(barcodes))
# barcodes=list(set(barcodes))
# print(len(barcodes))
# print(barcodes[0:100])

# drop all barcodes (rows) that are not present in the controls
df0 = df.drop(df[df.sum_controls == 0].index)
print("drop barcodes that have 0 reads in the summed controls: ", len(df0.index))

# drop all samples which do not belong to the experiment exp300821
df0_filtered_exp300821 = df0.drop((x for x in df.columns.tolist() if 'exp300821' not in x), axis=1)
df0_filtered_exp300821.to_csv('Run210929_filtered_exp300821.csv', sep=';', index=True)




# df1 = df.drop(df[df.sum_controls < 5].index)
# print("drop barcodes that have <5 reads in the summed controls: ", len(df1.index))

# df2 = df.drop(df[df.sum_controls < 10].index)
# print("drop barcodes that have <10 reads in the summed controls: ", len(df2.index))

# df3 = df.drop(df[df.sum_controls < 50].index)
# print("drop barcodes that have <50 reads in the summed controls: ", len(df3.index))

# df4 = df.drop(df[df.sum_controls < 100].index)
# print("drop barcodes that have <100 reads in the summed controls: ", len(df4.index))

####################################################################################
# ### store the retained/filtered seq.barcodes from filtered df0 into a set
# barcodes_filtered = set(df0.iloc[:,0].tolist())

# ## store barcodes from run_211011 into another set
# for fname in glob.glob("Run211011_raw data.csv"):
#     print(fname)
#     df_211011 = pd.read_csv(fname, sep=';')
#     print (df_211011) # [311688 rows x 49 columns]
# barcodes_211011 = df_211011.iloc[:,0].tolist()
# print(len(barcodes_211011))
# barcodes_211011=list(set(barcodes_211011))
# print(len(barcodes_211011))
# barcodes_211011 = set(barcodes_211011)

# intersection_filtered_210929_and_211011 = barcodes_210929_filtered.intersection(barcodes_211011)
# print(len(intersection_filtered_210929_and_211011))
##########################################################################################







stop = timeit.default_timer()
print(stop - start)  

density_plots.py 
#!/usr/bin/python
# -*- coding: utf-8 -*-
# jupyter notebook --browser='C:/Program Files (x86)/Google/Chrome/Application/chrome.exe %s'

import timeit
start = timeit.default_timer()


import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# tips = sns.load_dataset("tips")
# print(type(tips))
# print(tips)
 
df = pd.read_csv("barcodes_in_controls_stdev_4.csv", sep='	', header=0, index_col=0)
print(df)

# df["(max-min)/avg"].plot.kde()
# plt.savefig("variability_per_barcode_within_controls_across_XPs.pdf")

# df.groupby('experiment')['(max-min)/avg'].plot.density(legend=True)
# plt.title('Density Plots for # of reads variability (max-min)/avg of barcodes in controls')
# plt.savefig("variability_per_barcode_within_controls_across_XPs_density_histograms.pdf")

# bw_methodstr, scalar or callable, optional
# The method used to calculate the estimator bandwidth. This can be ‘scott’, ‘silverman’, a scalar constant or a callable. 
# If None (default), ‘scott’ is used. 
# See scipy.stats.gaussian_kde for more information.
# https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.gaussian_kde.html#scipy.stats.gaussian_kde

#### each experiment #####
# ax = sns.violinplot(x = 'experiment', y = '(max-min)/avg', data = df.reset_index(),
# 	# scale="count", inner="stick",
# 	bw=0.05)
# ax.set_xticklabels(ax.get_xticklabels(),rotation = 45,Fontsize=8)
# ax.figure.tight_layout()
# ax.set_title('Violin Plots for barcode variability (max-min)/avg across controls
(8 experiments, bw=0.05)', pad=5)

#### all experiments #####
ax = sns.violinplot(x = '(max-min)/avg', data = df.reset_index(),
	# scale="count", inner="stick",
	bw=0.05)
ax.set_xticklabels(ax.get_xticklabels(),rotation = 45,Fontsize=8)
ax.figure.tight_layout()
# ax.set_title('Violin Plots for barcode variability (max-min)/avg across controls
(all 8 experiments included, bw=0.05)', pad=5)


# plt.show()
plt.savefig("violin_plot_barcode_variability_all_xp_bw0.05.pdf")


stop = timeit.default_timer()
print(stop - start) 

drop_junk_barcodes.py 
#!/usr/bin/python
# -*- coding: utf-8 -*-
# jupyter notebook --browser='C:/Program Files (x86)/Google/Chrome/Application/chrome.exe %s'


import timeit
start = timeit.default_timer()

import pandas as pd
import glob

with open("/home/nina/Projects/2022_Barcodes/NetBioMed2022/data/combined_runs_filtered.csv", 'a+') as file_write:
    with open("/home/nina/Projects/2022_Barcodes/NetBioMed2022/data/combined_runs.csv") as file_read:
        data = file_read.readlines()
        # get indexes of the control columns
        file_write.write(data[0])
        for line in data[:1] :
            line = line.replace("
","").split(";")
            indices_control, indices_zeros = [], []
            for i in range(len(line)):
                if 'Contro' in line[i]:
                    indices_control.append(i)
                if 'Temps' in line[i]:
                    indices_zeros.append(i)
            print(len(indices_control), len(indices_zeros))
        for line in data[1:] :
            line = line.replace("
","").split(";")
            # reads_control = [line[i] for i in indices_control if float(line[i]) >= 5]
            # reads_zeros = [line[i] for i in indices_zeros if float(line[i]) >= 5]
            line_fill = [ '0' if i == "" else i for i in line]          
            reads_control = [line_fill[i] for i in indices_control if float(line_fill[i]) >= 1]
            reads_zeros = [line_fill[i] for i in indices_zeros if float(line_fill[i]) >= 1]
            if (len(reads_control) >= 5 and len(reads_zeros) >= 5) :
            # if (len(reads_control) == 32 and len(reads_zeros) == 14) :
                # print('OK', reads)
                file_write.write(';'.join(line_fill))
                file_write.write('
')
            # else :
            #     print('KO', reads)



stop = timeit.default_timer()
print(stop - start)  

drop_junk_barcodes_2023.py 
#!/usr/bin/python
# -*- coding: utf-8 -*-
# jupyter notebook --browser='C:/Program Files (x86)/Google/Chrome/Application/chrome.exe %s'


import timeit
start = timeit.default_timer()

import pandas as pd
import glob

with open("A:\Downloads\Projects\workFromHome\Projects\drug_screening\2022_Barcodes\data\combined_runs_filtered.csv", 'a+') as file_write:
    with open("A:\Downloads\Projects\workFromHome\Projects\drug_screening\2022_Barcodes\data\combined_runs.csv") as file_read:
        data = file_read.readlines()
        # get indexes of the control columns
        file_write.write(data[0])
        for line in data[:1] :
            line = line.replace("
","").split(";")
            indices_control, indices_zeros = [], []
            for i in range(len(line)):
                if 'Contro' in line[i]:
                    indices_control.append(i)
                if 'Temps' in line[i]:
                    indices_zeros.append(i)
            print(len(indices_control), len(indices_zeros))
        for line in data[1:] :
            line = line.replace("
","").split(";")
            # reads_control = [line[i] for i in indices_control if float(line[i]) >= 5]
            # reads_zeros = [line[i] for i in indices_zeros if float(line[i]) >= 5]
            line_fill = [ '0' if i == "" else i for i in line]          
            reads_control = [line_fill[i] for i in indices_control if float(line_fill[i]) >= 1]
            reads_zeros = [line_fill[i] for i in indices_zeros if float(line_fill[i]) >= 1]
            if (len(reads_control) >= 5 and len(reads_zeros) >= 5) :
            # if (len(reads_control) == 32 and len(reads_zeros) == 14) :
                # print('OK', reads)
                file_write.write(';'.join(line_fill))
                file_write.write('
')
            # else :
            #     print('KO', reads)



stop = timeit.default_timer()
print(stop - start)  

extract_true_barcode_list_2023.py 
#!/usr/bin/python
# -*- coding: utf-8 -*-
# jupyter notebook --browser='C:/Program Files (x86)/Google/Chrome/Application/chrome.exe %s'


import timeit
start = timeit.default_timer()

import pandas as pd
import glob

with open("A:\Downloads\Projects\workFromHome\Projects\drug_screening\2022_Barcodes\data\combined_runs_filtered_avoid_zero_reads_scaled.csv") as file_read:
	with open("A:\Downloads\Projects\workFromHome\Projects\drug_screening\2022_Barcodes\data\true_barcode_list.txt", "a+") as file_write:
		file_write.write("lentibarcode_id	lentibarcode_sequence
")
		data = file_read.readlines()
		i = 1
		for line in data[1:] :
			line = line.replace("
","").split(";")
			barcode = line[0]
			file_write.write("%s	%s
" % (i, barcode))
			i = i+1

stop = timeit.default_timer()
print(stop - start)   

filter_deseq2_results_2023.py 
#!/usr/bin/python
# -*- coding: utf-8 -*-
# jupyter notebook --browser='C:/Program Files (x86)/Google/Chrome/Application/chrome.exe %s'


import timeit
start = timeit.default_timer()

import pandas as pd
import glob

# for fname in glob.glob("condition_*.csv"):
for fname in glob.glob("condition_*.csv"):
    print(fname)
    # df = pd.read_csv(fname, sep=';', header=0, index_col=0)
    df = pd.read_csv(fname, sep=';', header=0, index_col=0, decimal=',')
    # df['pvalue'] = df['pvalue'].astype(float)
    # df['log2FoldChange'] = df['log2FoldChange'].astype(float)
    # filtered_df = df[df["pvalue"] < 0.05]
    # filtered_df = df[df["pvalue"] < 0.1]
    # filtered_df = filtered_df[filtered_df["log2FoldChange"] > 1]
    # filtered_df = filtered_df[(filtered_df["log2FoldChange"] < -1) | (filtered_df["log2FoldChange"] > 1)]
    # filtered_df.to_csv('pval_log2fc_filtered_%s' % fname[10:], sep=';', index=True)
    # filtered_df.to_csv('pval01/pval_filtered01_%s' % fname[10:], sep=';', index=True)
    # filtered_df.to_csv('pval_filtered\filtered_%s' % fname, sep=';', index=True)
    # filtered_df.to_csv('pval_log2fc_filtered_pos\%s' % fname[14:], sep=';', index=True)

    df.to_csv('unfiltered_%s' % fname[10:], sep=';', index=True)




stop = timeit.default_timer()
print(stop - start)  

filter_desq2_results.py 
#!/usr/bin/python
# -*- coding: utf-8 -*-
# jupyter notebook --browser='C:/Program Files (x86)/Google/Chrome/Application/chrome.exe %s'


import timeit
start = timeit.default_timer()

import pandas as pd
import glob

# for fname in glob.glob("condition_*.csv"):
for fname in glob.glob("condition_*.csv"):
    print(fname)
    # df = pd.read_csv(fname, sep=',', header=0, index_col=0)
    df = pd.read_csv(fname, sep=';', header=0, index_col=0)
    filtered_df = df[df["pvalue"] < 0.05]
    # filtered_df = filtered_df[filtered_df["log2FoldChange"] > 1]
    filtered_df = filtered_df[(filtered_df["log2FoldChange"] < -1) | (filtered_df["log2FoldChange"] > 1)]
    filtered_df.to_csv('pval_log2fc_filtered_%s' % fname[14:], sep=';', index=True)
    # filtered_df.to_csv('pval_filtered\filtered_%s' % fname, sep=';', index=True)
    # filtered_df.to_csv('pval_log2fc_filtered_pos\%s' % fname[14:], sep=';', index=True)





stop = timeit.default_timer()
print(stop - start)  

fold_change.py 
#!/usr/bin/python
# -*- coding: utf-8 -*-
# jupyter notebook --browser='C:/Program Files (x86)/Google/Chrome/Application/chrome.exe %s'

import timeit
start = timeit.default_timer()

from scipy.cluster.hierarchy import linkage, dendrogram, fcluster
from scipy.spatial.distance import squareform
import glob
import statistics


import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
 

experiments = ['exp200921','exp130921','exp300821','exp040821','exp181021','exp151121','exp271221','exp010821']
with open("fold_change_result6_fillna_control_renamed_filtered4.csv", 'a+') as file_write :
	with open("result6_fillna_control_renamed_filtered4.csv") as file_read:
		data = file_read.readlines()
		file_write.write(data[0])
		for experiment in experiments :
			print(experiment)
			barcodes_dict={}
			controls_indexes = []
			for line in data[:1] :
				line = line.replace("
","").split(";")
				for condition in line :
					if ("Contro" in condition) and (experiment in condition) :
						# print(condition, line.index(condition))
						controls_indexes.append(int(line.index(condition)))
			for line in data[1:] :
				line = line.replace("
","").split(";")
				controls_reads = []
				for index in controls_indexes :
					controls_reads.append(float(line[index]))
				controls_avg = statistics.mean(controls_reads)
				# print(controls_reads, controls_avg)
				
				barcode = line[0]
				for read in line[1:] :
					if barcode not in barcodes_dict :
						if controls_avg != 0 :
							barcodes_dict[barcode] = [(float(read)/controls_avg)*100]
						else :
							barcodes_dict[barcode] = ['NaN']
					else :
						if controls_avg != 0 :
							barcodes_dict[barcode].append((float(read)/controls_avg)*100)
						else :
							barcodes_dict[barcode].append('NaN')
			# print(barcodes_dict)	
	for barcode in barcodes_dict :
		file_write.write("%s;%s
" % (barcode,';'.join([str(i) for i in barcodes_dict[barcode]])))			


stop = timeit.default_timer()
print(stop - start) 

fold_change_averaged_over_replicates.py 
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
 

# experiments = ['exp200921','exp130921','exp300821','exp040821','exp181021','exp151121','exp271221','exp010821']
# experiments = ['exp200921']

### we use a list of experiments that go in the same order as in the original file for simplicity ie : 
experiments = ['exp010821', 'exp040821', 'exp300821', 'exp130921', 'exp200921', 'exp181021', 'exp151121','exp271221',]

# with open("fold_change_result6_fillna_control_renamed_filtered4.csv", 'a+') as file_write :
# with open("fold_change_result6_fillna_control_renamed_filtered4_averaged_over_repliactes.csv", 'a+') as file_write :


### we build as many files (fc) as there are experiments, and we will agreggate them afterwards
with open("result6_fillna_control_renamed_filtered4.csv") as file_read:
	data = file_read.readlines()
	conditions_list = []
	barcodes_list = []
	## get conditions list in the same order as in the original data
	for line in data[:1] :
		line = line.replace("
","").split(";")
		for condition in line:
			conditions_list.append(condition)
	print(conditions_list)
	
	# for experiment in experiments : ## we are now looping the experiments, but not necessarily in the same order as in the original file
	for experiment in experiments : ## we are now looping the experiments, in the same order as in the original file
		print(experiment)
		with open("fold_change_result6_fillna_control_renamed_filtered4_%s.csv" % experiment, 'a+') as file_write :
			barcodes_dict={}
			controls_indexes = []
			conditions_indexes = []
			conditions_list = []
			### for each barcode, we want to get the number of reads in the controls of that experiment
			### we first keep in memory the positions of the corresponding columns
			for line in data[:1] :
				line = line.replace("
","").split(";")
				for condition in line :
					if ("Contro" not in condition) and (experiment in condition) :
						conditions_indexes.append(int(line.index(condition)))
						conditions_list.append(condition)
					if ("Contro" in condition) and (experiment in condition) :
						# print(condition, line.index(condition))
						controls_indexes.append(int(line.index(condition)))
			print(controls_indexes)
			print(conditions_indexes, len(conditions_indexes))
			file_write.write(";%s
" % (';'.join(conditions_list)))

			### then, we go over each barcode, and calculate the average of the number of reads within the controls
			for line in data[1:] :
				line = line.replace("
","").split(";")
				barcode = line[0]
				controls_reads = []
				
				for index in controls_indexes :
					controls_reads.append(float(line[index]))
				controls_avg = statistics.mean(controls_reads)
				# print(controls_reads, controls_avg)
				
				for index in conditions_indexes :
					read = line[index]
					if barcode not in barcodes_dict :
						if controls_avg != 0 :
							barcodes_dict[barcode] = [(float(read)/controls_avg)*100]
						else :
							barcodes_dict[barcode] = ['NaN']
					else :
						if controls_avg != 0 :
							barcodes_dict[barcode].append((float(read)/controls_avg)*100)
						else :
							barcodes_dict[barcode].append('NaN')
			# print(barcodes_dict)	


			for barcode in barcodes_dict :
				file_write.write("%s;%s
" % (barcode,';'.join([str(i) for i in barcodes_dict[barcode]])))		
						
	













	# print(len(barcodes_list))
	# print(len(list(set(barcodes_list)))) ## OK same numbers
	
	


	# df_dict = barcodes_dict
	# l = list(df_dict.keys())
	# print(len(l))
	# npl = np.array(l) 
	# nplT = npl.T         
	# doc_df = pd.DataFrame(df_dict, index=[nplT[0]], columns = [nplT[1]])
	# print(doc_df)


	# for (barcode_from_dict,condition) in barcodes_dict :
	# 	for barcode_from_list in barcodes_list :
	# 		if barcode_from_list == barcode_from_dict :
	# 			file_write.write("%s;%s
" % (barcode,';'.join([str(i) for i in barcodes_dict[barcode]])))			


	# for barcode in barcodes_dict :
		# file_write.write("%s;%s
" % (barcode,';'.join([str(i) for i in barcodes_dict[barcode]])))			


stop = timeit.default_timer()
print(stop - start) 

fold_change_experiments_separately.py 
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
 

# experiments = ['exp200921','exp130921','exp300821','exp040821','exp181021','exp151121','exp271221','exp010821']
# experiments = ['exp200921']

### we use a list of experiments that go in the same order as in the original file for simplicity ie : 
experiments = ['exp010821', 'exp040821', 'exp300821', 'exp130921', 'exp200921', 'exp181021', 'exp151121','exp271221',]

# with open("fold_change_result6_fillna_control_renamed_filtered4.csv", 'a+') as file_write :
# with open("fold_change_result6_fillna_control_renamed_filtered4_averaged_over_repliactes.csv", 'a+') as file_write :


### we build as many files (fc) as there are experiments, and we will agreggate them afterwards
with open("result6_fillna_control_renamed_filtered4.csv") as file_read:
	data = file_read.readlines()
	conditions_list = []
	barcodes_list = []
	## get conditions list in the same order as in the original data
	for line in data[:1] :
		line = line.replace("
","").split(";")
		for condition in line:
			conditions_list.append(condition)
	print(conditions_list)
	
	# for experiment in experiments : ## we are now looping the experiments, but not necessarily in the same order as in the original file
	for experiment in experiments : ## we are now looping the experiments, in the same order as in the original file
		print(experiment)
		with open("fold_change_result6_fillna_control_renamed_filtered4_%s.csv" % experiment, 'a+') as file_write :
			barcodes_dict={}
			controls_indexes = []
			conditions_indexes = []
			conditions_list = []
			### for each barcode, we want to get the number of reads in the controls of that experiment
			### we first keep in memory the positions of the corresponding columns
			for line in data[:1] :
				line = line.replace("
","").split(";")
				for condition in line :
					if ("Contro" not in condition) and (experiment in condition) :
						conditions_indexes.append(int(line.index(condition)))
						conditions_list.append(condition)
					if ("Contro" in condition) and (experiment in condition) :
						# print(condition, line.index(condition))
						controls_indexes.append(int(line.index(condition)))
			print(controls_indexes)
			print(conditions_indexes, len(conditions_indexes))
			file_write.write(";%s
" % (';'.join(conditions_list)))

			### then, we go over each barcode, and calculate the average of the number of reads within the controls
			for line in data[1:] :
				line = line.replace("
","").split(";")
				barcode = line[0]
				controls_reads = []
				
				for index in controls_indexes :
					controls_reads.append(float(line[index]))
				controls_avg = statistics.mean(controls_reads)
				# print(controls_reads, controls_avg)
				
				for index in conditions_indexes :
					read = line[index]
					if barcode not in barcodes_dict :
						if controls_avg != 0 :
							barcodes_dict[barcode] = [(float(read)/controls_avg)*100]
						else :
							barcodes_dict[barcode] = ['NaN']
					else :
						if controls_avg != 0 :
							barcodes_dict[barcode].append((float(read)/controls_avg)*100)
						else :
							barcodes_dict[barcode].append('NaN')
			# print(barcodes_dict)	


			for barcode in barcodes_dict :
				file_write.write("%s;%s
" % (barcode,';'.join([str(i) for i in barcodes_dict[barcode]])))		
						

stop = timeit.default_timer()
print(stop - start) 

fold_changes_dendro.py 
#!/usr/bin/python
# -*- coding: utf-8 -*-
# jupyter notebook --browser='C:/Program Files (x86)/Google/Chrome/Application/chrome.exe %s'

import timeit
start = timeit.default_timer()

from scipy.cluster.hierarchy import linkage, dendrogram, fcluster
from scipy.spatial.distance import squareform


import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
 

df = pd.read_csv("fold_changes_agreggated.csv", sep=';', header=0, index_col=0)
print (df) #

matrix = df.corr(
    method = 'pearson',  # The method of correlation
    min_periods = 1).round(2)     # Min number of observations required)

matrix.to_csv('fc_correlations.csv', sep='	', index=True)

# ##################  DENDROGRAM ########################################
# plt.figure(figsize=(12,5))
# dissimilarity = 1 - abs(matrix)
# Z = linkage(squareform(dissimilarity), 'complete')

# with plt.rc_context({'lines.linewidth': 0.5}):
# 	dendrogram(Z, labels=df.columns, orientation='top', leaf_rotation=90, leaf_font_size = 1)

# plt.savefig("dendro.pdf", format='pdf')
# # plt.show()
# ########################################################################

# ############# CORRPLOT + DENDROGRAM ####################################
# # sns.clustermap(matrix, method="complete", cmap='RdBu', annot=True, 
# #                annot_kws={"size": 7}, vmin=-1, vmax=1, figsize=(15,12));
# # plt.show()
# ########################################################################

stop = timeit.default_timer()
print(stop - start) 

fold_changes_join_data.py 
#!/usr/bin/python
# -*- coding: utf-8 -*-
# jupyter notebook --browser='C:/Program Files (x86)/Google/Chrome/Application/chrome.exe %s'


import timeit
start = timeit.default_timer()

import pandas as pd
import glob
df_list = []

for fname in glob.glob("fold_change_result6_fillna_control_renamed_filtered4_*.csv"):
    print(fname)
    df = pd.read_csv(fname, sep=';', header=0, index_col=0)
    print (df) #
    df_list.append(df)

print(len(df_list))

# for i in range(0,len(df_list)-1) :
# 	print(i) 
# 	result = pd.concat([df_list[i], df_list[i+1]], axis=1)

result = pd.concat(df_list, axis=1)

result.to_csv('fold_changes_agreggated.csv', sep=';', index=True)










stop = timeit.default_timer()
print(stop - start)  

heatmap_all_barcodes.py 
#!/usr/bin/python
# -*- coding: utf-8 -*-
# jupyter notebook --browser='C:/Program Files (x86)/Google/Chrome/Application/chrome.exe %s'

import timeit
start = timeit.default_timer()

import scipy.stats as stats 
import numpy as np
from statistics import mean
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# df = pd.read_csv("result6_fillna_control_renamed_filtered4_millions.csv", sep=';', header=0, index_col=0)



experiments = ['exp010821', 'exp040821', 'exp300821', 'exp130921', 'exp200921', 'exp181021', 'exp151121','exp271221',]

for exp in experiments :
    plt.figure()
    df = pd.read_csv("result6_fillna_control_renamed_filtered4_millions.csv", sep=';', header=0, index_col=0)

    # df.drop([col for col in df.columns if 'exp010821' not in col],axis=1,inplace=True)
    df.drop([col for col in df.columns if exp not in col],axis=1,inplace=True)
    df = np.log2(df)

    df.sort_values(by=list(df.columns), axis=0, ascending=True, inplace=False, kind='quicksort', na_position='last')
    print(df)

    xticks = df.columns

    # sns.color_palette("Blues", as_cmap=True)
    ax = plt.axes()

    res = sns.heatmap(df, annot=False,yticklabels=False, xticklabels=xticks, cmap="Blues", ax=ax)
    ax.set_title(exp)
    res.set_xticklabels(res.get_xticklabels(), fontsize = 4, rotation=90)
    res.figure.tight_layout()


    # plt.savefig("result6_fillna_control_renamed_filtered4_heatmap_exp010821.pdf")
    plt.savefig("result6_fillna_control_renamed_filtered4_millions_heatmap_log2_%s.png" % exp)



stop = timeit.default_timer()
print(stop - start) 

heatmap_from_fold_changes.py 
#!/usr/bin/python
# -*- coding: utf-8 -*-
# jupyter notebook --browser='C:/Program Files (x86)/Google/Chrome/Application/chrome.exe %s'

import timeit
start = timeit.default_timer()

import scipy.stats as stats 
import numpy as np
from statistics import mean
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


df = pd.read_csv("fold_changes_agreggated.csv", sep=';', header=0, index_col=0)
print(df)


### average the fold changes over the 4 replicates, for each condition...

sns.heatmap(df, annot=True)
plt.savefig("fold_changes_agreggated_heatmap.pdf")


stop = timeit.default_timer()
print(stop - start) 

histograms_plots.py 
#!/usr/bin/python
# -*- coding: utf-8 -*-
# jupyter notebook --browser='C:/Program Files (x86)/Google/Chrome/Application/chrome.exe %s'

import timeit
start = timeit.default_timer()


import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import pylab as pl

 
df = pd.read_csv("barcodes_in_controls_stdev_4.csv", sep='	', header=0, index_col=0)
print(df)
df = df.reset_index()

# axes = df['(max-min)/avg'].hist(by=df['experiment'],bins=100)
axes = df['(max-min)/avg'].hist(bins=100)

# plt.title('Frequency Plots for barcode variability (max-min)/avg 
 across controls in each experiment
(100 bins)')
plt.savefig('histo_plots_all_xps.pdf')

stop = timeit.default_timer()
print(stop - start) 

igraph_networks.py 
#!/usr/bin/python
# -*- coding: utf-8 -*-
# jupyter notebook --browser='C:/Program Files (x86)/Google/Chrome/Application/chrome.exe %s'

import timeit
start = timeit.default_timer()

from scipy.cluster.hierarchy import linkage, dendrogram, fcluster
from scipy.spatial.distance import squareform


import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import networkx as nx
from igraph import *



df1 = pd.read_csv("result6_fillna_control_renamed_filtered4_exp130921_T0_renamed_replicates_averages.csv", sep=';', header=0, index_col=0)

print (df1) #
df1 = df1.astype(float)
df1 = df1.sort_index(axis=1, ascending=True, key=lambda x: x.str.lower())

matrix1 = df1.corr(
    method = 'pearson',  # The method of correlation
    min_periods = 1).round(2)     # Min number of observations required)


 
# Transform it in a links data frame (3 columns only):
links = matrix1.stack().reset_index()
links.columns = ['source', 'target', 'correlation']
 



# Keep only correlation over a threshold and remove self correlation (cor(A,A)=1)
links_filtered=links.loc[ (links['correlation'] > 0.75) & (links['source'] != links['target']) ]
# Build your graph
G = Graph.DataFrame(edges, directed=False, vertices=links['correlation'])




# plt.show()








stop = timeit.default_timer()
print(stop - start) 

join_data.py 
#!/usr/bin/python
# -*- coding: utf-8 -*-
# jupyter notebook --browser='C:/Program Files (x86)/Google/Chrome/Application/chrome.exe %s'


import timeit
start = timeit.default_timer()

import pandas as pd
import glob
path = "A:\Downloads\Projects\workFromHome\Projects\drug_screening\data\*.csv"
df_list = []
run = 1
for fname in glob.glob("result6.csv"):
    print(fname)
    df = pd.read_csv(fname, sep=';', header=0, index_col=0)
    #print (df) #
    df_list.append(df)
    run = run + 1
    if run == 3 :
        break

result = pd.concat(df_list, axis=1)
result.to_csv('result7.csv', sep=';', index=True)










stop = timeit.default_timer()
print(stop - start)  

merge_csv_files.py 
#!/usr/bin/python
# -*- coding: utf-8 -*-
# jupyter notebook --browser='C:/Program Files (x86)/Google/Chrome/Application/chrome.exe %s'


import timeit
start = timeit.default_timer()

import os
import glob
import pandas as pd
os.chdir("A:\Downloads\Projects\workFromHome\Projects\drug_screening\2022_Barcodes\data\")
extension = 'csv'
all_filenames = [i for i in glob.glob('Run*.{}'.format(extension))]
#combine all files in the list
combined_csv = pd.concat([pd.read_csv(f, sep=";", index_col=0) for f in all_filenames ], axis=1)
#export to csv
combined_csv.to_csv( "combined_runs.csv", index=True, encoding='utf-8-sig', sep=";")



# path = "A:\Downloads\Projects\workFromHome\Projects\drug_screening\2022_Barcodes\data\combine*.csv"

# for fname in glob.glob(path):
#     print(fname)
#     df = pd.read_csv(fname, sep=',', header=0, index_col=0)
#     print (len(df.columns)) #









stop = timeit.default_timer()
print(stop - start)  

merge_deseq2_results_2023.py 
#!/usr/bin/python
# -*- coding: utf-8 -*-
# jupyter notebook --browser='C:/Program Files (x86)/Google/Chrome/Application/chrome.exe %s'


import timeit
start = timeit.default_timer()

import pandas as pd
import glob

colnames = []
# for fname in glob.glob("condition_*.csv"):
# for fname in glob.glob("pval_log2fc_filtered_*.csv"):
for fname in glob.glob("unfiltered_*.csv"):

    print(fname)
    colname = '_'.join(fname.split('_')[1:3])
    colnames.append(colname)
    print(colname)

merged_df = pd.DataFrame(columns = colnames)
# print(merged_df)


# # for fname in glob.glob("condition_*.csv"):
# for fname in glob.glob("pval_log2fc_filtered_*.csv"):
for fname in glob.glob("unfiltered_*.csv"):
    print(fname)  
    colname = '_'.join(fname.split('_')[1:3])
    df = pd.read_csv(fname, sep=';', header=0, index_col=0)
    logfc = df['log2FoldChange']
    merged_df[colname] = logfc

print(merged_df)

#df2 = merged_df.T.drop_duplicates().T

# merged_df = merged_df.dropna()

# print(merged_df)

#print(df2)


# df2.to_csv('merged_logfc_pval_filtered_deseq2_2023.csv', sep=';', index=True)
# merged_df.to_csv('merged_logfc_pval_filtered_deseq2_2023.csv', sep=';', index=True)
merged_df.to_csv('merged_unfiltered_deseq2_2023.csv', sep=';', index=True)





stop = timeit.default_timer()
print(stop - start)  

merge_deseq2_results_with_expdate_2023.py 
#!/usr/bin/python
# -*- coding: utf-8 -*-
# jupyter notebook --browser='C:/Program Files (x86)/Google/Chrome/Application/chrome.exe %s'


import timeit
start = timeit.default_timer()

import pandas as pd
import glob

colnames = []
# for fname in glob.glob("condition_*.csv"):
for fname in glob.glob("pval_log2fc_filtered_*.csv"):
# for fname in glob.glob("unfiltered_*.csv"):

    # print(fname)
    colname = fname.split('_')
    # print(colname)
    # colname_sublist = [colname[1], colname[2], colname[5][0:-4]] # for unfiltered
    colname_sublist = [colname[3], colname[4], colname[7][0:-4]] # for filtered
    # print(colname_sublist)
    colname_joined = '_'.join(colname_sublist)
    colnames.append(colname_joined)
    print(colname_joined)

merged_df = pd.DataFrame(columns = colnames)
# print(merged_df)


# # for fname in glob.glob("condition_*.csv"):
for fname in glob.glob("pval_log2fc_filtered_*.csv"):
# for fname in glob.glob("unfiltered_*.csv"):
    # print(fname)  
    colname = fname.split('_')
    # colname_sublist = [colname[1], colname[2], colname[5][0:-4]] # for unfiltered
    colname_sublist = [colname[3], colname[4], colname[7][0:-4]] # for filtered
    colname_joined = '_'.join(colname_sublist)
    df = pd.read_csv(fname, sep=';', header=0, index_col=0)
    logfc = df['log2FoldChange']
    merged_df[colname_joined] = logfc

print(merged_df)

#df2 = merged_df.T.drop_duplicates().T

# merged_df = merged_df.dropna()

# print(merged_df)

#print(df2)


#merged_df.to_csv('merged_logfc_pval_filtered_deseq2_with_expdate_2023.csv', sep=';', index=True)
# merged_df.to_csv('merged_unfiltered_deseq2_with_expdate_2023.csv', sep=';', index=True)





stop = timeit.default_timer()
print(stop - start)  

multiply_counts_by_100.py 
#!/usr/bin/python
# -*- coding: utf-8 -*-
# jupyter notebook --browser='C:/Program Files (x86)/Google/Chrome/Application/chrome.exe %s'


import timeit
start = timeit.default_timer()

import pandas as pd
import glob

with open("A:\Downloads\Projects\workFromHome\Projects\drug_screening\2022_Barcodes\data\combined_runs_filtered_avoid_zero_reads_all_replaced.csv") as file_read:
    df = pd.read_csv(file_read, sep=';', header=0, index_col=0)
    print(df)
    df = df.fillna(0).apply(lambda x: 100 * x)

    df.to_csv('A:\Downloads\Projects\workFromHome\Projects\drug_screening\2022_Barcodes\data\combined_runs_filtered_avoid_zero_reads_all_replaced_x100.csv', sep=';', index=True, float_format='%.0f')

stop = timeit.default_timer()
print(stop - start)   

networks.py 
#!/usr/bin/python
# -*- coding: utf-8 -*-
# jupyter notebook --browser='C:/Program Files (x86)/Google/Chrome/Application/chrome.exe %s'

import timeit
start = timeit.default_timer()

from scipy.cluster.hierarchy import linkage, dendrogram, fcluster
from scipy.spatial.distance import squareform


import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import networkx as nx



matrix1 = pd.read_csv("correl_matrix1_merged_pval_filtered_deseq2.csv", sep=';', header=0, index_col=0)


 
# Transform it in a links data frame (3 columns only):
links = matrix1.stack().reset_index()
links.columns = ['source', 'target', 'correlation']
 



# Keep only correlation over a threshold and remove self correlation (cor(A,A)=1)
# links_filtered=links.loc[ (links['correlation'] > 0) & (links['source'] != links['target']) ]
# links_filtered=links.loc[ (links['source'] != links['target']) ]
links_filtered=links.loc[ (links['correlation'] > 0) & (links['source'] != links['target']) ]
# Build your graph
G=nx.from_pandas_edgelist(links_filtered, 'source', 'target', 'correlation')
# G=nx.from_pandas_edgelist(links, 'source', 'target', 'correlation')
my_nodes = list(G.nodes)
print(len(my_nodes))
# print(my_nodes)



mapping = {}
for node in my_nodes :
	print(node)
	mapping[node] = '_'.join(node.split("_")[0:2])
	print(mapping[node])

print(len(mapping))

G = nx.relabel_nodes(G, mapping)

# pos = nx.circular_layout(G, scale = 2, dim = 2)#, seed=7)  # positions for all nodes - seed for reproducibility
pos = nx.spring_layout(G, k=0.1)#, scale=2)

# pos = nx.random_layout(G)

np.random.shuffle(my_nodes)

# nodes

color_dict = {}
with open("color_mapping.tsv", 'r') as file_read :
	data = file_read.readlines()
	for line in data :
		line = line.replace("
", "").split("	")
		color = int(line[0])
		if line[1] == 'Fluor_006u' :
			drug = '5Fluor_006u'
		elif line[1] == 'Azacyt_1,5u' :
			drug = 'Azacyt_1.5u'
		elif line[1] == 'Bafilo_1,2n' :
			drug = 'Bafilo_1.2n'
		else :
			drug = line[1]
		color_dict[drug] = color
print(color_dict)
colors = []
for node in G.nodes :
	colors.append(((color_dict[node]+1)/100))
print(colors)

# nx.draw_networkx_nodes(G, pos, node_size=100, node_color='yellow', alpha = 0.5)
nx.draw_networkx_nodes(G, pos, node_size=300, node_color=colors)#, alpha = 0.5)

# edges
elarge = [(u, v) for (u, v, d) in G.edges(data=True) if d['correlation'] > 0.8]
# esmall = [(u, v) for (u, v, d) in G.edges(data=True) if 0.5 < d['correlation'] <= 0.9]
nx.draw_networkx_edges(G, pos, edgelist=elarge, width=0.1, edge_color='black')#, alpha=0.5)
# nx.draw_networkx_edges(G, pos, edgelist=esmall, width=0.5, alpha=0.5, edge_color="b", style="dashed")

# labels

label_dict = {}
for node in G.nodes :
	label = node.split("_")[0]
	label_dict[node] = label
print(label_dict)
nx.draw_networkx_labels(G, pos, labels = label_dict, font_size=4, font_family="sans-serif")





 
# # Plot the network:
# # nx.draw(G, pos=nx.circular_layout(G), with_labels=True, node_color='orange', node_size=25, edge_color='pink', linewidths=1, font_size=6)
# # nx.draw(G, pos=nx.fruchterman_reingold_layout(G), with_labels=True, node_color='orange', node_size=250, edge_color='pink', linewidths=1, font_size=6)
# nx.draw(G, pos=nx.spring_layout(G), with_labels=True, node_color='orange', node_size=250, edge_color='pink', linewidths=1, font_size=6)

# plt.show()

plt.savefig("network_correlation_colored.pdf", format='pdf')
plt.savefig("network_correlation_colored.png", format='png')








# df2 = pd.read_csv("result6_fillna_control_renamed_filtered4_exp130921_T0_renamed_replicates_averages_fc_zeros.csv", sep=';', header=0, index_col=0)
# print (df2) #
# df2 = df2.astype(float)
# df2 = df2.sort_index(axis=1, ascending=True, key=lambda x: x.str.lower())

# matrix2 = df2.corr(
#     method = 'pearson',  # The method of correlation
#     min_periods = 1).round(2)     # Min number of observations required)



stop = timeit.default_timer()
print(stop - start) 

normalize_data_within_samples.py 
#!/usr/bin/python
# -*- coding: utf-8 -*-
# jupyter notebook --browser='C:/Program Files (x86)/Google/Chrome/Application/chrome.exe %s'


import timeit
start = timeit.default_timer()

import pandas as pd
import glob
# for fname in glob.glob("*.csv"):
for fname in glob.glob("result6_fillna_control_renamed_filtered4.csv"):
    print(fname)
    df = pd.read_csv(fname, sep=';', header=0, index_col=0)
    print ("df") #
    print (df) #

    percents_df = df.apply(lambda x: 100 * x / float(x.sum()))
    print("percents_df")
    print(percents_df)

    percents_df.to_csv('result6_fillna_control_renamed_filtered4_percents.csv', sep=';', index=True)

    



stop = timeit.default_timer()
print(stop - start)  

normalize_data_within_samples_millions.py 
#!/usr/bin/python
# -*- coding: utf-8 -*-
# jupyter notebook --browser='C:/Program Files (x86)/Google/Chrome/Application/chrome.exe %s'


import timeit
start = timeit.default_timer()

import pandas as pd
import glob
# for fname in glob.glob("*.csv"):
for fname in glob.glob("result6_fillna_control_renamed_filtered4.csv"):
    print(fname)
    df = pd.read_csv(fname, sep=';', header=0, index_col=0)
    print ("df") #
    print (df) #

    percents_df = df.apply(lambda x: 1000000 * x / float(x.sum()))
    print("percents_df")
    print(percents_df)

    percents_df.to_csv('result6_fillna_control_renamed_filtered4_millions.csv', sep=';', index=True)

    



stop = timeit.default_timer()
print(stop - start)  

normalize_per_millions_2023.py 
#!/usr/bin/python
# -*- coding: utf-8 -*-
# jupyter notebook --browser='C:/Program Files (x86)/Google/Chrome/Application/chrome.exe %s'


import timeit
start = timeit.default_timer()

import pandas as pd
import glob

with open("A:\Downloads\Projects\workFromHome\Projects\drug_screening\2022_Barcodes\data\combined_runs_filtered_avoid_zero_reads.csv") as file_read:
    df = pd.read_csv(file_read, sep=';', header=0, index_col=0)
    print(df)
    df = df.fillna(0).apply(lambda x: 1000000 * x / float(x.sum()))

    df.to_csv('A:\Downloads\Projects\workFromHome\Projects\drug_screening\2022_Barcodes\data\combined_runs_filtered_avoid_zero_reads_scaled.csv', sep=';', index=True, float_format='%.6f')

stop = timeit.default_timer()
print(stop - start)   

replace_NaN.py 
#!/usr/bin/python
# -*- coding: utf-8 -*-
# jupyter notebook --browser='C:/Program Files (x86)/Google/Chrome/Application/chrome.exe %s'


import timeit
start = timeit.default_timer()

import pandas as pd
import glob
# path = "A:\Downloads\Projects\workFromHome\Projects\drug_screening\data\joined\*.csv"
# for fname in glob.glob("*.csv"):
for fname in glob.glob("joined\result7.csv"):
    print(fname)
    df = pd.read_csv(fname, sep=';', header=0, index_col=0)
    # print (df)
    df_filled = df.fillna(0)
    df_filled.to_csv('result7_fillna.csv', sep=';', index=True)






stop = timeit.default_timer()
print(stop - start)  

replace_with_commas_2023.py 
#!/usr/bin/python
# -*- coding: utf-8 -*-
# jupyter notebook --browser='C:/Program Files (x86)/Google/Chrome/Application/chrome.exe %s'


import timeit
start = timeit.default_timer()

import pandas as pd
import glob

with open("A:\Downloads\Projects\workFromHome\Projects\drug_screening\2022_Barcodes\data\combined_runs_filtered_avoid_zero_reads_scaled.csv") as file_read:
	with open('A:\Downloads\Projects\workFromHome\Projects\drug_screening\2022_Barcodes\data\combined_runs_filtered_avoid_zero_reads_scaled_commas.csv', 'a+') as file_write:
		data = file_read.readlines()
		file_write.write(data[0])
		for line in data[1:] :
			line = line.replace(".", ",")
			file_write.write(line)

stop = timeit.default_timer()
print(stop - start)   

replicates_handling.py 
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
 

with open("result6_fillna_control_renamed_filtered4.csv") as file_read:
	data = file_read.readlines()
	### parse headers
	drugs_list, exp_list = [], []
	drugs_index_dict, exp_index_dict = {}, {}
	barcode_dict = {}
	index_content_dict = {}
	drug_dict = {}
	exp_dict = {}
	for line in data[:1] :
		line = line.replace("
","").split(";")
		for condition in line[1:] :
			conditions = condition.split("_")
			drug = conditions[0][:-1]
			replicate = conditions[0][-1]
			print(drug,replicate)
			drugs_list.append(drug)
			exp = conditions[2]
			exp_list.append(exp)

			if drug not in drugs_index_dict :
				drugs_index_dict[drug] = [line.index(condition)]
			else :
				drugs_index_dict[drug].append(line.index(condition))
			if exp not in exp_index_dict :
				exp_index_dict[exp] = [line.index(condition)]
			else :
				exp_index_dict[exp].append(line.index(condition))

			index_content_dict[line.index(condition)] = (drug,replicate,exp)

			################# fill the drug_dict ############################
			if drug not in drug_dict :
				drug_dict[drug] = {}
				
			if exp not in drug_dict[drug] :
				drug_dict[drug][exp] = {}
			
			if replicate not in drug_dict[drug][exp] :
				drug_dict[drug][exp][replicate] = [line.index(condition)]
			else :
				drug_dict[drug][exp][replicate].append(line.index(condition))
			#################################################################

			################# fill the exp_dict #############################
			if exp not in exp_dict :
				exp_dict[exp] = {}
				
			if drug not in exp_dict[exp] :
				exp_dict[exp][drug] = {}
			
			if replicate not in exp_dict[exp][drug] :
				exp_dict[exp][drug][replicate] = [line.index(condition)]
			else :
				exp_dict[exp][drug][replicate].append(line.index(condition))			

				


	drugs_list = list(set(drugs_list))
	exp_list = list(set(exp_list))

print("drugs_list", drugs_list)
print("exp_list", exp_list)	

for drug in drugs_index_dict :
	print(drug, drugs_index_dict[drug])

for exp in exp_index_dict :
	print(exp, exp_index_dict[exp])

# for index in index_content_dict :
# 	# print(index, index_content_dict[index])
# 	print(index, index_content_dict[index][0], index_content_dict[index][1], index_content_dict[index][2])

for drug in drug_dict :
	print(drug, drug_dict[drug])

for exp in exp_dict :
	print(exp, exp_dict[exp])

#### build file displaying the reads numbers as fold changes compared to 
with open("result6_fillna_control_renamed_filtered4.csv") as file_read:
	data = file_read.readlines()
	for line in data[1:] :
		line = line.replace("
","").split(";")
		barcode = line[0]
		reads = line[1:]


stop = timeit.default_timer()
print(stop - start) 

swarmplot_reads_variability_within_controls.py 
#!/usr/bin/python
# -*- coding: utf-8 -*-
# jupyter notebook --browser='C:/Program Files (x86)/Google/Chrome/Application/chrome.exe %s'

import timeit
start = timeit.default_timer()


import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# tips = sns.load_dataset("tips")
# print(type(tips))
# print(tips)
 
df = pd.read_csv("barcodes_in_controls_stdev_4.csv", sep='	', header=0, index_col=0)
print(df)

# ax = sns.swarmplot(x="barcode", y="(max-min)/avg", hue="experiment", data=df)
# ax = sns.swarmplot(x=df["(max-min)/avg"])
# plt.show()

sns.set_theme(style="whitegrid")
ax = sns.swarmplot(x=df["(max-min)/avg"])
plt.savefig("test.pdf")

stop = timeit.default_timer()
print(stop - start) 

test_controls_variability.py 
#!/usr/bin/python
# -*- coding: utf-8 -*-
# jupyter notebook --browser='C:/Program Files (x86)/Google/Chrome/Application/chrome.exe %s'

import timeit
start = timeit.default_timer()

from scipy.cluster.hierarchy import linkage, dendrogram, fcluster
from scipy.spatial.distance import squareform
import glob
import statistics


import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
 


with open("result6_fillna_control_renamed_filtered4.csv") as file_read:
	data = file_read.readlines()
	conditions, experiments, runs = [],[],[]
	for line in data[:1] :
		line = line.replace("
","").split(";")
		for condition in line[1:] :
			condition = condition.split("_")
			# print(('	').join(condition))


			conditions.append('_'.join(condition[0:2]))
			experiments.append(condition[2])
			runs.append(condition[3])
conditions = list(set(conditions))
experiments = list(set(experiments))
runs = list(set(runs))

# print("conditions")
# print(conditions)
# print("experiments")
# print(experiments)
# print("runs")
# print(runs)

# print("conditions")
# for i in sorted(conditions,key=str.lower) :
# 	print(i)
# print("experiments")
# for i in experiments :
# 	print(i)
# print("runs")
# for i in runs :
# 	print(i)

with open("barcodes_in_controls_stdev_4.csv", 'a+') as file_write:
	file_write.write("experiment	barcode	cntrl1	cntrl2	cntrl3	cntrl4	avg	stdev	(max-min)/avg
")
	with open("result6_fillna_control_renamed_filtered4.csv") as file_read:
		data = file_read.readlines()
		for experiment in experiments :
			# print(experiment)
			controls_indexes = []
			for line in data[:1] :
				line = line.replace("
","").split(";")
				for condition in line :
					if ("Contro" in condition) and (experiment in condition) :
						# print(condition, line.index(condition))
						controls_indexes.append(int(line.index(condition)))
			for line in data[1:] :
				line = line.replace("
","").split(";")
				controls_reads = []
				for index in controls_indexes :
					controls_reads.append(float(line[index]))
				controls_avg = round(statistics.mean(controls_reads), 2)
				controls_stdev = round(statistics.stdev(controls_reads), 2)
				# avg_std = round(controls_avg / controls_stdev, 2)
				if controls_avg != 0 :
					max_min = round(((max(controls_reads) - min(controls_reads)) / controls_avg), 2)
				else :
					max_min = '0'
				# print(experiment,line[0], '	'.join([str(i) for i in controls_reads]), controls_avg, controls_stdev,avg_std, max_min)
				# print(experiment,line[0], '	'.join([str(i) for i in controls_reads]), controls_avg, controls_stdev, max_min)
				# print(experiment, controls_avg, controls_stdev,avg_std, max_min)
				file_write.write("%s	%s	%s	%s	%s	%s
" % (experiment,line[0], '	'.join([str(i) for i in controls_reads]), controls_avg, controls_stdev, max_min))
stop = timeit.default_timer()
print(stop - start) 

test_controls_variability_from_fold_change_values.py 
#!/usr/bin/python
# -*- coding: utf-8 -*-
# jupyter notebook --browser='C:/Program Files (x86)/Google/Chrome/Application/chrome.exe %s'

import timeit
start = timeit.default_timer()

from scipy.cluster.hierarchy import linkage, dendrogram, fcluster
from scipy.spatial.distance import squareform
import glob
import statistics


import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
 


with open("result6_fillna_control_renamed_filtered4.csv") as file_read:
	data = file_read.readlines()
	conditions, experiments, runs = [],[],[]
	for line in data[:1] :
		line = line.replace("
","").split(";")
		for condition in line[1:] :
			condition = condition.split("_")
			# print(('	').join(condition))


			conditions.append('_'.join(condition[0:2]))
			experiments.append(condition[2])
			runs.append(condition[3])
conditions = list(set(conditions))
experiments = list(set(experiments))
runs = list(set(runs))

# print("conditions")
# print(conditions)
# print("experiments")
# print(experiments)
# print("runs")
# print(runs)

# print("conditions")
# for i in sorted(conditions,key=str.lower) :
# 	print(i)
# print("experiments")
# for i in experiments :
# 	print(i)
# print("runs")
# for i in runs :
# 	print(i)

with open("barcodes_in_controls_stdev_4_from_fold_changes.csv", 'a+') as file_write:
	file_write.write("experiment	barcode	cntrl1	cntrl2	cntrl3	cntrl4	avg	stdev	(max-min)/avg
")
	with open("fold_change_result6_fillna_control_renamed_filtered4.csv") as file_read:
		data = file_read.readlines()
		for experiment in experiments :
			# print(experiment)
			controls_indexes = []
			for line in data[:1] :
				line = line.replace("
","").split(";")
				for condition in line :
					if ("Contro" in condition) and (experiment in condition) :
						# print(condition, line.index(condition))
						controls_indexes.append(int(line.index(condition)))
			for line in data[1:] :
				line = line.replace("
","").split(";")
				controls_reads = []
				for index in controls_indexes :
					controls_reads.append(float(line[index]))
				controls_avg = round(statistics.mean(controls_reads), 2)
				controls_stdev = round(statistics.stdev(controls_reads), 2)
				# avg_std = round(controls_avg / controls_stdev, 2)
				if controls_avg != 0 :
					max_min = round(((max(controls_reads) - min(controls_reads)) / controls_avg), 2)
				else :
					max_min = '0'
				# print(experiment,line[0], '	'.join([str(i) for i in controls_reads]), controls_avg, controls_stdev,avg_std, max_min)
				# print(experiment,line[0], '	'.join([str(i) for i in controls_reads]), controls_avg, controls_stdev, max_min)
				# print(experiment, controls_avg, controls_stdev,avg_std, max_min)
				file_write.write("%s	%s	%s	%s	%s	%s
" % (experiment,line[0], '	'.join([str(i) for i in controls_reads]), controls_avg, controls_stdev, max_min))
stop = timeit.default_timer()
print(stop - start) 

unique_barcodes.py 
#!/usr/bin/python
# -*- coding: utf-8 -*-
# jupyter notebook --browser='C:/Program Files (x86)/Google/Chrome/Application/chrome.exe %s'

import timeit
start = timeit.default_timer()

from scipy.cluster.hierarchy import linkage, dendrogram, fcluster
from scipy.spatial.distance import squareform
import glob


import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
 
# df = pd.read_csv("result6_fillna_control_renamed_filtered6.csv", sep=';', header=0, index_col=0)

controls, time_zeros, compounds = [],[],[]

for fname in glob.glob("result6_fillna_control_renamed_filtered6.csv"):
    print(fname)
    
    df = pd.read_csv(fname, sep=';', header=0, index_col=0)
    print (df.columns) #

    for col in df.columns:
        print(col)
#     for col in df.columns :
#         if 'Contro' in col :
#             print (col)
#             controls.append(col)
#         if 'Temps' in col :
#         	time_zeros.append(col)
#         else :
#         	compounds.append(col)

# print(len(controls))
# print(len(time_zeros))
# print(len(compounds))
# print(set(controls).intersection(set(time_zeros)))







stop = timeit.default_timer()
print(stop - start) 
