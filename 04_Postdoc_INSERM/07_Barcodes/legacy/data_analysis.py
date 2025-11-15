#!/usr/bin/python
# -*- coding: utf-8 -*-
# jupyter notebook --browser='C:/Program Files (x86)/Google/Chrome/Application/chrome.exe %s'


import timeit
start = timeit.default_timer()

import pandas as pd
import glob
path = "A:\\Downloads\\Projects\\workFromHome\\Projects\\drug_screening\\data\\*.csv"
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