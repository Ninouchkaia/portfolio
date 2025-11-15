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
