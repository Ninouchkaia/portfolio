#!/usr/bin/python
# -*- coding: utf-8 -*-
# jupyter notebook --browser='C:/Program Files (x86)/Google/Chrome/Application/chrome.exe %s'


import timeit
start = timeit.default_timer()

import os
import glob
import pandas as pd
os.chdir("A:\\Downloads\\Projects\\workFromHome\\Projects\\drug_screening\\2022_Barcodes\\data\\")
extension = 'csv'
all_filenames = [i for i in glob.glob('Run*.{}'.format(extension))]
#combine all files in the list
combined_csv = pd.concat([pd.read_csv(f, sep=";", index_col=0) for f in all_filenames ], axis=1)
#export to csv
combined_csv.to_csv( "combined_runs.csv", index=True, encoding='utf-8-sig', sep=";")



# path = "A:\\Downloads\\Projects\\workFromHome\\Projects\\drug_screening\\2022_Barcodes\\data\\combine*.csv"

# for fname in glob.glob(path):
#     print(fname)
#     df = pd.read_csv(fname, sep=',', header=0, index_col=0)
#     print (len(df.columns)) #









stop = timeit.default_timer()
print(stop - start)  