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


