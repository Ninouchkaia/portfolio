#!/usr/bin/python
# -*- coding: utf-8 -*-
# jupyter notebook --browser='C:/Program Files (x86)/Google/Chrome/Application/chrome.exe %s'


import timeit
start = timeit.default_timer()

import pandas as pd
import glob
path = "A:\\Downloads\\Projects\\workFromHome\\Projects\\drug_screening\\data\\*.csv"
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