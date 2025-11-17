#!/usr/bin/python
# -*- coding: utf-8 -*-
# jupyter notebook --browser='C:/Program Files (x86)/Google/Chrome/Application/chrome.exe %s'


import timeit
start = timeit.default_timer()

import pandas as pd
import glob

with open("A:\\Downloads\\Projects\\workFromHome\\Projects\\drug_screening\\2022_Barcodes\\data\\combined_runs_filtered.csv", 'a+') as file_write:
    with open("A:\\Downloads\\Projects\\workFromHome\\Projects\\drug_screening\\2022_Barcodes\\data\\combined_runs.csv") as file_read:
        data = file_read.readlines()
        # get indexes of the control columns
        file_write.write(data[0])
        for line in data[:1] :
            line = line.replace("\n","").split(";")
            indices_control, indices_zeros = [], []
            for i in range(len(line)):
                if 'Contro' in line[i]:
                    indices_control.append(i)
                if 'Temps' in line[i]:
                    indices_zeros.append(i)
            print(len(indices_control), len(indices_zeros))
        for line in data[1:] :
            line = line.replace("\n","").split(";")
            # reads_control = [line[i] for i in indices_control if float(line[i]) >= 5]
            # reads_zeros = [line[i] for i in indices_zeros if float(line[i]) >= 5]
            line_fill = [ '0' if i == "" else i for i in line]          
            reads_control = [line_fill[i] for i in indices_control if float(line_fill[i]) >= 1]
            reads_zeros = [line_fill[i] for i in indices_zeros if float(line_fill[i]) >= 1]
            if (len(reads_control) >= 5 and len(reads_zeros) >= 5) :
            # if (len(reads_control) == 32 and len(reads_zeros) == 14) :
                # print('OK', reads)
                file_write.write(';'.join(line_fill))
                file_write.write('\n')
            # else :
            #     print('KO', reads)



stop = timeit.default_timer()
print(stop - start)  