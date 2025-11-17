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
			line = line.replace("\n","").split(';')
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
			file_write.write("%s;%s\n" % (barcode, ';'.join(fc_class_list)))



stop = timeit.default_timer()
print(stop - start) 