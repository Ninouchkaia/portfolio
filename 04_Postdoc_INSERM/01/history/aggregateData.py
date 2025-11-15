#!/usr/bin/python
# -*- coding: utf-8 -*-
# jupyter notebook --browser='C:/Program Files (x86)/Google/Chrome/Application/chrome.exe %s'


import timeit
start = timeit.default_timer()



import sys
import os
import csv
# with open("outputs_ABM_9_2.txt", 'a+') as file_write:
with open(sys.argv[2], 'a+') as file_write:
	# rootdir = 'A:\\Downloads\\Projects\\workFromHome\\Projects\\ABM2021\\explorations\\ABM_9_2\\ABM_9_2'
	rootdir = sys.argv[1]
	for subdir, dirs, files in os.walk(rootdir):
		for file in files :
			print(file)
			path_to_file = rootdir + "\\" + file
			if file == 'population1.csv' :
				with open(path_to_file, 'r') as file_read:
					data = file_read.readlines()
					for line in data :
						file_write.write(line)
			else :
				with open(path_to_file, 'r') as file_read:
					data = file_read.readlines()
					for line in data[1:] :
						file_write.write(line)

stop = timeit.default_timer()
print(stop - start)  
		


