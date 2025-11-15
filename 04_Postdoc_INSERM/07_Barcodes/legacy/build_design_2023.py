#!/usr/bin/python
# -*- coding: utf-8 -*-
# jupyter notebook --browser='C:/Program Files (x86)/Google/Chrome/Application/chrome.exe %s'

import timeit
start = timeit.default_timer()


 
with open("design_6_2023.csv", 'w') as file_write :
	file_write.write(";run;exp;replicate;condition\n")
	with open("A:\\Downloads\\Projects\\workFromHome\\Projects\\drug_screening\\2022_Barcodes\\data\\combined_runs_filtered_avoid_zero_reads_x100.csv") as file_read:
		data = file_read.readlines()
		for line in data[:1] :
			line = line.replace("\n","").split(";")
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
					file_write.write("control\n")
					# file_write.write("control_%s_%s\n" % (exp,run))
				else :
					file_write.write("%s_%s\n" % (drug,dosage) )
					# file_write.write("%s_%s_%s_%s\n" % (drug,dosage,exp,run) )


			


stop = timeit.default_timer()
print(stop - start) 