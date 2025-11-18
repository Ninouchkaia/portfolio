#!/usr/bin/python
# -*- coding: utf-8 -*-
# jupyter notebook --browser='C:/Program Files (x86)/Google/Chrome/Application/chrome.exe %s'

import timeit
start = timeit.default_timer()


 
with open("design4.csv", 'w') as file_write :
	file_write.write(";exp;replicate;condition\n")
	with open("A:\\Downloads\\Projects\\workFromHome\\Projects\\drug_screening\\2022_Barcodes\\data\\combined_runs_normalized_filtered.csv") as file_read:
		data = file_read.readlines()
		for line in data[:1] :
			line = line.replace("\n","").split(";")
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
					file_write.write("control\n")
				else :
					file_write.write("%s_%s\n" % (drug,dosage) )


			


stop = timeit.default_timer()
print(stop - start) 