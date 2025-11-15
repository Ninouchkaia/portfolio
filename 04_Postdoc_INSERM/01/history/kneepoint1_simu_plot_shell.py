#!/usr/bin/python
# -*- coding: utf-8 -*-
# jupyter notebook --browser='C:/Program Files (x86)/Google/Chrome/Application/chrome.exe %s'


import timeit
start = timeit.default_timer()

patients_list_with_mono = ['CRE1704-1.1%','LAU1405-2.5%','DES2105-1.25%','ORE1706-0.68%','CAS1802-1.04%', 'GER160522-0.45%','REI230522-0.95%','PUJ240522-0.34%','LAR300522-0.21%','CAZ310522-3.48%']
patients_list = []

for patient in patients_list_with_mono :
	patient_name = patient.split("-")[0]
	patients_list.append(patient_name)
print(patients_list)

with open("plot_kneepoint1_simu_by_patient_with_scores.sh", 'w') as file_write :
	file_write.write("#!/bin/bash\n")
	for patient in patients_list : 
		file_write.write("python scripts/plot_kneepoint1_sim_vs_exp_with_scores.py %s kneepoint1_class1_simu_%s\n" % (patient,patient))



stop = timeit.default_timer()
print(stop - start)  




