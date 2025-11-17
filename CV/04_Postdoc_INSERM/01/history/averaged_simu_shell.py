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

with open("patient_command_averaged_class1_simu.sh", 'w') as file_write :
	file_write.write("#!/bin/bash\n")
	for patient in patients_list : 
		file_write.write("/I/Program\ Files/NetLogo\ 6.1.0/netlogo-headless.bat --model %s/ABM_2D_%s.nlogo --setup-file class1_averaged.xml --experiment averaged_simu_%s --table %s/BehaviorSpace/averaged_class1_simu_%s.csv --threads 4\n\n" % (patient,patient,patient,patient,patient))



stop = timeit.default_timer()
print(stop - start)  



