#!/usr/bin/python
# -*- coding: utf-8 -*-
# jupyter notebook --browser='C:/Program Files (x86)/Google/Chrome/Application/chrome.exe %s'


import timeit
start = timeit.default_timer()

import sys
import os
import pandas as pd

# my_file = sys.argv[1]
# output = sys.argv[2]

### usage : python scripts/make_best_sets_all_patients.py best_via_set 
### usage : python scripts/make_best_sets_all_patients.py best_conc_set 
### usage : python scripts/make_best_sets_all_patients.py knee_point_set 

type_of_set = sys.argv[1] ## 'best_via_set', 'best_conc_set', 'knee_point_set'

patients_list_with_mono = ['CRE1704-1.1%','LAU1405-2.5%','DES2105-1.25%','ORE1706-0.68%','CAS1802-1.04%', 'GER160522-0.45%','REI230522-0.95%','PUJ240522-0.34%','LAR300522-0.21%','CAZ310522-3.48%']
patients_list = []

for patient in patients_list_with_mono :
	patient_name = patient.split("-")[0]
	patients_list.append(patient_name)
print(patients_list)

best_set_all_patients = pd.DataFrame(columns = patients_list)
for patient in patients_list :
	df = pd.read_csv('%s\\best_param_sets_ABM_2D_%s.tsv' % (patient,patient), sep='\t', header=0, index_col=0)
	# print(df)
	best_set = df.loc[type_of_set]
	# print(best_via_set)
	best_set_all_patients[patient] = best_set.T

# print(best_via_set_all_patients)

best_set_all_patients.to_csv('%ss_all_patients.tsv' % type_of_set, index=True, sep='\t')


stop = timeit.default_timer()
print(stop - start)  


