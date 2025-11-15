#!/usr/bin/python
# -*- coding: utf-8 -*-
# jupyter notebook --browser='C:/Program Files (x86)/Google/Chrome/Application/chrome.exe %s'


import timeit
start = timeit.default_timer()

import pandas as pd

patients_list_with_mono = ['CRE1704-1.1%','LAU1405-2.5%','DES2105-1.25%','ORE1706-0.68%','CAS1802-1.04%', 'GER160522-0.45%','REI230522-0.95%','PUJ240522-0.34%','LAR300522-0.21%','CAZ310522-3.48%']

patients_list = []
for patient in patients_list_with_mono :
	patient_name = patient.split("-")[0]
	patients_list.append(patient_name)
print(patients_list)

df_summary_knee_points = pd.DataFrame(columns = patients_list)

for patient in patients_list :
	print(patient)
	my_file = "%s\\best_param_sets_ABM_2D_%s.tsv" % (patient,patient)
	df = pd.read_table(my_file, index_col=0)
	df_knee_point = df.loc['knee_point_set']
	print(df_knee_point)
	df_summary_knee_points[patient] = df_knee_point

print(df_summary_knee_points)

df_summary_knee_points.to_csv("summary_knee_points.csv", index=True, sep=',')

stop = timeit.default_timer()
print(stop - start)  


