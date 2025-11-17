#!/usr/bin/python
# -*- coding: utf-8 -*-
# jupyter notebook --browser='C:/Program Files (x86)/Google/Chrome/Application/chrome.exe %s'
import timeit
start = timeit.default_timer()


import os
import pandas as pd

patient_dict = {}

patients_list_with_mono = ['CRE1704-1.1%','LAU1405-2.5%','DES2105-1.25%','ORE1706-0.68%','CAS1802-1.04%', 'GER160522-0.45%','REI230522-0.95%','PUJ240522-0.34%','LAR300522-0.21%','CAZ310522-3.48%']


patients_list = []

for patient in patients_list_with_mono :
	patient_name = patient.split("-")[0]
	patient_mono = float(patient.split("-")[1][:-1])
	patients_list.append(patient_name)
	patient_dict[patient_name] = [patient_mono]

print(patients_list)


### fetch viability and concentration data for each patient from the fused dataframes

## viability
df_viability = pd.read_table('fused_10patients_viability_renamed.tsv', index_col='Day')
print(df_viability)
print(df_viability.columns)

### return %apo-init cells from viability at Day 0
for patient in patients_list :
	apo_init = 100 - df_viability.at[0,patient]
	print(df_viability.at[0,patient], "{:.2f}".format(apo_init))
	patient_dict[patient].append(float("{:.2f}".format(apo_init)))

### and store it into the patient_dict	
for patient in patient_dict :
	print(patient, patient_dict[patient][0], patient_dict[patient][1])



## concentration
df_concentration = pd.read_table('fused_10patients_concentration_renamed.tsv', index_col='Day')
print(df_concentration)
print(df_concentration.columns)

# patients_list = ['CRE1704']

for patient in patients_list :
	os.mkdir('%s' % patient)
	print(patient)
	viability = df_viability.loc[:,patient]
	concentration = df_concentration.loc[:,patient]
	df_patient = pd.merge(viability, concentration, on='Day', suffixes=('_viability','_concentration'))
	df_patient=df_patient.dropna(axis=0)
	print(df_patient)
	df_patient.to_csv('%s\\%s.csv' % (patient,patient), index=True)
	days_measured = [str(x * 24) for x in list(df_patient.index.values)]
	time_series = "%s" % ' '.join(days_measured)
	print(time_series)

	with open('%s\\ABM_2D_%s.nlogo' % (patient,patient), 'w') as file_write :
		with open('ABM_2D_10patients_16.nlogo', 'r') as file_read :
			data = file_read.readlines()
			
			for line in data :
				if '  if (member? ticks [0 24 48 72 144 168 192 216 240 311])' in line :
					print(line)
					adapted_line = '  if (member? ticks [%s])' % time_series
					file_write.write('%s\n' % adapted_line)
				else :
					file_write.write(line)

	with open('%s\\ABM_2D_%s.oms' % (patient,patient), 'w') as file_write :
		with open("ABM_2D_10patients_18.oms", 'r') as file_read :
			data = file_read.readlines()
			for line in data :
				if '    NetLogo6Task(workDirectory / "ABM_2D_10patients_16.nlogo", launch, embedWorkspace = false, switch3d = false, seed = mySeed) set( ' in line :
					adapted_line = '    NetLogo6Task(workDirectory / \"ABM_2D_%s.nlogo\", launch, embedWorkspace = false, switch3d = false, seed = mySeed) set( ' % patient
					file_write.write('%s\n' % adapted_line)
				
				elif '    dataFile1 := (workDirectory / "filtered_fused_10patients.csv"),' in line :
					adapted_line = '    dataFile1 := (workDirectory / \"%s.csv\"),' % patient
					file_write.write('%s\n' % adapted_line)
				
				elif '        evaluation = modelTask(1.2, 4.32, simuDur) -- fitnessTask,' in line :
					adapted_line = '        evaluation = modelTask(%s, %s, simuDur) -- fitnessTask,' % (patient_dict[patient][0],patient_dict[patient][1])
					file_write.write('%s\n' % adapted_line)

				elif '    ) on ifb hook(workDirectory / "ABM_2D_10_patients_18", 1) // 1 save each 1 pop' in line :
					adapted_line = '    ) on ifb hook(workDirectory / \"ABM_2D_%s\", 1) // 1 save each 1 pop' % patient
					file_write.write('%s\n' % adapted_line)
				else :
					file_write.write(line)


stop = timeit.default_timer()
print(stop - start) 