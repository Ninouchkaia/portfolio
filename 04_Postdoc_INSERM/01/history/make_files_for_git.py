#!/usr/bin/python
# -*- coding: utf-8 -*-
# jupyter notebook --browser='C:/Program Files (x86)/Google/Chrome/Application/chrome.exe %s'
import timeit
start = timeit.default_timer()


import os
import pandas as pd

patient_dict = {}

patients_list_with_mono = ['CRE1704-1.1%','LAU1405-2.5%','DES2105-1.25%','ORE1706-0.68%','CAS1802-1.04%','REI230522-0.95%','PUJ240522-0.34%','LAR300522-0.21%','CAZ310522-3.48%']


patients_list = []

for patient in patients_list_with_mono :
	patient_name = patient.split("-")[0]
	patient_mono = float(patient.split("-")[1][:-1])
	patients_list.append(patient_name)
	patient_dict[patient_name] = [patient_mono]

print(patients_list)

os.mkdir("patient_data_for_git")
for patient in patients_list :
	patient_alias = patients_list.index(patient) + 1
	os.mkdir('patient_data_for_git\\patient_%s' % patient_alias)
	print(patient)
	print(patient_alias)
	
	with open('patient_data_for_git\\patient_%s\\pareto_front_patient_%s.txt' % (patient_alias,patient_alias), 'w') as file_write :
		with open('%s\\pareto_ABM_2D_%s.txt' % (patient,patient), 'r') as file_read :
			data = file_read.readlines()
			data[0] = "delta_fitness_via,delta_fitness_conc, apoCellsMovementProba, needSigCellsMvtProba, layersAroundNLC, antiApoBoost, monoPhagoEff, nlcPhagoEff, macroPhagoEff, macroKillEff, cllSensingDistance, monocyteSensingDistance, nlcSensingDistance, macrophageSensingDistance, nlcThreshold, signalInitMean, signalInitStd, monoDiffThreshold , monoDiffTimeStd, gammaLifeInitRate, gammaLifeInitShape\n"
		file_write.write(''.join(data))

	with open('patient_data_for_git\\patient_%s\\NSGAII_exploration_output_patient_%s.txt' % (patient_alias,patient_alias), 'w') as file_write :
		with open('%s\\outputs_ABM_2D_%s_duplicates_removed_filtered_only_samples_kept_50.0.txt' % (patient,patient), 'r') as file_read :
			data = file_read.readlines()
			data[0] = "apoCellsMovementProba,needSigCellsMvtProba,layersAroundNLC,antiApoBoost,monoPhagoEff,nlcPhagoEff,macroPhagoEff,macroKillEff,cllSensingDistance,monocyteSensingDistance,nlcSensingDistance,macrophageSensingDistance,nlcThreshold,signalInitMean,signalInitStd,monoDiffThreshold,monoDiffTimeStd,gammaLifeInitRate,gammaLifeInitShape,fitnessVia,fitnessConc,evolution$samples\n"
		file_write.write(''.join(data))
	

stop = timeit.default_timer()
print(stop - start) 