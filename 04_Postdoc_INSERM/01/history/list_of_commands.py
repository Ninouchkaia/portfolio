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

for patient in patients_list : 
	with open("patient_command_%s.sh" % patient, 'w') as file_write :
		file_write.write("#!/bin/bash\n")
		file_write.write("echo \"This is a shell script for patient %s\"\n" % patient)
		file_write.write("python scripts/aggregateData.py %s/ABM_2D_%s %s/outputs_ABM_2D_%s.txt\n\n" % (patient,patient,patient,patient))
		file_write.write("python scripts/remove_duplicates_generic_with_filtering_keeping_only_samples.py %s/outputs_ABM_2D_%s.txt fitnessVia 50\n\n" % (patient, patient))
		file_write.write("python scripts/copy_for_git_paretoFrontGenericStochastic.py %s/outputs_ABM_2D_%s_duplicates_removed_filtered_only_samples_kept_50.0.txt %s/pareto_ABM_2D_%s\n\n" % (patient,patient,patient,patient))
		file_write.write("python scripts/extract_param_sets_from_pareto_adapted.py %s/pareto_ABM_2D_%s.txt %s/best_param_sets_ABM_2D_%s.tsv\n\n" % (patient,patient,patient,patient))
		file_write.write("python scripts/parse_best_param_for_behavior_space_adapted.py %s/best_param_sets_ABM_2D_%s.tsv %s/netlogo_best_param_sets_ABM_2D_%s.txt\n\n" % (patient,patient,patient,patient))
		file_write.write("python scripts/make_behavior_space_experiment_file.py %s/best_param_sets_ABM_2D_%s.tsv %s/experiment_file.xml\n\n" % (patient,patient,patient))
		file_write.write("/I/Program\ Files/NetLogo\ 6.1.0/netlogo-headless.bat --model %s/ABM_2D_%s.nlogo --setup-file %s/experiment_file.xml --experiment stocha_best_via --table %s/BehaviorSpace/stocha_best_via.csv --threads 4\n\n" % (patient,patient,patient,patient))
		file_write.write("/I/Program\ Files/NetLogo\ 6.1.0/netlogo-headless.bat --model %s/ABM_2D_%s.nlogo --setup-file %s/experiment_file.xml --experiment stocha_knee_point --table %s/BehaviorSpace/stocha_knee_point.csv --threads 4\n\n" % (patient,patient,patient,patient))
		file_write.write("/I/Program\ Files/NetLogo\ 6.1.0/netlogo-headless.bat --model %s/ABM_2D_%s.nlogo --setup-file %s/experiment_file.xml --experiment stocha_best_conc --table %s/BehaviorSpace/stocha_best_conc.csv --threads 4\n\n" % (patient,patient,patient,patient))
		file_write.write("python scripts/plot_sim_vs_exp.py %s stocha_best_via\n\n" % patient)
		file_write.write("python scripts/plot_sim_vs_exp.py %s stocha_knee_point\n\n" % patient)
		file_write.write("python scripts/plot_sim_vs_exp.py %s stocha_best_conc\n\n" % patient)



stop = timeit.default_timer()
print(stop - start)  



# os.system('ls -l')

# import subprocess


# os.system("I:\\Program Files\\NetLogo\\ 6.1.0\\netlogo-headless.bat --model CAS1802\\ABM_2D_CAS1802.nlogo --setup-file CAS1802\\experiment_file.xml --experiment stocha_best_conc --table CAS1802\\BehaviorSpace\\stocha_best_conc.csv --threads 4")
# subprocess.call("I:\\Program Files\\NetLogo\\ 6.1.0\\netlogo-headless.bat --model CAS1802\\ABM_2D_CAS1802.nlogo --setup-file CAS1802\\experiment_file.xml --experiment stocha_best_conc --table CAS1802\\BehaviorSpace\\stocha_best_conc.csv --threads 4")
# subprocess.run(["I:/Program Files/NetLogo 6.1.0/netlogo-headless.bat" , "--model CAS1802/ABM_2D_CAS1802.nlogo",  "--setup-file CAS1802/experiment_file.xml", "--experiment stocha_best_conc" , "--table CAS1802/BehaviorSpace/stocha_best_conc.csv",  "--threads 4"])

# subprocess.run(["ls", "-lha"],shell=True)