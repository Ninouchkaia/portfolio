#!/usr/bin/python
# -*- coding: utf-8 -*-
# jupyter notebook --browser='C:/Program Files (x86)/Google/Chrome/Application/chrome.exe %s'


import timeit
start = timeit.default_timer()

import sys
import os

my_file = sys.argv[1]
output = sys.argv[2]

patient = my_file.split("/")[0]
with open("patient_dict.txt", 'r') as file_read :
	data = file_read.readlines()
	for line in data[1:] :
		line = line.replace("\n", "").split(" ")
		patient_name = line[0]
		if patient_name == patient :
			mono_init = line[1]
			apo_init = line[2]

os.mkdir('%s/BehaviorSpace' % patient)

with open(my_file, 'r') as file_read :
	with open("%s" % output, 'w') as file_write :
		data = file_read.readlines()
		n = 1
		for line in data[1:] :
			line = line.replace("\n", "").split("\t")
			gui_apo_mov = line[3]
			gui_need_sig_mov = line[4]
			gui_layers = line[5]
			gui_alpha = line[6]
			gui_mono_phago_eff = line[7]
			gui_NLC_phago_eff = line[8]
			gui_M_phago_eff = line[9]
			gui_M_kill_eff = line[10]
			gui_cll_sens_dist = line[11]
			gui_mono_sens_dist = line[12]
			gui_nlc_sens_dist = line[13]
			gui_macro_sens_dist = line[14]
			gui_nlc_threshold = line[15]
			gui_sig_init = line[16]
			gui_sig_init_std = line[17]
			gui_diff_mean = line[18]
			gui_diff_std = line[19]
			gui_life_init_gamma = line[20]
			gui_alpha_distrib = line[21]
			if (n == 1) :
				file_write.write("Best_via_set\n")
			elif (n == 2) :
				file_write.write("Knee_point_set\n")
			elif (n == 3) :
				file_write.write("Best_conc_set\n")
			file_write.write("[\"gui-prop-mono-init\" %s]\n" % mono_init)
			file_write.write("[\"gui-prop-apo-init\" %s]\n" % apo_init)				
			file_write.write("[\"gui-apo-mov\" %s]\n" % gui_apo_mov)
			file_write.write("[\"gui-need-sig-mov\" %s]\n" % gui_need_sig_mov)
			file_write.write("[\"gui-layers\" %s]\n" % gui_layers)
			file_write.write("[\"gui-alpha\" %s]\n" % gui_alpha)
			file_write.write("[\"gui-mono-phago-eff\" %s]\n" % gui_mono_phago_eff)
			file_write.write("[\"gui-NLC-phago-eff\" %s]\n" % gui_NLC_phago_eff)
			file_write.write("[\"gui-M-phago-eff\" %s]\n" % gui_M_phago_eff)
			file_write.write("[\"gui-M-kill-eff\" %s]\n" % gui_M_kill_eff)
			file_write.write("[\"gui-cll-sens-dist\" %s]\n" % gui_cll_sens_dist)
			file_write.write("[\"gui-mono-sens-dist\" %s]\n" % gui_mono_sens_dist)
			file_write.write("[\"gui-nlc-sens-dist\" %s]\n" % gui_nlc_sens_dist)
			file_write.write("[\"gui-macro-sens-dist\" %s]\n" % gui_macro_sens_dist)
			file_write.write("[\"gui-nlc-threshold\" %s]\n" % gui_nlc_threshold)
			file_write.write("[\"gui-sig-init\" %s]\n" % gui_sig_init)
			file_write.write("[\"gui-sig-init-std\" %s]\n" % gui_sig_init_std)
			file_write.write("[\"gui-diff-mean\" %s]\n" % gui_diff_mean)
			file_write.write("[\"gui-diff-std\" %s]\n" % gui_diff_std)
			file_write.write("[\"gui-life-init-gamma\" %s]\n" % gui_life_init_gamma)
			file_write.write("[\"gui-alpha-distrib\" %s]\n\n" % gui_alpha_distrib)

			n = n+1

print(n)
stop = timeit.default_timer()
print(stop - start)  


