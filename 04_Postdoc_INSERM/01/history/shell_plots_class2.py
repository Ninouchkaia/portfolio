#!/usr/bin/python
# -*- coding: utf-8 -*-
# jupyter notebook --browser='C:/Program Files (x86)/Google/Chrome/Application/chrome.exe %s'


import timeit
start = timeit.default_timer()

exp_list = ["perturb-gui-apo-mov","perturb-gui-need-sig-mov","perturb-gui-layers","perturb-gui-alpha","perturb-gui-mono-phago-eff","perturb-gui-NLC-phago-eff","perturb-gui-M-phago-eff","perturb-gui-M-kill-eff","perturb-gui-cll-sens-dist","perturb-gui-mono-sens-dist","perturb-gui-nlc-sens-dist","perturb-gui-macro-sens-dist","perturb-gui-nlc-threshold","perturb-gui-sig-init","perturb-gui-sig-init-std","perturb-gui-diff-mean","perturb-gui-diff-std","perturb-gui-life-init-gamma","perturb-gui-alpha-distrib"]


with open("plot_sensitivity_class_2_experiments.sh" , 'a+') as file_write :
	file_write.write("#!/bin/bash\n")
	for exp in exp_list : 

		file_write.write("echo \"This is a shell script for exp %s\"\n" % exp)
		file_write.write(" python sensitivity_analysis/plot_sensitivity_analysis_class2.py %s\n" % exp)



stop = timeit.default_timer()
print(stop - start)