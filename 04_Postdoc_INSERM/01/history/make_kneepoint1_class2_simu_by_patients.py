#!/usr/bin/python
# -*- coding: utf-8 -*-
# jupyter notebook --browser='C:/Program Files (x86)/Google/Chrome/Application/chrome.exe %s'


import timeit
start = timeit.default_timer()


#### $ python scripts/make_kneepoint1_class2_simu_by_patients.py kneepointset1_class2.tsv kneepoint1_class2_simu_by_patient.xml

import sys

my_file = sys.argv[1]
output = sys.argv[2]

patient_dict = {}
with open("patient_dict.txt", 'r') as file_read :
    data = file_read.readlines()
    for line in data[1:] :
      line = line.replace("\n", "").split(" ")
      patient_name = line[0]
      gui_prop_mono_init = line[1]
      gui_prop_apo_init = line[2]
      patient_dict[patient_name] = [gui_prop_mono_init,gui_prop_apo_init]

with open(my_file, 'r') as file_read :
  with open("%s" % output, 'w') as file_write :
    file_write.write("<experiments>\n")
    data = file_read.readlines()
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
      for patient_name in patient_dict :
        file_write.write("  <experiment name=\"kneepoint1_class2_simu_%s\" repetitions=\"12\" runMetricsEveryStep=\"true\">\n" % patient_name)
        file_write.write("    <setup>setup</setup>\n")
        file_write.write("    <go>go</go>\n")
        file_write.write("    <timeLimit steps=\"312\"/>\n")
        file_write.write("    <metric>getSeed</metric>\n")
        file_write.write("    <metric>getViability</metric>\n")
        file_write.write("    <metric>getRemainingCellRatio</metric>\n")

        file_write.write("    <enumeratedValueSet variable=\"gui-prop-mono-init\"><value value=\"%s\"/></enumeratedValueSet>\n" % patient_dict[patient_name][0])
        file_write.write("    <enumeratedValueSet variable=\"gui-prop-apo-init\"><value value=\"%s\"/></enumeratedValueSet>\n" % patient_dict[patient_name][1])       
        file_write.write("    <enumeratedValueSet variable=\"gui-apo-mov\"><value value=\"%s\"/></enumeratedValueSet>\n" % gui_apo_mov)
        file_write.write("    <enumeratedValueSet variable=\"gui-need-sig-mov\"><value value=\"%s\"/></enumeratedValueSet>\n" % gui_need_sig_mov)
        file_write.write("    <enumeratedValueSet variable=\"gui-layers\"><value value=\"%s\"/></enumeratedValueSet>\n" % gui_layers)
        file_write.write("    <enumeratedValueSet variable=\"gui-alpha\"><value value=\"%s\"/></enumeratedValueSet>\n" % gui_alpha)
        file_write.write("    <enumeratedValueSet variable=\"gui-mono-phago-eff\"><value value=\"%s\"/></enumeratedValueSet>\n" % gui_mono_phago_eff)
        file_write.write("    <enumeratedValueSet variable=\"gui-NLC-phago-eff\"><value value=\"%s\"/></enumeratedValueSet>\n" % gui_NLC_phago_eff)
        file_write.write("    <enumeratedValueSet variable=\"gui-M-phago-eff\"><value value=\"%s\"/></enumeratedValueSet>\n" % gui_M_phago_eff)
        file_write.write("    <enumeratedValueSet variable=\"gui-M-kill-eff\"><value value=\"%s\"/></enumeratedValueSet>\n" % gui_M_kill_eff)
        file_write.write("    <enumeratedValueSet variable=\"gui-cll-sens-dist\"><value value=\"%s\"/></enumeratedValueSet>\n" % gui_cll_sens_dist)
        file_write.write("    <enumeratedValueSet variable=\"gui-mono-sens-dist\"><value value=\"%s\"/></enumeratedValueSet>\n" % gui_mono_sens_dist)
        file_write.write("    <enumeratedValueSet variable=\"gui-nlc-sens-dist\"><value value=\"%s\"/></enumeratedValueSet>\n" % gui_nlc_sens_dist)
        file_write.write("    <enumeratedValueSet variable=\"gui-macro-sens-dist\"><value value=\"%s\"/></enumeratedValueSet>\n" % gui_macro_sens_dist)
        file_write.write("    <enumeratedValueSet variable=\"gui-nlc-threshold\"><value value=\"%s\"/></enumeratedValueSet>\n" % gui_nlc_threshold)
        file_write.write("    <enumeratedValueSet variable=\"gui-sig-init\"><value value=\"%s\"/></enumeratedValueSet>\n" % gui_sig_init)
        file_write.write("    <enumeratedValueSet variable=\"gui-sig-init-std\"><value value=\"%s\"/></enumeratedValueSet>\n" % gui_sig_init_std)
        file_write.write("    <enumeratedValueSet variable=\"gui-diff-mean\"><value value=\"%s\"/></enumeratedValueSet>\n" % gui_diff_mean)
        file_write.write("    <enumeratedValueSet variable=\"gui-diff-std\"><value value=\"%s\"/></enumeratedValueSet>\n" % gui_diff_std)
        file_write.write("    <enumeratedValueSet variable=\"gui-life-init-gamma\"><value value=\"%s\"/></enumeratedValueSet>\n" % gui_life_init_gamma)
        file_write.write("    <enumeratedValueSet variable=\"gui-alpha-distrib\"><value value=\"%s\"/></enumeratedValueSet>\n" % gui_alpha_distrib)
        file_write.write("  </experiment>\n")
    
    file_write.write("  </experiments>\n")



stop = timeit.default_timer()
print(stop - start)  
