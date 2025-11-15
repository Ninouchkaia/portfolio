#!/usr/bin/python
# -*- coding: utf-8 -*-
# jupyter notebook --browser='C:/Program Files (x86)/Google/Chrome/Application/chrome.exe %s'


import timeit
start = timeit.default_timer()


### usage $ python sensitivity_analysis/make_behavior_space_experiment_file.py best_param_sets_pareto_ABM_2D_9patients_0_50.tsv sensitivity_analysis_experiment_file.xml


import sys

my_file = sys.argv[1]
output = sys.argv[2]

gui_prop_mono_init = 1.28
gui_prop_apo_init = 4.55

ranges_dict = {}
### build the dict with param names and the value ranges to perturb
ranges_dict["gui-apo-mov"] = [0,2,10]
ranges_dict["gui-need-sig-mov"] = [0,2,10]
ranges_dict["gui-layers"] = [1,1,3]
ranges_dict["gui-alpha"] = [0,50,300]
ranges_dict["gui-mono-phago-eff"] = [0,10,100]
ranges_dict["gui-NLC-phago-eff"] = [0,10,100]
ranges_dict["gui-M-phago-eff"] = [0,10,100]
ranges_dict["gui-M-kill-eff"] = [0,1,5]
ranges_dict["gui-cll-sens-dist"] = [1,1,3]
ranges_dict["gui-mono-sens-dist"] = [1,1,3]
ranges_dict["gui-nlc-sens-dist"] = [1,1,3]
ranges_dict["gui-macro-sens-dist"] = [1,1,3]
ranges_dict["gui-nlc-threshold"] = [90,20,210]
ranges_dict["gui-sig-init"] = [0,12,72]
ranges_dict["gui-sig-init-std"] = [0,8,48]
ranges_dict["gui-diff-mean"] = [48,4,72]
ranges_dict["gui-diff-std"] = [0,8,48]
ranges_dict["gui-life-init-gamma"] = [50,125,2500]
ranges_dict["gui-alpha-distrib"] = [0.1,0.05,1.0]

knee_point_dict = {}
with open(my_file, 'r') as file_read :
  data = file_read.readlines()
  for line in data[1:] :
    # print(line)
    line = line.replace("\n", "").split("\t")
    # print(line[0])
    if line[0] == "knee_point_set" :
      knee_point_dict["gui-apo-mov"] = line[3]
      knee_point_dict["gui-need-sig-mov"] = line[4]
      knee_point_dict["gui-layers"] = line[5]
      knee_point_dict["gui-alpha"] = line[6]
      knee_point_dict["gui-mono-phago-eff"] = line[7]
      knee_point_dict["gui-NLC-phago-eff"] = line[8]
      knee_point_dict["gui-M-phago-eff"] = line[9]
      knee_point_dict["gui-M-kill-eff"] = line[10]
      knee_point_dict["gui-cll-sens-dist"] = line[11]
      knee_point_dict["gui-mono-sens-dist"] = line[12]
      knee_point_dict["gui-nlc-sens-dist"] = line[13]
      knee_point_dict["gui-macro-sens-dist"] = line[14]
      knee_point_dict["gui-nlc-threshold"] = line[15]
      knee_point_dict["gui-sig-init"] = line[16]
      knee_point_dict["gui-sig-init-std"] = line[17]
      knee_point_dict["gui-diff-mean"] = line[18]
      knee_point_dict["gui-diff-std"] = line[19]
      knee_point_dict["gui-life-init-gamma"] = line[20]
      knee_point_dict["gui-alpha-distrib"] = line[21]

# print(knee_point_dict)

with open("%s" % output, 'w') as file_write :
    file_write.write("<experiments>\n")      
    for param_name in ranges_dict :
      # print (param_name)
      ## put all the other param in a list (remaining_params) and iterate over it
      remaining_params = list(ranges_dict.keys())
      remaining_params.remove(param_name)
      # print(remaining_params)
      print("perturb-%s" % param_name)
      file_write.write("  <experiment name=\"perturb-%s\" repetitions=\"3\" runMetricsEveryStep=\"true\">\n" % param_name)
      file_write.write("    <setup>setup</setup>\n")
      file_write.write("    <go>go</go>\n")
      file_write.write("    <timeLimit steps=\"312\"/>\n")
      file_write.write("    <metric>getSeed</metric>\n")
      file_write.write("    <metric>getViability</metric>\n")
      file_write.write("    <metric>getRemainingCellRatio</metric>\n")

      file_write.write("    <enumeratedValueSet variable=\"gui-prop-mono-init\"><value value=\"%s\"/></enumeratedValueSet>\n" % gui_prop_mono_init)
      file_write.write("    <enumeratedValueSet variable=\"gui-prop-apo-init\"><value value=\"%s\"/></enumeratedValueSet>\n" % gui_prop_apo_init)  

      for other_param_name in remaining_params :
        file_write.write("    <enumeratedValueSet variable=\"%s\"><value value=\"%s\"/></enumeratedValueSet>\n" % (other_param_name, knee_point_dict[other_param_name]))

      file_write.write("    <steppedValueSet variable=\"%s\" first=\"%s\" step=\"%s\" last=\"%s\"/>\n" % (param_name, ranges_dict[param_name][0], ranges_dict[param_name][1], ranges_dict[param_name][2]))
      file_write.write("  </experiment>\n")

    file_write.write("  </experiments>\n")



stop = timeit.default_timer()
print(stop - start)