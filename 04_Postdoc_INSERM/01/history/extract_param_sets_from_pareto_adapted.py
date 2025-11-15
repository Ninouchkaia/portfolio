#!/usr/bin/python
# -*- coding: utf-8 -*-
# jupyter notebook --browser='C:/Program Files (x86)/Google/Chrome/Application/chrome.exe %s'


import timeit
start = timeit.default_timer()

import sys
import numpy as np
import math
import pandas as pd

my_file = sys.argv[1]
output = sys.argv[2]

pareto_front = pd.read_csv("%s" % my_file)
pareto_front['distances'] = np.sqrt(pareto_front['delta_fitness_via'] * pareto_front['delta_fitness_via'] + pareto_front['delta_fitness_conc'] * pareto_front['delta_fitness_conc'])

best_via_set = pareto_front[pareto_front.delta_fitness_via == pareto_front.delta_fitness_via.min()]
best_conc_set = pareto_front[pareto_front.delta_fitness_conc == pareto_front.delta_fitness_conc.min()]
knee_point_set = pareto_front[pareto_front.distances == pareto_front.distances.min()]

best_params = best_via_set
best_params = best_params.append(knee_point_set)
best_params = best_params.append(best_conc_set)
best_params = best_params.drop('distances', 1)

best_params.insert(loc=0, value=['best_via_set', 'knee_point_set', 'best_conc_set'], column="set")

best_params.to_csv('%s' % output, sep='\t', index=False)  



stop = timeit.default_timer()
print(stop - start)  