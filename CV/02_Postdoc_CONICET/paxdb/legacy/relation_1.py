#!/usr/bin/python
# -*- coding: utf-8 -*-

import timeit
start = timeit.default_timer()
from math import fsum, log
import mpmath

d_A = 12
d_C = 741
d_D = 114
d_d = 77
d_F = 208
d_G = 12
d_H = 536
d_I = 65
d_K = 242
d_L = 55
d_M = 446
d_N = 147
d_P = 61
d_Q = 130
d_R = 109
d_S = 70
d_T = 112
d_V = 47
d_W = 892
d_Y = 350

aa_cost_list=[d_A,d_C,d_D,d_d,d_F,d_G,d_H,d_I,d_K,d_L,d_M,d_N,d_P,d_Q,d_R,d_S,d_T,d_V,d_W,d_Y]

decay_dict = {}
with open("yeast_table_with_decay.txt", 'rU') as file_open :
	decay_data = file_open.readlines()
	for line in decay_data[1:] :
		line = line.split("\t")
		string_id = line[0]
		decay = float(line[6].replace("\n", ""))
		if decay < 0 :
			decay = 300
		decay_dict[string_id] = decay

SUM_total = 0
with open("4932-S.cerevisiae_tabletemp2_update.txt", 'rU') as file_open :
	yeast_data = file_open.readlines()
	front_line = yeast_data[0].split("\t")

	with open("yeast_relation_1_A.txt", 'a+') as file_write :
		for i in range(0,23):
			file_write.write("%s\t" % front_line[i])
		for i in range(26,46):
			file_write.write("ln(%s)\t" % front_line[i])
		file_write.write("Protein_cost E\n")

		for line in yeast_data[1:] :
			line = line.split("\t")
			string_id = line[0]

			if string_id in decay_dict :
				protein_cost = decay_dict[string_id]
				exp_SUM_list = []
				for aa_cost in aa_cost_list :
					exp_SUM_list.append(mpmath.exp(-protein_cost-aa_cost))
				for i in range(0,23):
					file_write.write("%s\t" % line[i])
				for i in range(26,46):
					if float(line[i]) > 0 :
						ln_q = log(float(line[i]))
					else :
						ln_q = 0
					file_write.write("%s\t" % ln_q)
				file_write.write("%s\n" % protein_cost)
				SUM = fsum(exp_SUM_list)
				SUM_total = SUM_total + SUM
				intercept = log(SUM_total)

	with open("yeast_relation_1_B.txt", 'a+') as file_write :
		for i in range(0,23):
			file_write.write("%s\t" % front_line[i])
		for i in range(26,46):
			file_write.write("ln(%s)\t" % front_line[i])
		file_write.write("Protein_cost E\t")
		# file_write.write("intercept = ln(Sum_j(exp-(E+d_j))\t")
		for aa_cost in aa_cost_list[:-1] :
			file_write.write("E+%s\t" % aa_cost)
		file_write.write("E+%s\n" % aa_cost_list[-1])

		for line in yeast_data[1:] :
			line = line.split("\t")
			string_id = line[0]
			if string_id in decay_dict :
				protein_cost = decay_dict[string_id]
				exp_SUM_list = []
				for aa_cost in aa_cost_list :
					exp_SUM_list.append(mpmath.exp(-protein_cost-aa_cost))
				for i in range(0,23):
					file_write.write("%s\t" % line[i])
				for i in range(26,46):
					if float(line[i]) > 0 :
						ln_q = log(float(line[i]))
					else :
						ln_q = 0
					file_write.write("%s\t" % ln_q)
				file_write.write("%s\t" % protein_cost)
				# file_write.write("%s\t" % intercept)
				for aa_cost in aa_cost_list[:-1] :
					value = protein_cost + aa_cost
					file_write.write("%s\t" % value)
				value = protein_cost + aa_cost_list[-1]
				file_write.write("%s\n" % value)









				


stop = timeit.default_timer()
print stop - start 