#!/usr/bin/python
# -*- coding: utf-8 -*-

import timeit
start = timeit.default_timer()
import os
import sys
import re
import math
from scipy.stats.stats import pearsonr
import matplotlib
import pylab


# genome_list = ['3702-A.thaliana']
genome_list = ["10090-M.musculus", "160490-S.pyogenes", "3702-A.thaliana", "449447-M.aeruginosa", "4896-S.pombe", "4932-S.cerevisiae", "511145-E.coli", "593117-Spectral_counting_T.gammatolerans", "6239-C.elegans", "7227-D.melanogaster", "83332-M.tuberculosis", "9031-Spectral_counting_G.gallus", "9606-H.sapiens", "9913-B.taurus", "99287-Spectral_counting_S.typhimurium"]
#took off B.subtili and L.interrogans

for genome in genome_list :
	print genome
	file_open = open("abundance_datasets/tables/tables_2bis/%s_table2bis.txt" % genome, 'rU')
	abundance_data = file_open.readlines()
	file_open.close()
	p_tot = 0
	for line in abundance_data[6:] :
		line = line.split("\t")
		string_id = line[0]
		aa_count = "\t".join(line[2:22])
		aa_sum = float(line[22])
		abundance = float(line[23])
		p_prime = aa_sum * abundance
		p_tot = p_tot + p_prime

	file_open = open("abundance_datasets/tables/tables_4/temp/%s_tabletemp1.txt" % genome, 'rU')
	temp_data = file_open.readlines()
	file_open.close()

	x = []
	y = []
	for line in temp_data[1:] :
		line = line.replace("\n", "")
		line = line.split("\t")
		p_prime = float(line[23])
		p = p_prime/p_tot
		abundance = line[22]
		q = [0]*20
		r = [0]*20
		r_ln = [0]*20
		H = 0
		for k in range(0,20) :
			q[k] = (float(line[k+1]) * float(abundance)/p_tot)
			r[k] = q[k]/p
			if r[k] == 0.0 :
				r_ln[k] = 0
			else :
				r_ln[k] = r[k]*(math.log(r[k]))
			H = H + r_ln[k]
		r_e = [0]*20
		r_e[0] = r[0] * 12
		r_e[1] = r[1] * 741
		r_e[2] = r[2] * 127
		r_e[3] = r[3] * 77
		r_e[4] = r[4] * 156
		r_e[5] = r[5] * 12
		r_e[6] = r[6] * 536
		r_e[7] = r[7] * 65
		r_e[8] = r[8] * 242
		r_e[9] = r[9] * 55
		r_e[10] = r[10] * 446
		r_e[11] = r[11] * 162
		r_e[12] = r[12] * 41
		r_e[13] = r[13] * 147
		r_e[14] = r[14] * 109
		r_e[15] = r[15] * 70
		r_e[16] = r[16] * 112
		r_e[17] = r[17] * 47
		r_e[18] = r[18] * 892
		r_e[19] = r[19] * 350

		diff = -(H) - sum(r_e) 

		x.append(diff)
		y.append(math.log(p))

	print len(x), len(y)
	print pearsonr(x, y)
	matplotlib.pyplot.scatter(x,y)

	matplotlib.pyplot.show()		










stop = timeit.default_timer()
print stop - start 