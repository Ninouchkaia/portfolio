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

#tabletemp1 and tabletemp2 have to be first generated with script called newdef_protein.py

# this script will plot the relationship between ln(q_k) and e_k 

# genome_list = ['3702-A.thaliana']
genome_list = ["10090-M.musculus", "160490-S.pyogenes", "3702-A.thaliana", "449447-M.aeruginosa", "4896-S.pombe", "4932-S.cerevisiae", "511145-E.coli", "593117-Spectral_counting_T.gammatolerans", "6239-C.elegans", "7227-D.melanogaster", "83332-M.tuberculosis", "9031-Spectral_counting_G.gallus", "9606-H.sapiens", "9913-B.taurus", "99287-Spectral_counting_S.typhimurium"]
#took off B.subtili and L.interrogans

for genome in genome_list :
	print genome
	file_open = open("abundance_datasets/tables/tables_2bis/%s_table2bis.txt" % genome, 'rU')
	abundance_data = file_open.readlines()
	file_open.close()
	
	file_write = open("abundance_datasets/tables/tables_4/temp/%s_tabletemp0.txt" % genome, 'a+')
	file_write.write("\t%s\n" % ( "\t".join( ((abundance_data[0]).split("\t"))[2:] ) )  )
	file_write.write("%s" % ("\t".join( ((abundance_data[4]).split("\t"))[1:]))) 
	e_A = 12
	e_C = 741
	e_D = 127
	e_E = 77
	e_F = 156
	e_G = 12
	e_H = 536
	e_I = 65
	e_K = 242
	e_L = 55
	e_M = 446
	e_N = 162
	e_P = 41
	e_Q = 147
	e_R = 109
	e_S = 70
	e_T = 112
	e_V = 47
	e_W = 892
	e_Y = 350

	y = [e_A,e_C,e_D,e_E,e_F,e_G,e_H,e_I,e_K,e_L,e_M,e_N,e_P,e_Q,e_R,e_S,e_T,e_V,e_W,e_Y]
	
	ln_q_A = math.log(float((abundance_data[4].split("\t"))[2]))
	ln_q_C = math.log(float((abundance_data[4].split("\t"))[3]))
	ln_q_D = math.log(float((abundance_data[4].split("\t"))[4]))
	ln_q_E = math.log(float((abundance_data[4].split("\t"))[5]))
	ln_q_F = math.log(float((abundance_data[4].split("\t"))[6]))
	ln_q_G = math.log(float((abundance_data[4].split("\t"))[7]))
	ln_q_H = math.log(float((abundance_data[4].split("\t"))[8]))
	ln_q_I = math.log(float((abundance_data[4].split("\t"))[9]))
	ln_q_K = math.log(float((abundance_data[4].split("\t"))[10]))
	ln_q_L = math.log(float((abundance_data[4].split("\t"))[11]))
	ln_q_M = math.log(float((abundance_data[4].split("\t"))[12]))
	ln_q_N = math.log(float((abundance_data[4].split("\t"))[13]))
	ln_q_P = math.log(float((abundance_data[4].split("\t"))[14]))
	ln_q_Q = math.log(float((abundance_data[4].split("\t"))[15]))
	ln_q_R = math.log(float((abundance_data[4].split("\t"))[16]))
	ln_q_S = math.log(float((abundance_data[4].split("\t"))[17]))
	ln_q_T = math.log(float((abundance_data[4].split("\t"))[18]))
	ln_q_V = math.log(float((abundance_data[4].split("\t"))[19]))
	ln_q_W = math.log(float((abundance_data[4].split("\t"))[20]))
	ln_q_Y = math.log(float((abundance_data[4].split("\t"))[21]))

	x = ln_q_A,ln_q_C,ln_q_D,ln_q_E,ln_q_F,ln_q_G,ln_q_H,ln_q_I,ln_q_K,ln_q_L,ln_q_M,ln_q_N,ln_q_P,ln_q_Q,ln_q_R,ln_q_S,ln_q_T,ln_q_V,ln_q_W,ln_q_Y

	file_write.write("ln(q_k)\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (ln_q_A,ln_q_C,ln_q_D,ln_q_E,ln_q_F,ln_q_G,ln_q_H,ln_q_I,ln_q_K,ln_q_L,ln_q_M,ln_q_N,ln_q_P,ln_q_Q,ln_q_R,ln_q_S,ln_q_T,ln_q_V,ln_q_W,ln_q_Y))
	file_write.write("energies\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t" % (e_A,e_C,e_D,e_E,e_F,e_G,e_H,e_I,e_K,e_L,e_M,e_N,e_P,e_Q,e_R,e_S,e_T,e_V,e_W,e_Y))

	print pearsonr(x, y)
	matplotlib.pyplot.scatter(x,y)

	matplotlib.pyplot.show()

	






stop = timeit.default_timer()
print stop - start 