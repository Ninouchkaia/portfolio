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

genome_list = ["4932-S.cerevisiae"]
# genome_list = ['3702-A.thaliana']
# genome_list = ["10090-M.musculus", "160490-S.pyogenes", "3702-A.thaliana", "449447-M.aeruginosa", "4896-S.pombe", "4932-S.cerevisiae", "511145-E.coli", "593117-Spectral_counting_T.gammatolerans", "6239-C.elegans", "7227-D.melanogaster", "83332-M.tuberculosis", "9031-Spectral_counting_G.gallus", "9606-H.sapiens", "9913-B.taurus", "99287-Spectral_counting_S.typhimurium"]
#took off B.subtili and L.interrogans

for genome in genome_list :
	print genome
	# file_open = open("abundance_datasets/tables/tables_2bis/%s_table2bis.txt" % genome, 'rU')
	# abundance_data = file_open.readlines()
	# file_open.close()
	# file_write = open("abundance_datasets/tables/tables_update/%s_tabletemp1_update.txt" % genome, 'a+')
	# file_write.write("string_id\t%s\tall amino-acids\tabundance\tp_prime\n" % ( "\t".join( ((abundance_data[0]).split("\t"))[2:-1] ) )  )
	# p_tot = 0
	# for line in abundance_data[6:] :
	# 	line = line.split("\t")
	# 	string_id = line[0]
	# 	aa_count = "\t".join(line[2:22])
	# 	aa_sum = float(line[22])
	# 	abundance = float(line[23])
	# 	p_prime = aa_sum * abundance
	# 	p_tot = p_tot + p_prime
	# 	file_write.write("%s\t%s\t%s\t%s\t%s\n" % (string_id, aa_count, aa_sum, abundance, p_prime))
	# file_write.close()


	# file_open = open("abundance_datasets/tables/tables_update/%s_tabletemp1_update.txt" % genome, 'rU')
	# temp_data = file_open.readlines()
	# file_open.close()
	# file_write = open("abundance_datasets/tables/tables_update/%s_tabletemp2_update.txt" % genome, 'a+')
	# file_write.write(temp_data[0].replace("\n", ""))
	# file_write.write("\tp (=p_prime/p_tot)\t\tq_A\tq_C\tq_D\tq_E\tq_F\tq_G\tq_H\tq_I\tq_K\tq_L\tq_M\tq_N\tq_P\tq_Q\tq_R\tq_S\tq_T\tq_V\tq_W\tq_Y\t")
	# file_write.write("\tr_A\tr_C\tr_D\tr_E\tr_F\tr_G\tr_H\tr_I\tr_K\tr_L\tr_M\tr_N\tr_P\tr_Q\tr_R\tr_S\tr_T\tr_V\tr_W\tr_Y\t")
	# file_write.write("\tr_A*ln(r_A)\tr_C*ln(r_C)\tr_D*ln(r_D)\tr_E*ln(r_E)\tr_F*ln(r_F)\tr_G*ln(r_G)\tr_H*ln(r_H)\tr_I*ln(r_I)\tr_K*ln(r_K)\tr_L*ln(r_L)\tr_M*ln(r_M)\tr_N*ln(r_N)\tr_P*ln(r_P)\tr_Q*ln(r_Q)\tr_R*ln(r_R)\tr_S*ln(r_S)\tr_T*ln(r_T)\tr_V*ln(r_V)\tr_W*ln(r_W)\tr_Y*ln(r_Y)\t")
	# file_write.write("\tr_A*e_A\tr_C*e_C\tr_D*e_D\tr_E*e_E\tr_F*e_F\tr_G*e_G\tr_H*e_H\tr_I*e_I\tr_K*e_K\tr_L*e_L\tr_M*e_M\tr_N*e_N\tr_P*e_P\tr_Q*e_Q\tr_R*e_R\tr_S*e_S\tr_T*e_T\tr_V*e_V\tr_W*e_W\tr_Y*e_Y\t")
	# file_write.write("\tSUM(r_k*ln(r_k)) --> -H(amino-acids)\t")
	# file_write.write("\tSUM(r_k*e_k)\t")
	# file_write.write("\tdiff (H(amino-acids)-SUM(r_k*e_k))\t")
	# file_write.write("ln(p)\n")
	# for line in temp_data[1:] :
	# 	line = line.replace("\n", "")
	# 	file_write.write(line)
	# 	line = line.split("\t")
	# 	p_prime = float(line[23])
	# 	p = p_prime/p_tot
	# 	file_write.write("\t%s\t" % p)
	# 	abundance = line[22]
	# 	q = [0]*20
	# 	r = [0]*20
	# 	r_ln = [0]*20
	# 	H = 0
	# 	for k in range(0,20) :
	# 		q[k] = (float(line[k+1]) * float(abundance)/p_tot)
	# 		r[k] = q[k]/p
	# 		if r[k] == 0.0 :
	# 			r_ln[k] = 0
	# 		else :
	# 			r_ln[k] = r[k]*(math.log(r[k]))
	# 		H = H + r_ln[k]
	# 	r_e = [0]*20
	# 	r_e[0] = r[0] * 12
	# 	r_e[1] = r[1] * 741
	# 	r_e[2] = r[2] * 127
	# 	r_e[3] = r[3] * 77
	# 	r_e[4] = r[4] * 156
	# 	r_e[5] = r[5] * 12
	# 	r_e[6] = r[6] * 536
	# 	r_e[7] = r[7] * 65
	# 	r_e[8] = r[8] * 242
	# 	r_e[9] = r[9] * 55
	# 	r_e[10] = r[10] * 446
	# 	r_e[11] = r[11] * 162
	# 	r_e[12] = r[12] * 41
	# 	r_e[13] = r[13] * 147
	# 	r_e[14] = r[14] * 109
	# 	r_e[15] = r[15] * 70
	# 	r_e[16] = r[16] * 112
	# 	r_e[17] = r[17] * 47
	# 	r_e[18] = r[18] * 892
	# 	r_e[19] = r[19] * 350

	# 	SUM_rk_ek = sum(r_e)

	# 	diff = -(H) - SUM_rk_ek


				
		
	# 	for k in range(0,20) :
	# 		q[k] = str(q[k])
	# 		r[k] = str(r[k])
	# 		r_ln[k] = str(r_ln[k])
	# 		r_e[k] = str(r_e[k])

		
	# 	q = "\t".join(q)
	# 	file_write.write("\t%s\t" % q)

	# 	r = "\t".join(r)
	# 	file_write.write("\t%s\t" % r)

	# 	r_ln = "\t".join(r_ln)
	# 	file_write.write("\t%s\t" % r_ln)

	# 	r_e = "\t".join(r_e)
	# 	file_write.write("\t%s\t" % r_e)

	# 	file_write.write("\t%s\t" % H)

	# 	file_write.write("\t%s\t" % SUM_rk_ek)

	# 	file_write.write("\t%s\t" % diff)

	# 	file_write.write("%s\n" % math.log(p))

		


	# file_write.close()

	file_open = open("/Users/verstrat/temp_qb/paxdb/abundance_datasets/tables/tables_update/%s_tabletemp2_update.txt" % genome, 'rU')
	temp_data = file_open.readlines()
	file_open.close()
	file_write = open("/Users/verstrat/temp_qb/paxdb/yeast/%s_table.txt" % genome, 'a+')
	for line in temp_data : 
		line = line.split("\t")
		file_write.write("%s\t%s\t%s\t%s\t%s" % (line[0], line[21], line[110], line[112], line[115]))









stop = timeit.default_timer()
print stop - start 