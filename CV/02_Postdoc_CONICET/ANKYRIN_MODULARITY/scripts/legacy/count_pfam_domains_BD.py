#!/usr/bin/python
# -*- coding: utf-8 -*-

import timeit
start = timeit.default_timer()

import re
from Bio import SeqIO # to parse the fasta file
from commands import getoutput # permet d'obtenir l'output d'une commande bash.
import sys
import os
import glob
from collections import Counter


my_fasta = sys.argv[1]
uniprot_list_total = []
with open(my_fasta, 'rU') as fasta_handle :
	for record in SeqIO.parse(fasta_handle, "fasta") :
		uniprot_id = record.id[3:9]
		uniprot_list_total.append(uniprot_id)




unp_list, domain_list, Pfam_type_list = [], [], []
for uniprot_id in uniprot_list_total :
	filename = "pfam_in_binding_partners_2038/%s_pfam.txt" % uniprot_id
	if os.path.isfile(filename) == True :
		with open(filename, 'rU') as file_open : 
			data = file_open.readlines()
			for line in data :
				line = line.split("\t")
				domain_name = line[6]
				domain_list.append(domain_name)
				unp_id = line[0]
				unp_list.append(unp_id)
		
print 'len(unp_list)', len(unp_list)
print "len(set(unp_list))",len(set(unp_list))
print 'len(domain_list)', len(domain_list)
print "len(set(domain_list))",len(set(domain_list))

count_domain_names = Counter(domain_list)
# print count_domain_names

# count_Pfam_types = Counter(Pfam_type_list)
# for i in count_Pfam_types :
# 	print count_Pfam_types[i], i, "matches, " 


#build dict of total swissprot counts for all SMART domains
Pfam_count_dict, Pfam_description_dict = {},{}
total = 0
with open("FrequenciesPfam_mapped.txt", 'rU') as file_read :
	data = file_read.readlines()
	for line in data[1:] :
		line = line.split("\t")
		domain = line[1]
		description = line[-2]
		domain_count_in_swissprot = (line[-1].replace("\n",""))
		Pfam_count_dict[domain] = float(domain_count_in_swissprot)
		total = total + Pfam_count_dict[domain]
		Pfam_description_dict[domain] = description
print len(Pfam_count_dict) #14'831
print total #28738352.0

# num_prot_in_swissprot = 18523877
# num_prot_in_subfamily = 510
# with open("hypergeometric_test_Ank_proteins_510_3.txt", 'a+') as file_write :
# 	file_write.write("Domain name\t# in Ank proteins\t# in Uniprot\tp-value <\tp-value >\n")
# 	for domain_name in count_domain_names :
# 		if domain_name == 'DUF3424' :
# 			pass
# 		else :
# 			domain_count_in_subfamily = float(count_domain_names[domain_name])
# 			domain_count_in_unp = float(Pfam_count_dict[domain_name])
# 			arg3 = float(num_prot_in_swissprot-domain_count_in_unp)
# 			p_val_inf = str(stats.hypergeom.cdf(domain_count_in_subfamily + 1 , domain_count_in_unp, arg3, num_prot_in_subfamily))
# 			p_val_sup = str(stats.hypergeom.cdf(domain_count_in_subfamily - 1 , domain_count_in_unp, arg3 , num_prot_in_subfamily))
# 			file_write.write("%s\t%s\t%s\t%s\t%s\n" % (domain_name, domain_count_in_subfamily, domain_count_in_unp, p_val_inf, p_val_sup))

import math

num_prot_in_swissprot = 18523877
num_prot_in_subfamily = 2038
with open("Pfam_domains_in_BD_2038.txt", 'a+') as file_write :
	file_write.write("Domain name\t# in Binding Partners\t# in Uniprot\texpected # in Binding Partners\tlog(obs/exp)\n")
	for domain_name in count_domain_names :
		# if domain_name == 'DUF3424' :
		# 	pass
		# else :
		domain_count_in_subfamily = float(count_domain_names[domain_name])
		domain_count_in_unp = float(Pfam_count_dict[domain_name])
		exp_count_in_BD = float(domain_count_in_unp * num_prot_in_subfamily / num_prot_in_swissprot)
		obs_exp = domain_count_in_subfamily/exp_count_in_BD
		file_write.write("%s\t%s\t%s\t%s\t%s\n" % (domain_name, domain_count_in_subfamily, domain_count_in_unp, exp_count_in_BD, math.log(obs_exp)))
import glob
for filename in glob.iglob("Pfam_domains_in_BD_2038.txt") :
	print filename
	with open(filename, 'rU') as file_read :
		with open("%s_replaced.txt" % filename[:-4], 'a+') as file_write :
			data = file_read.readlines()
			for line in data :
				line = line.replace(".",",")
				file_write.write("%s" % line)




stop = timeit.default_timer()
print stop - start