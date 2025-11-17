#!/usr/bin/python
# -*- coding: utf-8 -*-

import timeit
start = timeit.default_timer()

import sys
import os
from collections import Counter
import re
from Bio import SeqIO
import math


fasta_file_all_ank = sys.argv[1]
fasta_file_BD = sys.argv[2]


# build unp_list from fasta and extract the UNP_IDs correpsonding to ankyrin (from 1234)

unp_from_all_ank, unp_from_BD,  = [],[]
with open(fasta_file_all_ank, 'rU') as fasta_handle :
	for record in SeqIO.parse(fasta_handle, "fasta") :
		unp_id = record.id[3:9]
		unp_from_all_ank.append(unp_id) 

with open(fasta_file_BD, 'rU') as fasta_handle :
	for record in SeqIO.parse(fasta_handle, "fasta") :
		unp_id = record.id[3:9]
		unp_from_BD.append(unp_id) 

unp_from_BD_without_ank = []
for unp_id in unp_from_BD :
	if unp_id in unp_from_all_ank :
		pass
	else :
		unp_from_BD_without_ank.append(unp_id)

print len(unp_from_BD_without_ank) #2028 OK

#build dict of total swissprot counts for all PFAM domains
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

unp_list, domain_list, Pfam_type_list = [], [], []
with open("pfam_in_binding_partners_2038/CONCATEN.txt", 'rU') as file_open : 
	data = file_open.readlines()
	for line in data :
		line = line.split("\t")
		unp_id = line[0]
		if unp_id in unp_from_BD_without_ank :
			unp_list.append(unp_id)
			domain_name = line[6]
			domain_list.append(domain_name)

		
print 'len(unp_list)', len(unp_list)
print "len(set(unp_list))",len(set(unp_list))
print 'len(domain_list)', len(domain_list)
print "len(set(domain_list))",len(set(domain_list))

count_domain_names = Counter(domain_list)
print len(count_domain_names)


num_prot_in_swissprot = 545388
num_prot_in_subfamily = 2038 #len(set(unp_list))
with open("%s_pfam_domains_counts_WITHOUT_ank_IDs_2038.txt" % fasta_file_BD[:-5], 'a+') as file_write :
	file_write.write("Domain name\t# in Ank binding partners\t# in Uniprot\texpected # in Ank binding partners\tlog(obs/exp)\n")
	for domain_name in count_domain_names :
		# if domain_name == 'DUF3424' :
		# 	pass
		# else :
		domain_count_in_subfamily = float(count_domain_names[domain_name])
		domain_count_in_unp = float(Pfam_count_dict[domain_name])
		exp_count_in_subfamily = float(domain_count_in_unp * num_prot_in_subfamily / num_prot_in_swissprot)
		obs_exp = domain_count_in_subfamily/exp_count_in_subfamily
		file_write.write("%s\t%s\t%s\t%s\t%s\n" % (domain_name, domain_count_in_subfamily, domain_count_in_unp, exp_count_in_subfamily, math.log(obs_exp)))


# import glob
# for filename in glob.iglob("%s_pfam_domains_counts_without_ank_IDs.txt" % fasta_file_BD[:-5]) :
# 	print filename
# 	with open(filename, 'rU') as file_read :
# 		with open("%s_replaced.txt" % filename[:-4], 'a+') as file_write :
# 			data = file_read.readlines()
# 			for line in data :
# 				line = line.replace(".",",")
# 				file_write.write("%s" % line)




stop = timeit.default_timer()
print stop - start