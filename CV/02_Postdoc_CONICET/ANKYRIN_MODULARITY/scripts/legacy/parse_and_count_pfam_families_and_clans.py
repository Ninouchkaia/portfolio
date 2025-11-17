#!/usr/bin/python
# -*- coding: utf-8 -*-

import timeit
start = timeit.default_timer()

import sys
import os
from collections import Counter
import re
import glob


pfam_file = sys.argv[1]


# unp_list, domain_list, clan_list = [], [], []

# for filename in glob.iglob("pfam_in_%s/*.txt" % pfam_file[:-6]) :
# 	print filename
# 	with open(filename, 'rU') as file_open :
# 		data = file_open.readlines()
# 		for line in data :
# 			line = line.split("\t")
# 			domain_name = line[6]
# 			domain_list.append(domain_name)
# 			unp_id = line[0]
# 			unp_list.append(unp_id)
# 			clan_name = line[14].replace("\n","")
# 			clan_list.append(clan_name)
		
# print 'len(unp_list)', len(unp_list)
# print "len(set(unp_list))",len(set(unp_list))
# print 'len(domain_list)', len(domain_list)
# print "len(set(domain_list))",len(set(domain_list))
# print 'len(clan_list)', len(clan_list)
# print "len(set(clan_list))",len(set(clan_list))


unp_list, domain_list, clan_list = [], [], []

for filename in glob.iglob("pfam_in_%s/*.txt" % pfam_file[:-6]) :
	print filename
	with open(filename, 'rU') as file_open :
		data = file_open.readlines()
		for line in data :
			line = line.split("\t")
			clan_name = line[14].replace("\n","")
			if clan_name != "No_clan" :
				clan_list.append(clan_name)
			else :
				domain_name = line[6]
				domain_list.append(domain_name)
			unp_id = line[0]
			unp_list.append(unp_id)

		
print 'len(unp_list)', len(unp_list)
print "len(set(unp_list))",len(set(unp_list))
print 'len(domain_list)', len(domain_list)
print "len(set(domain_list))",len(set(domain_list))
print 'len(clan_list)', len(clan_list)
print "len(set(clan_list))",len(set(clan_list))


count_domain_names = Counter(domain_list)

count_clan_names = Counter(clan_list)

count_domain_and_clan_names = dict(count_domain_names.items() + count_clan_names.items())


print len(count_domain_names)
print len(count_clan_names)
print len(count_domain_and_clan_names)





#build dict of total swissprot counts for all PFAM domains
Pfam_count_dict = {}
with open("FrequenciesPfam_mapped_CLANS.txt", 'rU') as file_read :
	data = file_read.readlines()
	for line in data :
		line = line.split("\t")
		domain = line[0]
		domain_count_in_swissprot = (line[-1].replace("\n",""))
		Pfam_count_dict[domain] = float(domain_count_in_swissprot)
print len(Pfam_count_dict) #10784







import math
num_prot_in_swissprot = 545388
num_prot_in_subfamily = len(set(unp_list))

with open("Pfam_clans_and_families_counts_in_%s.txt" % pfam_file[:-6], 'a+') as file_write :
	file_write.write("Domain/CLAN name\t# in Ank binding partners\t# in Uniprot\texpected # in Ank binding partners\tlog(obs/exp)\n")
	for domain_name in count_domain_and_clan_names :
		domain_count_in_subfamily = float(count_domain_and_clan_names[domain_name])
		domain_count_in_unp = float(Pfam_count_dict[domain_name])
		exp_count_in_subfamily = float(domain_count_in_unp * num_prot_in_subfamily / num_prot_in_swissprot)
		obs_exp = domain_count_in_subfamily/exp_count_in_subfamily
		file_write.write("%s\t%s\t%s\t%s\t%s\n" % (domain_name, domain_count_in_subfamily, domain_count_in_unp, exp_count_in_subfamily, math.log(obs_exp)))


# import glob
# for filename in glob.iglob("%s_pfam_domains_counts.txt" % pfam_file[:-5]) :
# 	print filename
# 	with open(filename, 'rU') as file_read :
# 		with open("%s_replaced.txt" % filename[:-4], 'a+') as file_write :
# 			data = file_read.readlines()
# 			for line in data :
# 				line = line.replace(".",",")
# 				file_write.write("%s" % line)











stop = timeit.default_timer()
print stop - start