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
from collections import Counter, defaultdict


my_fasta = sys.argv[1]

num_homologs = int(sys.argv[2])

uniprot_list_total = []
with open(my_fasta, 'rU') as fasta_handle :
	for record in SeqIO.parse(fasta_handle, "fasta") :
		uniprot_id = record.id[3:9]
		uniprot_list_total.append(uniprot_id)

# uniprot_list_total = uniprot_list_total[34:56]

# for filename in glob.iglob("%s_MaxHomologs_100_queries_pfam_output/*" % my_fasta[:-6]) :
# 	uniprot_id = filename[-32:-26]






super_averaged_conservation_dict = {}
domain_name_list_total = []
for uniprot_id in uniprot_list_total :
	print uniprot_id
	filename_query = "%s_MaxHomologs_MAX%s_queries_pfam_output/pfam_domains_in_%s_MaxHomologs_%s_query.txt" % (my_fasta[:-6], num_homologs, uniprot_id, num_homologs)
	filename_homolog = "%s_MaxHomologs_MAX%s_homologs_pfam_output/pfam_domains_in_%s_MaxHomologs_%s_homolog.txt" % (my_fasta[:-6], num_homologs, uniprot_id, num_homologs)
	# print uniprot_id, "query"
	query_domain_dict = {}
	domain_name_list = []

	with open(filename_query, 'rU') as file_open :
		data = file_open.readlines()

		for line in data :
			if line[0] == "#" or line == "\n" or "FATAL" in line:
				pass
			else :
				line = line.split()
				homolog_ref = int(line[0][7:])
				domain_name = line[6]
				domain_name_list.append(domain_name)
				domain_name_list_total.append(domain_name)
				domain_start = line[1]
				domain_end = line[2]
				if homolog_ref not in query_domain_dict :
					query_domain_dict[homolog_ref] = [(domain_name, domain_start, domain_end)]
				else :
					query_domain_dict[homolog_ref].append((domain_name, domain_start, domain_end))
		# count = sum(len(domain_names) for domain_names in query_domain_dict.itervalues())				
		# print count, "domains in query"
		domain_name_list = list(set(domain_name_list))


	# print uniprot_id, "homologs"
	homolog_domain_dict = {}
	with open(filename_homolog, 'rU') as file_open :
		data = file_open.readlines()
		for line in data :
			if line[0] == "#" or line == "\n" or "FATAL" in line:
				pass
			else :
				line = line.split()
				homolog_ref = int((line[0].split("|"))[-1][7:])
				domain_name = line[6]
				domain_start = line[1]
				domain_end = line[2]
				if homolog_ref not in homolog_domain_dict :
					homolog_domain_dict[homolog_ref] = [(domain_name, domain_start, domain_end)]
				else :
					homolog_domain_dict[homolog_ref].append((domain_name, domain_start, domain_end))
	# count = sum(len(domain_names) for domain_names in homolog_domain_dict.itervalues())				
	# print count, "domains in homologs"

	
	with open("conservation_state_in_BD_%s.txt" % num_homologs, 'a+') as file_write :

		partial_conservations_dict = {}

		

		for homolog_ref in query_domain_dict :

			if homolog_ref not in homolog_domain_dict : #aucun domaine n'a été détecté dans la sequence homologue de l'alignement
				print "no Pfam domain were found in homolog of %s alignment # %s" % (uniprot_id, homolog_ref)
				for (domain_name, domain_start, domain_end) in query_domain_dict[homolog_ref] :
					file_write.write("%s\t%s\t%s\t%s\t%s\tnot_conserved\n" % (uniprot_id, homolog_ref, domain_name, domain_start, domain_end))
			
			elif homolog_ref in homolog_domain_dict : #at least one domain was detected in the corresponding alignment	
				count_domains_in_query_dict = Counter([domain_info[0] for domain_info in query_domain_dict[homolog_ref]])
				count_domains_in_homolog_dict = Counter([domain_info[0] for domain_info in homolog_domain_dict[homolog_ref]])
				# print homolog_ref
				# print query_domain_dict[homolog_ref]
				# print count_domains_in_query_dict
				# print homolog_domain_dict[homolog_ref]
				# print count_domains_in_homolog_dict			
				for domain_name in count_domains_in_query_dict :
					if domain_name in count_domains_in_homolog_dict :
						
						if int(count_domains_in_query_dict[domain_name]) <= int(count_domains_in_homolog_dict[domain_name]) :
							for (domain_name, domain_start, domain_end) in query_domain_dict[homolog_ref] :
								file_write.write("%s\t%s\t%s\t%s\t%s\tconserved\n" % (uniprot_id, homolog_ref, domain_name, domain_start, domain_end))
					
						elif int(count_domains_in_query_dict[domain_name]) > int(count_domains_in_homolog_dict[domain_name]) :
							n = int(count_domains_in_query_dict[domain_name]) - int(count_domains_in_homolog_dict[domain_name])
							m = int(count_domains_in_homolog_dict[domain_name])
							file_write.write("%s\t%s\t%s\t\t\tnot_conserved\n" % (uniprot_id, homolog_ref, domain_name) * n)
							file_write.write("%s\t%s\t%s\t\t\t\t\tconserved\n" % (uniprot_id, homolog_ref, domain_name) * m)




domain_name_list_total = []
with open("conservation_state_in_BD_%s.txt" % num_homologs, 'rU') as file_open :
	matched_data = file_open.readlines()
	for line in matched_data :
		line = line.split("\t")
		domain_name = line[2]
		domain_name_list_total.append(domain_name)
domain_name_list_total = list(set(domain_name_list_total))

conserved, not_conserved = {}, {}
for domain_name in domain_name_list_total :
	conserved[domain_name] = 0
	not_conserved[domain_name] = 0

with open("conservation_state_in_BD_%s.txt" % num_homologs, 'rU') as file_open :
	matched_data = file_open.readlines()
	for line in matched_data :
		line = line.split("\t")
		domain_name = line[2]
		conservation_state = line[-1].replace("\n","")
		if conservation_state == "conserved" :
			conserved[domain_name] += 1
		elif conservation_state == "not_conserved" :
			not_conserved[domain_name] += + 1
		else : 
			print "weird conservation state", conservation_state

conservation_degree_dict = {}
with open("domain_conservation_percentages_in_BD_%s.txt" % num_homologs, 'a+') as file_write :
	for domain_name in domain_name_list_total :

				if (conserved[domain_name]+not_conserved[domain_name]) == 0 :
					print "the domain %s was NEVER found in Binding Partners" % (domain_name)
					# pass
				else :							
					percentage = float(conserved[domain_name])/float(conserved[domain_name]+not_conserved[domain_name])
					conservation_degree_dict[domain_name] = percentage
					
					file_write.write("%s\t%.5f\n" % (domain_name, percentage))







	


stop = timeit.default_timer()
print stop - start 
