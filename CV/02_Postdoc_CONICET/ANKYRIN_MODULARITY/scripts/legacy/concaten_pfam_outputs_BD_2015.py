#!/usr/bin/python
# -*- coding: utf-8 -*-

import timeit
start = timeit.default_timer()

from Bio import SeqIO # to parse the fasta file
import sys
import os



my_fasta = sys.argv[1]

num_homologs = int(sys.argv[2])

uniprot_list_total = []
with open(my_fasta, 'rU') as fasta_handle :
	for record in SeqIO.parse(fasta_handle, "fasta") :
		uniprot_id = record.id[3:9]
		uniprot_list_total.append(uniprot_id)


cmd = "mkdir %s_MaxHomologs_MAX%s_homologs_pfam_output" % (my_fasta[:-6],num_homologs)
os.system(cmd)

for uniprot_id in uniprot_list_total:
	print uniprot_id
	with open("%s_MaxHomologs_MAX%s_homologs_pfam_output/pfam_domains_in_%s_MaxHomologs_%s_homolog.txt" % (my_fasta[:-6],num_homologs, uniprot_id, num_homologs), 'a+') as file_write:
		for max_homolog in range (100, num_homologs+100, 100) :
			print max_homolog
			with open("%s_MaxHomologs_%s_homologs_pfam_output/pfam_domains_in_%s_MaxHomologs_%s_homolog.txt" % (my_fasta[:-6],max_homolog, uniprot_id, max_homolog), 'rU') as file_open:
				data = file_open.readlines()
				for line in data :
					file_write.write("%s" % line)



cmd = "mkdir %s_MaxHomologs_MAX%s_queries_pfam_output" % (my_fasta[:-6],num_homologs)
os.system(cmd)

for uniprot_id in uniprot_list_total:
	print uniprot_id
	with open("%s_MaxHomologs_MAX%s_queries_pfam_output/pfam_domains_in_%s_MaxHomologs_%s_query.txt" % (my_fasta[:-6],num_homologs, uniprot_id, num_homologs), 'a+') as file_write:
		
		for max_homolog in range (100, num_homologs+100, 100) :
			print max_homolog
			with open("%s_MaxHomologs_%s_queries_pfam_output/pfam_domains_in_%s_MaxHomologs_%s_query.txt" % (my_fasta[:-6],max_homolog, uniprot_id, max_homolog), 'rU') as file_open:
				data = file_open.readlines()
				for line in data :
					file_write.write("%s" % line)

	


stop = timeit.default_timer()
print stop - start 
