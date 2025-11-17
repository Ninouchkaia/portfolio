#!/usr/bin/python
# -*- coding: utf-8 -*-

import timeit
start = timeit.default_timer()

import sys
import os
import re
from Bio import SeqIO
import subprocess


my_fasta = sys.argv[1]

threshold = int(sys.argv[2])

uniprot_list_total = []
with open(my_fasta, 'rU') as fasta_handle :
	for record in SeqIO.parse(fasta_handle, "fasta") :
		uniprot_id = record.id[3:9]
		uniprot_list_total.append(uniprot_id)


cmd1 = "mkdir %s_MaxHomologs_%s_queries_pfam_output" % (my_fasta[:-6], threshold)
cmd2 = "mkdir %s_MaxHomologs_%s_homologs_pfam_output" % (my_fasta[:-6], threshold)
os.system(cmd1)
os.system(cmd2)


for uniprot_id in uniprot_list_total :

	print uniprot_id
	sys.stdout.flush() # http://stackoverflow.com/questions/8537932/why-is-my-nohup-out-empty

	query_fasta = "%s_MaxHomologs_%s_queries_fastas/check_pfam_in_%s_MaxHomologs_%s_query_seq.fasta" % (my_fasta[:-6], threshold, uniprot_id, threshold)
	homolog_fasta = "%s_MaxHomologs_%s_homologs_fastas/check_pfam_in_%s_MaxHomologs_%s_homolog_hit_seq.fasta" % (my_fasta[:-6], threshold, uniprot_id, threshold)

	nohup_query_output = "%s_MaxHomologs_%s_queries_pfam_output/pfam_domains_in_%s_MaxHomologs_%s_query.txt" % (my_fasta[:-6], threshold, uniprot_id, threshold)
	nohup_homolog_output = "%s_MaxHomologs_%s_homologs_pfam_output/pfam_domains_in_%s_MaxHomologs_%s_homolog.txt" % (my_fasta[:-6], threshold, uniprot_id, threshold)

	cmd_query = "nohup perl pfam_scan.pl -cpu 4 -fasta %s -dir /home/tannat/Desktop/Nina/domain_retrieve/Pfam_local_search/PfamScan/binding_partners &" % (query_fasta)
	# cmd_query = cmd_query.split(" ")

	cmd_homolog = "nohup perl pfam_scan.pl -cpu 4 -fasta %s -dir /home/tannat/Desktop/Nina/domain_retrieve/Pfam_local_search/PfamScan/binding_partners &" % (homolog_fasta)
	# cmd_homolog = cmd_homolog.split(" ")
	
	pfam_scan_query_output = subprocess.Popen(cmd_query, stderr=subprocess.STDOUT, stdout=subprocess.PIPE, shell=True)
	stdoutdata_query = pfam_scan_query_output.communicate()[0]
	with open(nohup_query_output, 'a+') as file_write :
		file_write.write("%s" % stdoutdata_query)

	pfam_scan_homolog_output = subprocess.Popen(cmd_homolog, stderr=subprocess.STDOUT, stdout=subprocess.PIPE, shell=True)
	stdoutdata_homolog = pfam_scan_homolog_output.communicate()[0]
	with open(nohup_homolog_output, 'a+') as file_write :
		file_write.write("%s" % stdoutdata_homolog)






stop = timeit.default_timer()
print stop - start

