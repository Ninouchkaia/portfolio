#!/usr/bin/python
# -*- coding: utf-8 -*-

import timeit
start = timeit.default_timer()

from Bio import SeqIO # to parse the fasta file
import sys
from collections import Counter


fasta_file_all_ank = sys.argv[1]
fasta_file_ank_522 = sys.argv[2]
fasta_file_BD = sys.argv[3]

unp_from_all_ank, unp_from_BD, unp_from_ank_522 = [],[],[]
with open(fasta_file_all_ank, 'rU') as fasta_handle :
	for record in SeqIO.parse(fasta_handle, "fasta") :
		unp_id = record.id[3:9]
		unp_from_all_ank.append(unp_id) 

with open(fasta_file_ank_522, 'rU') as fasta_handle :
	for record in SeqIO.parse(fasta_handle, "fasta") :
		unp_id = record.id[3:9]
		unp_from_ank_522.append(unp_id) 

with open(fasta_file_BD, 'rU') as fasta_handle :
	for record in SeqIO.parse(fasta_handle, "fasta") :
		unp_id = record.id[3:9]
		unp_from_BD.append(unp_id) 


print len(unp_from_all_ank), len(unp_from_ank_522), len(unp_from_BD)



for unp_id in unp_from_BD :
	if unp_id in unp_from_all_ank :
		print unp_id

stop = timeit.default_timer()
print stop - start 