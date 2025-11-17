#!/usr/bin/python
# -*- coding: utf-8 -*-

import timeit
start = timeit.default_timer()

from Bio import SeqIO
import sys


ank_proteins_from_uniprot_fasta = "all-ank-20130926.fasta"
sequence_dict_from_fasta = {}
ank_from_fasta_list_total = []
with open(ank_proteins_from_uniprot_fasta, 'rU') as fasta_handle :
	for record in SeqIO.parse(fasta_handle, "fasta") :
		uniprot_id = record.id[3:9]
		ank_from_fasta_list_total.append(uniprot_id)
		uniprot_seq = str(record.seq)
		sequence_dict_from_fasta[uniprot_id] = uniprot_seq

ank_proteins_from_rocio_fasta = "uniprotAnkUnique.fasta"
sequence_dict_from_rocio_fasta = {}
ank_from_rocio_fasta_list_total = []
with open(ank_proteins_from_rocio_fasta, 'rU') as fasta_handle :
	for record in SeqIO.parse(fasta_handle, "fasta") :
		uniprot_id = record.id[3:9]
		ank_from_rocio_fasta_list_total.append(uniprot_id)
		uniprot_seq = str(record.seq)
		sequence_dict_from_rocio_fasta[uniprot_id] = uniprot_seq

print len(ank_from_fasta_list_total), len(ank_from_rocio_fasta_list_total)


for uniprot_id in ank_from_fasta_list_total :
	if uniprot_id not in ank_from_rocio_fasta_list_total :
		print uniprot_id
		with open("check_ank_presence_24.fasta", 'a+') as file_write :
			sequence = sequence_dict_from_fasta[uniprot_id]
			file_write.write(">sp|%s\n" % uniprot_id)
			file_write.write("%s\n" % sequence)



stop = timeit.default_timer()
print stop - start 
