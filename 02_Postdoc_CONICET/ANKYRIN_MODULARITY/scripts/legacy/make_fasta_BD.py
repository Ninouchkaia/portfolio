#!/usr/bin/python
# -*- coding: utf-8 -*-

import timeit
start = timeit.default_timer()

import sys
import os
from collections import Counter
import re
from Bio import SeqIO


my_fasta = sys.argv[1]
uniprot_list_total = []
with open(my_fasta, 'rU') as fasta_handle :
	for record in SeqIO.parse(fasta_handle, "fasta") :
		uniprot_id = record.id[3:9]
		uniprot_list_total.append(uniprot_id)


with open("check_pfam_in_%s" % my_fasta, 'a+') as file_write_query : 
	with open("check_pfam_in_%s_homolog_hits.fasta" % (my_fasta[:-6]), 'a+') as file_write_hit : 

		for uniprot_id in uniprot_list_total :

			domain_dict = {}
			filename = "pfam_in_%s/%s_pfam.txt" % (my_fasta[:-6], uniprot_id)
			if os.path.isfile(filename) == True :
				print uniprot_id

				with open(filename, 'rU') as file_open :
					data = file_open.readlines()
					for line in data :
						line = line.split("\t")
						domain_start = int(line[1])
						domain_end = int(line[2])
						domain_id = line[5]
						if (domain_start, domain_end) not in domain_dict :
							domain_dict[(domain_start,domain_end)] = [domain_id]
						else :
							domain_dict[(domain_start,domain_end)].append(domain_id)

				align_ids = []
				align_boundaries_dict = {}
				with open("/home/nina/scripts/paper_elm_2014/elms/elm_conservation/elm_conservation_data/binding_partners/conservation_identity_30/alignments_30identity/%s_blast.output_alignments_30identity.txt" % uniprot_id, 'rU') as file_open :
					data = file_open.readlines()
						
					for line in data : 
						line = line.split("\t")

						if line[0][:5] == "query" :
							align_id = int(line[0][6:])

							if align_id <= 1000 :
								align_start = int(line[2])
								align_end = int(line[3])


								for (domain_start,domain_end) in domain_dict :
									domain_id = domain_dict[(domain_start,domain_end)]

									if (align_start <= domain_start) and (align_end >= domain_end) :
										align_ids.append(align_id)
										align_boundaries_dict[align_id] = (align_start, align_end)
										seq_query = (line[4].replace("\n","")).replace("-","")
										file_write_query.write(">%s_%s\t%s\t%s\n" % (uniprot_id, align_id, align_start, align_end))
										file_write_query.write("%s\n" % seq_query)
										break
									else :
										continue								


					for line in data :
						line = line.split("\t")

						if (line[0][:3] == "hit"): 
							align_id = int(line[0][4:])

							if align_id in align_ids :
								seq_hit = (line[4].replace("\n","")).replace("-","")
								homolog_id = line[1]
								align_start = align_boundaries_dict[align_id][0]
								align_end = align_boundaries_dict[align_id][1]
								file_write_hit.write(">%s|%s_%s\t%s\t%s\n" % (homolog_id, uniprot_id, align_id, align_start, align_end))
								file_write_hit.write("%s\n" % seq_hit)






stop = timeit.default_timer()
print stop - start

