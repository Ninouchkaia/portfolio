#!/usr/bin/python
# -*- coding: utf-8 -*-

import timeit # pour mesurer le temps que met le script a courir
start = timeit.default_timer()
import sys
import re
from commands import getoutput # permet d'obtenir l'output d'une commande bash.


#creation des listes

mapped_list1 = open("/home/nina/scripts/paxdb/%s" % sys.argv[1])
#for example : 'mouse_mapping.txt_results.txt'
mapped_list1 = mapped_list1.readlines()

mapped_list2 = open("/home/nina/scripts/paxdb/%s" % sys.argv[2])
#for example : 'mouse_mapping.txt_results2.txt'
mapped_list2 = mapped_list2.readlines()

file_write = open("/home/nina/scripts/paxdb/temp_mapping.txt", 'a+')

i=0
a=0
b=0
for identifier1 in mapped_list1 :
	a=a+1
	#print "START"
	#print identifier1
	identifier1 = identifier1.split("\t")
	string_id = identifier1[0]
	to_be_mapped = identifier1[1]
	to_be_mapped = to_be_mapped.replace('\n', '') # remove '\n' only
	#print to_be_mapped
	#print string_id, to_be_mapped
	for identifier2 in mapped_list2 :
		b=b+1
		#print identifier2
		identifier2 = identifier2.split("\t")
		mapped_by_uniprot = identifier2[0]
		#print mapped_by_uniprot
		if to_be_mapped[0:6] == mapped_by_uniprot[0:6] :
		#if to_be_mapped == mapped_by_uniprot:
			print "%s\t%s\t%s\t%s" % (string_id, to_be_mapped, mapped_by_uniprot, identifier2[1])
			i=i+1
			file_write.write("%s\t%s\t%s\t%s" % (string_id, to_be_mapped, mapped_by_uniprot, identifier2[1]))
	#print "END"

print i
print a 
print b

stop = timeit.default_timer()
print stop - start 