#!/usr/bin/python
# -*- coding: utf-8 -*-

import timeit
start = timeit.default_timer()
import sys
import os
import glob
import re
from commands import getoutput # permet d'obtenir l'output d'une commande bash.
from numpy import prod
# from Bio import SeqIO # to parse the fasta file
import collections
import math





domain_file = sys.argv[1]
elm_file = sys.argv[2]


################################################################################################################################################
####################################################      From the domain list : map elms        ###############################################
################################################################################################################################################

						################################################################################################
						######################      no enrichment threshold     ########################################
						################################################################################################


						# Pour chaque domaine trouvé dans famille domain_file :
						# --> assigner le clans
						# --> assigner l'elm interagissant avec le domaine ou clan ssi il se trouve au moins une fois dans la famille elm_file



# domains_that_have_binding_elms_dict = {} #mapping extended to clans
# with open(domain_file,'rU') as file_open : 
# 	data = file_open.readlines()
# 	for line in data[1:] :
# 		line = line.split("\t")
# 		domain_name = line[0]
# 		domain_or_clan = line[1]
# 		binding_elm = line[2]

# 		if binding_elm != 'NULL':
# 			binding_elms = binding_elm.split(", ")
# 			binding_elms = list(set(binding_elms))
# 			domains_that_have_binding_elms_dict[(domain_name, domain_or_clan)] = binding_elms

# print len(domains_that_have_binding_elms_dict)

# elms_found_in_family = []
# with open(elm_file, 'ru') as file_open :
# 	data = file_open.readlines()
# 	for line in data[1:] :
# 		line = line.split("\t")
# 		elm_name = line[0]
# 		binding_domain_or_clan = line[1].split(", ")
# 		binding_domain_or_clan = list(set(binding_domain_or_clan))
# 		elms_found_in_family.append((elm_name,binding_domain_or_clan))
# print len(elms_found_in_family)


# with open("matched_elm_domain_pairs.txt", 'a+') as file_write :
# 	file_write.write("%s\t%s\t%s\t%s\n" % ("domain_found_in_family", "corresponding_clan", "binding_elm_found_in_counter_family", "domain_or_clan_to_which_this_elm_binds"))
# 	for domain_tuple in domains_that_have_binding_elms_dict :
# 		for elm_name in domains_that_have_binding_elms_dict[domain_tuple] :
# 			for elm_tuple in elms_found_in_family :
# 				if elm_name in elm_tuple :
# 					domain_found_in_family = domain_tuple[0]
# 					corresponding_clan = domain_tuple[1]
# 					binding_elm_found_in_counter_family = elm_tuple[0]
# 					domain_or_clan_to_which_this_elm_binds = elm_tuple[1]

# 					file_write.write("%s\t%s\t%s\t%s\n" % (domain_found_in_family, corresponding_clan, binding_elm_found_in_counter_family, domain_or_clan_to_which_this_elm_binds))




# 						################################################################################################
# 						######################      enrichment threshold : Zscore > 0    ###############################
# 						################################################################################################
						

# 						# Pour chaque domaine avec Zscore > 0 trouvé dans famille domain_file :
# 						# --> assigner le clans
# 						# --> assigner l'elm interagissant avec le domaine ou clan 
# 						# ssi il se trouve dans la famille elm_file
# 						# et ssi il a un Zscore > 0


# domains_that_have_binding_elms_dict = {} #mapping extended to clans
# with open(domain_file,'rU') as file_open : 
# 	data = file_open.readlines()
# 	for line in data[1:] :
# 		line = line.split("\t")
# 		domain_name = line[0]
# 		domain_or_clan = line[1]
# 		binding_elm = line[2]
# 		Zscore = float(line[4])

# 		if Zscore > 0 :

# 			if binding_elm != 'NULL':
# 				binding_elms = binding_elm.split(", ")
# 				binding_elms = list(set(binding_elms))
# 				domains_that_have_binding_elms_dict[(domain_name, domain_or_clan)] = binding_elms

# print len(domains_that_have_binding_elms_dict)


# elms_found_in_family = []
# with open(elm_file, 'ru') as file_open :
# 	data = file_open.readlines()
# 	for line in data[1:] :
# 		line = line.split("\t")
# 		elm_name = line[0]
# 		Zscore = float(line[8])

# 		if Zscore > 0 :

# 			binding_domain_or_clan = line[1].split(", ")
# 			binding_domain_or_clan = list(set(binding_domain_or_clan))
# 			elms_found_in_family.append((elm_name,binding_domain_or_clan))

# print len(elms_found_in_family)


# with open("matched_elm_domain_pairs.txt", 'a+') as file_write :
# 	file_write.write("%s\t%s\t%s\t%s\n" % ("domain_found_in_family", "corresponding_clan", "binding_elm_found_in_counter_family", "domain_or_clan_to_which_this_elm_binds"))
# 	for domain_tuple in domains_that_have_binding_elms_dict :
# 		for elm_name in domains_that_have_binding_elms_dict[domain_tuple] :
# 			for elm_tuple in elms_found_in_family :
# 				if elm_name in elm_tuple :
# 					domain_found_in_family = domain_tuple[0]
# 					corresponding_clan = domain_tuple[1]
# 					binding_elm_found_in_counter_family = elm_tuple[0]
# 					domain_or_clan_to_which_this_elm_binds = elm_tuple[1]

# 					file_write.write("%s\t%s\t%s\t%s\n" % (domain_found_in_family, corresponding_clan, binding_elm_found_in_counter_family, domain_or_clan_to_which_this_elm_binds))




###############################################################################################################################################
###################################################      From the elm list : map domains       ###############################################
###############################################################################################################################################

						###############################################################################################
						#####################      no enrichment threshold     ########################################
						###############################################################################################


						# Pour chaque elm trouvé dans famille elm_file :
						# --> assigner le domain interagissant avec l'elm ssi il se trouve au moins une fois dans la famille domain_file


# elm_domain_mapping_dict = {}
# with open("elm_interaction_domains_modified.txt", 'rU') as file_open :
# 	data = file_open.readlines()
# 	for line in data :
# 		line = line.split("\t")
# 		elm_name = line[0]
# 		domain_name = line[2]
# 		if elm_name not in elm_domain_mapping_dict :
# 			elm_domain_mapping_dict[elm_name] = [domain_name]
# 		else : 
# 			elm_domain_mapping_dict[elm_name].append(domain_name)


# elms_that_have_binding_domain_or_clan_dict = {} 
# with open(elm_file,'rU') as file_open : 
# 	data = file_open.readlines()
# 	for line in data[1:] :
# 		line = line.split("\t")
# 		elm_name = line[0]
# 		binding_domain_or_clan = line[1] 

# 		if elm_name in elm_domain_mapping_dict :
# 			binding_domain_only = elm_domain_mapping_dict[elm_name] #liste!
# 		else :
# 			binding_domain_only = ['NULL']

# 		if binding_domain_or_clan != 'NULL':
# 			binding_domains_or_clans = binding_domain_or_clan.split(", ")
# 			binding_domains_or_clans = list(set(binding_domains_or_clans))
# 			binding_domains_only = list(set(binding_domain_only))
# 			elms_that_have_binding_domain_or_clan_dict[elm_name] = (binding_domains_only, binding_domains_or_clans)

# print len(elms_that_have_binding_domain_or_clan_dict)

# # for i in elms_that_have_binding_domain_or_clan_dict :
# # 	print i, elms_that_have_binding_domain_or_clan_dict[i]

# domains_found_in_family = {}
# with open(domain_file, 'ru') as file_open :
# 	data = file_open.readlines()
# 	for line in data[1:] :
# 		line = line.split("\t")
# 		domain_name = line[0]
# 		corresponding_domain_or_clan = line[1]
# 		binding_elms = line[2].split(", ")
# 		binding_elms = list(set(binding_elms))
# 		domains_found_in_family[domain_name] = (corresponding_domain_or_clan, binding_elms)
# print len(domains_found_in_family)

# # for i in domains_found_in_family :
# # 	print i, domains_found_in_family[i]


# with open("matched_domain_elm_pairs.txt", 'a+') as file_write :
# 	file_write.write("%s\t%s\t%s\t%s\n" % ("elm_found_in_family", "binding_domain_found_in_counter_family", "corresponding_clan", "elms_to_which_this_domain_or_clan_binds"))
	
# 	for elm_name in elms_that_have_binding_domain_or_clan_dict :
# 		# print "elm_name", elm_name
		
# 		binding_domains_only = elms_that_have_binding_domain_or_clan_dict[elm_name][0]
# 		# binding_domains_or_clans = elms_that_have_binding_domain_or_clan_dict[elm_name][1]
# 		# print "binding_domains_only", binding_domains_only
# 		# print "binding_domains_or_clans", binding_domains_or_clans

# 		for domain_name in binding_domains_only :
			
# 			if domain_name in domains_found_in_family :
			
# 				elm_found_in_family = elm_name
# 				binding_domain_found_in_counter_family = domain_name
# 				# print "domain_name", domain_name

# 				corresponding_clan = domains_found_in_family[domain_name][0]
# 				# print "corresponding_clan_bis", corresponding_clan_bis

# 				elms_to_which_this_domain_or_clan_binds = domains_found_in_family[domain_name][1]
# 				# print "elms_to_which_this_domain_or_clan_binds", elms_to_which_this_domain_or_clan_binds

# 				file_write.write("%s\t%s\t%s\t%s\n" % (elm_found_in_family, binding_domain_found_in_counter_family, corresponding_clan, elms_to_which_this_domain_or_clan_binds))



				




# 						################################################################################################
# 						######################      enrichment threshold : Zscore > 0    ###############################
# 						################################################################################################
						

# 						# Pour chaque elm avec Zscore > 0 trouvé dans famille elm_file :
# 						# --> assigner le domain interagissant avec l'elm
# 						# ssi il se trouve dans la famille domain_file
# 						# et ssi il a un Zscore > 0






elm_domain_mapping_dict = {}
with open("elm_interaction_domains_modified.txt", 'rU') as file_open :
	data = file_open.readlines()
	for line in data :
		line = line.split("\t")
		elm_name = line[0]
		domain_name = line[2]
		if elm_name not in elm_domain_mapping_dict :
			elm_domain_mapping_dict[elm_name] = [domain_name]
		else : 
			elm_domain_mapping_dict[elm_name].append(domain_name)


elms_that_have_binding_domain_or_clan_dict = {} 
with open(elm_file,'rU') as file_open : 
	data = file_open.readlines()
	for line in data[1:] :
		line = line.split("\t")
		elm_name = line[0]
		binding_domain_or_clan = line[1] 
		Zscore = float(line[8])

		if Zscore > 0 :

			if elm_name in elm_domain_mapping_dict :
				binding_domain_only = elm_domain_mapping_dict[elm_name] #liste!
			else :
				binding_domain_only = ['NULL']

			if binding_domain_or_clan != 'NULL':
				binding_domains_or_clans = binding_domain_or_clan.split(", ")
				binding_domains_or_clans = list(set(binding_domains_or_clans))
				binding_domains_only = list(set(binding_domain_only))
				elms_that_have_binding_domain_or_clan_dict[elm_name] = (binding_domains_only, binding_domains_or_clans)

print len(elms_that_have_binding_domain_or_clan_dict)

# for i in elms_that_have_binding_domain_or_clan_dict :
# 	print i, elms_that_have_binding_domain_or_clan_dict[i]

domains_found_in_family = {}
with open(domain_file, 'ru') as file_open :
	data = file_open.readlines()
	for line in data[1:] :
		line = line.split("\t")
		domain_name = line[0]
		corresponding_domain_or_clan = line[1]
		Zscore = float(line[4])

		if Zscore > 0 :
			binding_elms = line[2].split(", ")
			binding_elms = list(set(binding_elms))
			domains_found_in_family[domain_name] = (corresponding_domain_or_clan, binding_elms)
print len(domains_found_in_family)

# for i in domains_found_in_family :
# 	print i, domains_found_in_family[i]


with open("matched_domain_elm_pairs_Zscore_positive.txt", 'a+') as file_write :
	file_write.write("%s\t%s\t%s\t%s\n" % ("elm_found_in_family", "binding_domain_found_in_counter_family", "corresponding_clan", "elms_to_which_this_domain_or_clan_binds"))
	
	for elm_name in elms_that_have_binding_domain_or_clan_dict :
		# print "elm_name", elm_name
		
		binding_domains_only = elms_that_have_binding_domain_or_clan_dict[elm_name][0]
		# binding_domains_or_clans = elms_that_have_binding_domain_or_clan_dict[elm_name][1]
		# print "binding_domains_only", binding_domains_only
		# print "binding_domains_or_clans", binding_domains_or_clans

		for domain_name in binding_domains_only :
			
			if domain_name in domains_found_in_family :
			
				elm_found_in_family = elm_name
				binding_domain_found_in_counter_family = domain_name
				# print "domain_name", domain_name

				corresponding_clan = domains_found_in_family[domain_name][0]
				# print "corresponding_clan_bis", corresponding_clan_bis

				elms_to_which_this_domain_or_clan_binds = domains_found_in_family[domain_name][1]
				# print "elms_to_which_this_domain_or_clan_binds", elms_to_which_this_domain_or_clan_binds

				file_write.write("%s\t%s\t%s\t%s\n" % (elm_found_in_family, binding_domain_found_in_counter_family, corresponding_clan, elms_to_which_this_domain_or_clan_binds))















































# filename = sys.argv[1]
# liste = []
# with open(filename, 'rU') as file_open :
# 	data = file_open.readlines()
# 	for line in data[1:] :
# 		line = line.split("\t")
# 		name = line[0]
# 		liste.append(name)
# liste = list(set(liste))
# print len(liste)




















































stop = timeit.default_timer()
print stop - start 