#!/usr/bin/python
# -*- coding: utf-8 -*-
# jupyter notebook --browser='C:/Program Files (x86)/Google/Chrome/Application/chrome.exe %s'


import timeit
start = timeit.default_timer()

from scipy import stats
import numpy as np
from scipy.stats import shapiro,norm,normaltest
from scipy.special import erf, ndtr
import math
import os


# terms = ['Disease','Drug','GO','Prot','Symptom']
# outputFiles
# for filename in os.listdir("A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\Further_Analysis\\full_analysis_for_source_to_proteins_undirected\\pvalues_from_zscores\\outputFiles"):
# 	print(filename)
# 	with open("A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\Further_Analysis\\full_analysis_for_source_to_proteins_undirected\\pvalues_from_zscores\\outputFiles\\%s" % filename, 'r') as file_read :
# 		data = file_read.readlines()
# 		with open("A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\Further_Analysis\\full_analysis_for_source_to_proteins_undirected\\pvalues_from_zscores\\outputFiles\\%s_manlio_pval.tsv" % filename[:-4], 'a+') as file_write :
# 			file_write.write(data[0])
# 			for line in data[1:] :
# 				line = line.split('\t')
# 				my_node = line[0]
# 				sample_size = float(line[1])	
# 				min_deg = float(line[2])
# 				max_deg = float(line[3])
# 				mean_deg = float(line[4])
# 				med_deg = float(line[5])	
# 				sd = float(line[6])	
# 				zscore = float(line[7])
# 				deg_covid = float(line[8])	
				
# 				isNormalShapiro	= line[9]
# 				print(isNormalShapiro)
# 				if isNormalShapiro == 'True' :
# 					p_value_shapiro = float(line[10])
# 				else : 
# 					p_value_shapiro = 1 / (zscore*zscore)	
				
# 				isNormalDagostino = line[11]
# 				if isNormalDagostino == 'True':	
# 					p_value_dagostino = float(line[12].replace('\n',''))
# 				else :
# 					p_value_dagostino = 1 / (zscore*zscore)
				

# 				file_write.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (my_node,sample_size,min_deg,max_deg,mean_deg,med_deg,sd,zscore,deg_covid,isNormalShapiro,p_value_shapiro,isNormalDagostino,p_value_dagostino))


# for filename in os.listdir("A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\Further_Analysis\\full_analysis_for_source_to_proteins_undirected\\pvalues_from_zscores\\outputFiles\\abs"):
# 	print(filename)
# 	with open("A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\Further_Analysis\\full_analysis_for_source_to_proteins_undirected\\pvalues_from_zscores\\outputFiles\\%s" % filename, 'r') as file_read :
		
# 		with open("A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\Further_Analysis\\full_analysis_for_source_to_proteins_undirected\\pvalues_from_zscores\\outputFiles\\abs\\aggregated_nodes_abs.tsv", 'a+') as file_write :
# 			data = file_read.readlines()
# 			# file_write.write(data[0])
# 			for line in data[1:] :
# 				line = line.split('\t')
# 				my_node = line[0]
# 				sample_size = float(line[1])	
# 				min_deg = float(line[2])
# 				max_deg = float(line[3])
# 				mean_deg = float(line[4])
# 				med_deg = float(line[5])	
# 				sd = float(line[6])	
# 				zscore = float(line[7])
# 				deg_covid = float(line[8])	
				
# 				isNormalShapiro	= line[9]
# 				# print(isNormalShapiro)
# 				if isNormalShapiro == 'True' :
# 					p_value_shapiro = float(line[10])
# 				else : 
# 					p_value_shapiro = 1 / (zscore*zscore)	
				
# 				isNormalDagostino = line[11]
# 				if isNormalDagostino == 'True':	
# 					p_value_dagostino = float(line[12].replace('\n',''))
# 				else :
# 					p_value_dagostino = 1 / (zscore*zscore)
				

# 				file_write.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (my_node,sample_size,min_deg,max_deg,mean_deg,med_deg,sd,zscore,deg_covid,isNormalShapiro,p_value_shapiro,isNormalDagostino,p_value_dagostino))



# for filename in os.listdir("A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\Further_Analysis\\full_analysis_for_source_to_proteins_undirected\\pvalues_from_zscores\\outputFiles\\rel"):
# 	print(filename)
# 	with open("A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\Further_Analysis\\full_analysis_for_source_to_proteins_undirected\\pvalues_from_zscores\\outputFiles\\%s" % filename, 'r') as file_read :
		
# 		with open("A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\Further_Analysis\\full_analysis_for_source_to_proteins_undirected\\pvalues_from_zscores\\outputFiles\\rel\\aggregated_nodes_rel.tsv", 'a+') as file_write :
# 			data = file_read.readlines()
# 			# file_write.write(data[0])
# 			for line in data[1:] :
# 				line = line.split('\t')
# 				my_node = line[0]
# 				sample_size = float(line[1])	
# 				min_deg = float(line[2])
# 				max_deg = float(line[3])
# 				mean_deg = float(line[4])
# 				med_deg = float(line[5])	
# 				sd = float(line[6])	
# 				zscore = float(line[7])
# 				deg_covid = float(line[8])	
				
# 				isNormalShapiro	= line[9]
# 				# print(isNormalShapiro)
# 				if isNormalShapiro == 'True' :
# 					p_value_shapiro = float(line[10])
# 				else : 
# 					p_value_shapiro = 1 / (zscore*zscore)	
				
# 				isNormalDagostino = line[11]
# 				if isNormalDagostino == 'True':	
# 					p_value_dagostino = float(line[12].replace('\n',''))
# 				else :
# 					p_value_dagostino = 1 / (zscore*zscore)
				

# 				file_write.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (my_node,sample_size,min_deg,max_deg,mean_deg,med_deg,sd,zscore,deg_covid,isNormalShapiro,p_value_shapiro,isNormalDagostino,p_value_dagostino))



### add the node type

for filename in os.listdir("A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\Further_Analysis\\full_analysis_for_source_to_proteins_undirected\\pvalues_from_zscores\\outputFiles\\abs"):
	print(filename)
	with open("A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\Further_Analysis\\full_analysis_for_source_to_proteins_undirected\\pvalues_from_zscores\\outputFiles\\abs\\%s" % filename, 'r') as file_read :
		node_type = filename[:-7]
		print(node_type)
		if node_type == 'prot':
			my_node_type = 'Human PPI (target)'
		else :
			my_node_type = node_type
		with open("A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\Further_Analysis\\full_analysis_for_source_to_proteins_undirected\\pvalues_from_zscores\\outputFiles\\abs\\aggregated_nodes_abs_with_node_type.tsv", 'a+') as file_write :
			data = file_read.readlines()
			# file_write.write(data[0])
			for line in data[1:] :
				line = line.split('\t')
				my_node = line[0]
				sample_size = float(line[1])	
				min_deg = float(line[2])
				max_deg = float(line[3])
				mean_deg = float(line[4])
				med_deg = float(line[5])	
				sd = float(line[6])	
				zscore = float(line[7])
				deg_covid = float(line[8])	
				
				isNormalShapiro	= line[9]
				# print(isNormalShapiro)
				if isNormalShapiro == 'True' :
					p_value_shapiro = float(line[10])
				else : 
					p_value_shapiro = 1 / (zscore*zscore)	
				
				isNormalDagostino = line[11]
				if isNormalDagostino == 'True':	
					p_value_dagostino = float(line[12].replace('\n',''))
				else :
					p_value_dagostino = 1 / (zscore*zscore)
				

				file_write.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (my_node,my_node_type,sample_size,min_deg,max_deg,mean_deg,med_deg,sd,zscore,deg_covid,isNormalShapiro,p_value_shapiro,isNormalDagostino,p_value_dagostino))


# for filename in os.listdir("A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\Further_Analysis\\full_analysis_for_source_to_proteins_undirected\\pvalues_from_zscores\\outputFiles\\rel"):
# 	print(filename)
# 	with open("A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\Further_Analysis\\full_analysis_for_source_to_proteins_undirected\\pvalues_from_zscores\\outputFiles\\rel\\%s" % filename, 'r') as file_read :
# 		node_type = filename[:-7]
# 		print(node_type)
# 		if node_type == 'prot':
# 			my_node_type = 'Human PPI (target)'
# 		else :
# 			my_node_type = node_type
# 		with open("A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\Further_Analysis\\full_analysis_for_source_to_proteins_undirected\\pvalues_from_zscores\\outputFiles\\rel\\aggregated_nodes_rel_with_node_type.tsv", 'a+') as file_write :
# 			data = file_read.readlines()
# 			# file_write.write(data[0])
# 			for line in data[1:] :
# 				line = line.split('\t')
# 				my_node = line[0]
# 				sample_size = float(line[1])	
# 				min_deg = float(line[2])
# 				max_deg = float(line[3])
# 				mean_deg = float(line[4])
# 				med_deg = float(line[5])	
# 				sd = float(line[6])	
# 				zscore = float(line[7])
# 				deg_covid = float(line[8])	
				
# 				isNormalShapiro	= line[9]
# 				# print(isNormalShapiro)
# 				if isNormalShapiro == 'True' :
# 					p_value_shapiro = float(line[10])
# 				else : 
# 					p_value_shapiro = 1 / (zscore*zscore)	
				
# 				isNormalDagostino = line[11]
# 				if isNormalDagostino == 'True':	
# 					p_value_dagostino = float(line[12].replace('\n',''))
# 				else :
# 					p_value_dagostino = 1 / (zscore*zscore)
				

# 				file_write.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (my_node,my_node_type,sample_size,min_deg,max_deg,mean_deg,med_deg,sd,zscore,deg_covid,isNormalShapiro,p_value_shapiro,isNormalDagostino,p_value_dagostino))
