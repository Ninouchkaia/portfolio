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


# ### adjust p-values > 1 to 1

with open("A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\Further_Analysis\\full_analysis_for_source_to_proteins_undirected\\pvalues_from_zscores\\outputFiles\\final processing\\bootstrap_structural_degrees.tsv", 'r') as file_read :
	with open("A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\Further_Analysis\\full_analysis_for_source_to_proteins_undirected\\pvalues_from_zscores\\outputFiles\\final processing\\bootstrap_structural_degrees_adjusted.tsv", 'a+') as file_write :
		data = file_read.readlines()
		file_write.write(data[0])
		for line in data[1:] :
			line = line.split('\t')
			my_node = line[0]
			my_node_type = line[1]
			sample_size = float(line[2])	
			min_deg = float(line[3])
			max_deg = float(line[4])
			mean_deg = float(line[5])
			med_deg = float(line[6])	
			sd = float(line[7])	
			zscore = float(line[8])
			deg_covid = float(line[9])	
			isNormalShapiro	= line[10]
			p_value_shapiro = float(line[11])
			isNormalDagostino = line[12]
			p_value_dagostino = float(line[13].replace('\n',''))

			if sd == 0.0 :
				zscore = 0.0
				if isNormalShapiro == 'True' :
					p_value_shapiro = 0.5
				else : 
					p_value_shapiro = 1
			
				if isNormalDagostino == 'True':	
					p_value_dagostino = 0.5
				else :
					p_value_dagostino = 1

			if p_value_shapiro > 1 :
				p_value_shapiro = 1
			if p_value_dagostino > 1 :
				p_value_dagostino = 1

			file_write.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (my_node,my_node_type,sample_size,min_deg,max_deg,mean_deg,med_deg,sd,zscore,deg_covid,isNormalShapiro,p_value_shapiro,isNormalDagostino,p_value_dagostino))




with open("A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\Further_Analysis\\full_analysis_for_source_to_proteins_undirected\\pvalues_from_zscores\\outputFiles\\final processing\\bootstrap_structural_strength.tsv", 'r') as file_read :
	with open("A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\Further_Analysis\\full_analysis_for_source_to_proteins_undirected\\pvalues_from_zscores\\outputFiles\\final processing\\bootstrap_structural_strength_adjusted.tsv", 'a+') as file_write :
		data = file_read.readlines()
		file_write.write(data[0])
		for line in data[1:] :
			print(line)
			line = line.split('\t')
			my_node = line[0]
			my_node_type = line[1]
			sample_size = float(line[2])	
			min_deg = float(line[3].replace(',','.'))
			max_deg = float(line[4])
			mean_deg = float(line[5])
			med_deg = float(line[6])	
			sd = float(line[7])	
			zscore = float(line[8])
			deg_covid = float(line[9])	
			isNormalShapiro	= line[10]
			p_value_shapiro = float(line[11])
			isNormalDagostino = line[12]
			p_value_dagostino = float(line[13].replace('\n',''))

			if sd == 0.0 :
				zscore = 0.0
				if isNormalShapiro == 'True' :
					p_value_shapiro = 0.5
				else : 
					p_value_shapiro = 1
			
				if isNormalDagostino == 'True':	
					p_value_dagostino = 0.5
				else :
					p_value_dagostino = 1
			if p_value_shapiro > 1 :
				p_value_shapiro = 1
			if p_value_dagostino > 1 :
				p_value_dagostino = 1

			file_write.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (my_node,my_node_type,sample_size,min_deg,max_deg,mean_deg,med_deg,sd,zscore,deg_covid,isNormalShapiro,p_value_shapiro,isNormalDagostino,p_value_dagostino))
