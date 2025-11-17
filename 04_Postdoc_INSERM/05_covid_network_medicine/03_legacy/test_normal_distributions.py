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


## parse covidnet degree data :
# files = ["COVID19_GDDS_Disease_degrees_to_proteins.tsv", "COVID19_GDDS_GO_degrees_to_proteins.tsv", "COVID19_GDDS_Symptom_degrees_to_proteins.tsv"]
files = ["COVID19_GDDS_Human PPI (target)_degrees_to_proteins.tsv", "COVID19_GDDS_Drug_degrees_to_proteins.tsv"]

deg_in_covidnet_dict = {}
for file in files :
	print(file)
	with open(file, 'r') as file_read :
		data = file_read.readlines()
		for line in data :
			line = line.split('\t')
			my_node = line[0]
			absdeg = int(line[1])
			reldeg = float(line[2].replace('\n',''))
			deg_in_covidnet_dict[my_node] = (absdeg, reldeg)

print(len(deg_in_covidnet_dict)) #9820 (4176 + 3487 + 2157) OK
		

# ## parse the distributions
# files_absolute = [["covid_disease_degree_to_proteins_distribution_in_mock_networks.tsv", "COVID19_GDDS_Results_Disease_to_CovidProteins-withZeros.tsv"], 
# ["covid_GO_degree_to_proteins_distribution_in_mock_networks.tsv", "COVID19_GDDS_Results_GO_to_CovidProteins-withZeros.tsv"],
# ["covid_symptom_degree_to_proteins_distribution_in_mock_networks.tsv", "COVID19_GDDS_Results_Symptom_to_CovidProteins-withZeros.tsv"]]

#diseases
# files_absolute = [["covid_disease_degree_to_proteins_distribution_in_mock_networks.tsv", "COVID19_GDDS_Results_Disease_to_CovidProteins-withZeros.tsv"]]

#GO
# files_absolute = [["covid_GO_degree_to_proteins_distribution_in_mock_networks.tsv", "COVID19_GDDS_Results_GO_to_CovidProteins-withZeros.tsv"]]

#Symptom
# files_absolute = [["covid_symptom_degree_to_proteins_distribution_in_mock_networks.tsv", "COVID19_GDDS_Results_Symptom_to_CovidProteins-withZeros.tsv"]]

## protein
# files_absolute = [["covid_proteins_structural_degree_distribution_in_mock_networks.tsv", "COVID19_GDDS_Results_Human PPI (target)_to_CovidProteins-withZeros.tsv"]]

##drug
files_absolute = [["covid_drug_degree_to_proteins_distribution_in_mock_networks.tsv", "COVID19_GDDS_Results_Drug_to_CovidProteins-withZeros.tsv"]]


# files_relative = [["covid_disease_structural_strength_to_proteins_distribution_in_mock_networks.tsv", "COVID19_GDDS_Results_Covid_Disease_to_Proteins_Relative-withZeros.tsv"],
# ["covid_GO_structural_strength_to_proteins_distribution_in_mock_networks.tsv", "COVID19_GDDS_Results_Covid_GO_to_Proteins_Relative-withZeros.tsv"],
# ["covid_symptom_structural_strength_to_proteins_distribution_in_mock_networks.tsv", "COVID19_GDDS_Results_Covid_Symptom_to_Proteins_Relative-withZeros.tsv"]]

#protein
files_relative = [["covid_drug_structural_strength_to_proteins_distribution_in_mock_networks.tsv", "COVID19_GDDS_Results_Covid_Drug_to_Proteins_Relative-withZeros.tsv"]]



### for absolute degrees using Shapiro Normality test
for file_pair in files_absolute :
	print(file_pair)
	## fetch distributions from file
	with open(file_pair[0], 'r') as file_read :
		distributions = {} #build new dict when new node type (ie new file_pair)
		data = file_read.readlines()
		for line in data[1:] :
			line = line.split('\t')
			my_node = line[0]
			degree_list = [ int(deg) for deg in line[1:]]
			distributions[my_node] = np.array(degree_list)
		print(len(distributions))
	### write output file with all data
	with open('%s_NormTest_alpha_0.05.tsv' % file_pair[1][:-4], 'a+') as file_write :
		file_write.write("%s\t# of nets it is found in\tmin-degree\tmax-degree\tmean-degree\tmedian-degree\tsd-degree\tz-score\tdeg_covid\tisNormalShapiro\tp-value-shapiro\tisNormalDagostino\tp-value-dagostino\n" % file_pair[1].split('_')[3])
		### fetch results on degrees and zscores already computed
		with open(file_pair[1], 'r') as file_read :
			data = file_read.readlines()
			
			### parse data from results
			for line in data[1:] :
				line = line.split('\t')
				my_node = line[0]
				sample_size = float(line[1])
				min_deg = float(line[2])
				max_deg = float(line[3])
				mean_deg = float(line[4])
				med_deg = float(line[5])
				sd = float(line[6])
				zscore = float(line[8].replace('\n',''))
				deg_covid = float(deg_in_covidnet_dict[my_node][0]) #absdeg !
				if my_node not in distributions :
					distributions[my_node] = []
				distribution = distributions[my_node]
				
				### run shapiro test on stored distribution
				if(len(distributions[my_node])>3):
					stat, p = shapiro(distribution)
					alpha = 0.05
					if p < alpha :
						isNormalShapiro = False
					else :
						isNormalShapiro = True

					### calculate p-value
					if isNormalShapiro :
						if math.isnan(zscore) :
							p_value_shapiro = 0.5
						else :
							p_value_shapiro = 1 - ( erf(abs(zscore) / math.sqrt(2) ) )
					else :
						if (deg_covid == mean_deg) : ## zscore = 0, then p-value = 0.5
							p_value_shapiro = 0.5
						else :	
							p_value_shapiro = (sd * sd) / (sample_size * abs(deg_covid - mean_deg))	
				else :
					p_value_shapiro = float('NaN')

				### run D’Agostino’s K^2 test on stored distribution
				if(len(distributions[my_node])>8):
					stat, p = normaltest(distribution)
					alpha = 0.05
					if p < alpha :
						isNormalDagostino = False
					else :
						isNormalDagostino = True

					### calculate p-value
					if isNormalDagostino :
						if math.isnan(zscore) :
							p_value_shapiro = 0.5
						else :
							p_value_dagostino = 1 - ( erf(abs(zscore) / math.sqrt(2) ) )
					else :
						if (deg_covid == mean_deg) : ## zscore = 0, then p-value = 0.5
							p_value_dagostino = 0.5
						else :	
							p_value_dagostino = (sd * sd) / (sample_size * abs(deg_covid - mean_deg))	
				else :
					p_value_dagostino = float('NaN')
				### write the output
				file_write.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (my_node,sample_size,min_deg,max_deg,mean_deg,med_deg,sd,zscore,deg_covid,isNormalShapiro,p_value_shapiro,isNormalDagostino,p_value_dagostino))

### for relative degrees using Shapiro Normality test
for file_pair in files_relative :
	print(file_pair)
	## fetch distributions from file
	with open(file_pair[0], 'r') as file_read :
		distributions = {} #build new dict when new node type (ie new file_pair)
		data = file_read.readlines()
		for line in data[1:] :
			line = line.split('\t')
			my_node = line[0]
			degree_list = [ float(deg) for deg in line[1:]]
			distributions[my_node] = np.array(degree_list)
		print(len(distributions))
	### write output file with all data
	with open('%s_NormTest_alpha_0.05.tsv' % file_pair[1][:-4], 'a+') as file_write :
		file_write.write("%s\t# of nets it is found in\tmin-degree-rel\tmax-degree-rel\tmean-degree-rel\tmedian-degree-rel\tsd-degree-rel\tz-score-rel\tdeg_covid-rel\tisNormalShapiro\tp-value-shapiro\tisNormalDagostino\tp-value-dagostino\n" % file_pair[1].split('_')[3])
		### fetch results on degrees and zscores already computed
		with open(file_pair[1], 'r') as file_read :
			data = file_read.readlines()
			
			### parse data from results
			for line in data[1:] :
				line = line.split('\t')
				my_node = line[0]
				sample_size = float(line[1])
				min_deg_rel = float(line[2])
				max_deg_rel = float(line[3])
				mean_deg_rel = float(line[4])
				med_deg_rel = float(line[5])
				sd_rel = float(line[6])
				zscore_rel = float(line[8].replace('\n',''))
				deg_covid_rel = float(deg_in_covidnet_dict[my_node][1]) #reldeg !

				if my_node not in distributions :
					distributions[my_node] = []
				distribution = distributions[my_node]
				
				print(my_node)
				### run shapiro test on stored distribution
				if(len(distributions[my_node])>3):
					stat, p = shapiro(distribution)
					alpha = 0.05
					if p < alpha :
						isNormalShapiro = False
					else :
						isNormalShapiro = True

					### calculate p-value
					if isNormalShapiro :
						if math.isnan(zscore_rel) :
							p_value_shapiro = 0.5
						else :
							p_value_shapiro = 1 - ( erf(abs(zscore_rel) / math.sqrt(2) ) )
					else :
						if (deg_covid_rel == mean_deg_rel) : ## zscore_rel = 0, then p-value = 0.5
							p_value_shapiro = 0.5
						else :	
							p_value_shapiro = (sd_rel * sd_rel) / (sample_size * abs(deg_covid_rel - mean_deg_rel))	
				else :
					p_value_shapiro = float('NaN')

				### run D’Agostino’s K^2 test on stored distribution
				if(len(distributions[my_node])>8):
					stat, p = normaltest(distribution)
					alpha = 0.05
					if p < alpha :
						isNormalDagostino = False
					else :
						isNormalDagostino = True

					### calculate p-value
					if isNormalDagostino :
						if math.isnan(zscore_rel) :
							p_value_shapiro = 0.5
						else :
							p_value_dagostino = 1 - ( erf(abs(zscore_rel) / math.sqrt(2) ) )
					else :
						if (deg_covid_rel == mean_deg_rel) : ## zscore_rel = 0, then p-value = 0.5
							p_value_dagostino = 0.5
						else :	
							p_value_dagostino = (sd_rel * sd_rel) / (sample_size * abs(deg_covid_rel - mean_deg_rel))	
				else :
					p_value_dagostino = float('NaN')
				### write the output
				file_write.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (my_node,sample_size,min_deg_rel,max_deg_rel,mean_deg_rel,med_deg_rel,sd_rel,zscore_rel,deg_covid_rel,isNormalShapiro,p_value_shapiro,isNormalDagostino,p_value_dagostino))




stop = timeit.default_timer()
print(stop - start)  









































### for absolute degrees
# with open('pvalues_from_zscores\\covid_GO_degree_to_proteins_distribution_in_mock_networks.tsv', 'r') as file_read :
# 	data = file_read.readlines()
# 	for line in data[1:] :
# 		line = line.split('\t')
# 		GO = line[0]
# 		degree_list = [ int(deg) for deg in line[1:]]
# 		distributions[GO] = np.array(degree_list)

# with open('pvalues_from_zscores\\test_normal_covid_GO_degree_to_proteins_in_mock_networks.tsv', 'a+') as file_write :
# 	file_write.write("covid_GO\tisNormalDistribution\tpValueForNormalDistribution\n")
# 	for GO in distributions :
# 		distribution = distributions[GO]
# 		k2, p = stats.normaltest(distribution)
# 		alpha = 1e-3
# 		if p < alpha :
# 			isNormal = False
# 		else :
# 			isNormal = True
# 		file_write.write("%s\t%s\t%s\n" % (GO, isNormal, p))

## for relative degrees
# with open('pvalues_from_zscores\\covid_GO_structural_strength_to_proteins_distribution_in_mock_networks.tsv', 'r') as file_read :
# 	data = file_read.readlines()
# 	for line in data[1:] :
# 		line = line.split('\t')
# 		GO = line[0]
# 		degree_list = [ float(deg) for deg in line[1:]]
# 		distributions[GO] = np.array(degree_list)

# with open('pvalues_from_zscores\\test_normal_covid_GO_structural_strength_to_proteins_in_mock_networks.tsv', 'a+') as file_write :
# 	file_write.write("covid_GO\tisNormalDistribution\tpValueForNormalDistribution\n")
# 	for GO in distributions :
# 		distribution = distributions[GO]
# 		k2, p = stats.normaltest(distribution)
# 		alpha = 1e-3
# 		if p < alpha :
# 			isNormal = False
# 		else :
# 			isNormal = True
# 		file_write.write("%s\t%s\t%s\n" % (GO, isNormal, p))


# ### for absolute degrees using Kolmogorov Smirnov Statistic test
# with open('pvalues_from_zscores\\covid_GO_degree_to_proteins_distribution_in_mock_networks.tsv', 'r') as file_read :
# 	data = file_read.readlines()
# 	for line in data[1:] :
# 		line = line.split('\t')
# 		GO = line[0]
# 		degree_list = [ int(deg) for deg in line[1:]]
# 		distributions[GO] = np.array(degree_list)

# from scipy.stats import kstest, norm
# with open('pvalues_from_zscores\\KS_test_alpha0.05_covid_GO_degree_to_proteins_in_mock_networks.tsv', 'a+') as file_write :
# 	file_write.write("covid_GO\tisNormalDistribution\tpValueForNormalDistribution\n")
# 	for GO in distributions :
# 		distribution = distributions[GO]
# 		ks_statistic, p = kstest(distribution, 'norm')
# 		alpha = 5e-3
# 		if p == 0 :
# 			isNormal = True
# 		elif p < alpha :
# 			isNormal = False
# 		elif p > alpha :
# 			isNormal = True
# 		file_write.write("%s\t%s\t%s\n" % (GO, isNormal, p))


# ## for relative degrees using Kolmogorov Smirnov Statistic test
# with open('pvalues_from_zscores\\covid_GO_structural_strength_to_proteins_distribution_in_mock_networks.tsv', 'r') as file_read :
# 	data = file_read.readlines()
# 	for line in data[1:] :
# 		line = line.split('\t')
# 		GO = line[0]
# 		degree_list = [ float(deg) for deg in line[1:]]
# 		distributions[GO] = np.array(degree_list)

# from scipy.stats import kstest, norm
# with open('pvalues_from_zscores\\KS_test_alpha0.05_covid_GO_structural_strength_to_proteins_in_mock_networks.tsv', 'a+') as file_write :
# 	file_write.write("covid_GO\tisNormalDistribution\tpValueForNormalDistribution\n")
# 	for GO in distributions :
# 		distribution = distributions[GO]
# 		ks_statistic, p = kstest(distribution, 'norm')
# 		alpha = 5e-3
# 		if p == 0 :
# 			isNormal = True
# 		elif p < alpha :
# 			isNormal = False
# 		elif p > alpha :
# 			isNormal = True
# 		file_write.write("%s\t%s\t%s\n" % (GO, isNormal, p))






# import math
# import scipy

# p = scipy.stats.norm.sf(2)*2
# q = 1 - scipy.special.erf(2/math.sqrt(2))

# print(p,q)

























# distributions_aggregated = []
# from scipy.stats import shapiro
# ### for absolute degrees using Shapiro Wilk Statistic test
# with open('pvalues_from_zscores\\covid_GO_degree_to_proteins_distribution_in_mock_networks.tsv', 'r') as file_read :
# 	data = file_read.readlines()
# 	for line in data[1:] :
# 		line = line.split('\t')
# 		GO = line[0]
# 		degree_list = [ int(deg) for deg in line[1:]]
# 		# distributions[GO] = np.array(degree_list)
# 		distributions_aggregated.append(degree_list)

# # print(distributions_aggregated)
# distributions_aggregated = [j for i in distributions_aggregated for j in i]
# # print(distributions_aggregated)

# from scipy.stats import kstest, norm
# with open('pvalues_from_zscores\\test2.tsv', 'a+') as file_write :
# 	file_write.write("covid_GO\tisNormalDistribution\tpValueForNormalDistribution\n")
# 	ks_statistic, p = kstest(distributions_aggregated[1000:1250], 'norm')
# 	alpha = 1e-3
# 	if p == 0 :
# 		isNormal = True
# 	elif p < alpha :
# 		isNormal = False
# 	elif p > alpha :
# 		isNormal = True
# 	file_write.write("%s\t%s\t%s\n" % (GO, isNormal, p))

# import math
# import scipy

# # p=scipy.stats.norm.sf(0.1166)*2
# p = 1 - erf(10.1166/math.sqrt(2))
# # p = 1 - ndtr(-0.1166)

# print(p)




# import matplotlib.pyplot as plt
# import numpy as np

# plt.hist(distributions_aggregated, bins=100)
# plt.show()







# from scipy.stats import kstest, norm
# my_data = norm.rvs(size=1000)
# print(type(my_data))
# ks_statistic, p_value = kstest(my_data, 'norm')
# print(ks_statistic, p_value)







# print("p = {:g}".format(p))

# if p < alpha:  # null hypothesis: x comes from a normal distribution
#     # print("The null hypothesis can be rejected")
#     print("NOT normal")
# else:
#     print("normal")
