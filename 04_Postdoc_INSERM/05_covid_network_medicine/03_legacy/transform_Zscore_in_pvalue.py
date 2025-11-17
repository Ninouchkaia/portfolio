#!/usr/bin/python
# -*- coding: utf-8 -*-
# jupyter notebook --browser='C:/Program Files (x86)/Google/Chrome/Application/chrome.exe %s'

isNormalAbsDict = {}
with open("pvalues_from_zscores\\KS_test_alpha0.05_covid_GO_degree_to_proteins_in_mock_networks.tsv" , 'r') as file_read :
	data = file_read.readlines()
	for line in data[1:] :
		line = line.split("\t")
		GO = line[0]
		isNormal = line[1]
		isNormalAbsDict[GO] = isNormal

isNormalRelDict = {}
with open("pvalues_from_zscores\\KS_test_alpha0.05_covid_GO_structural_strength_to_proteins_in_mock_networks.tsv" , 'r') as file_read :
	data = file_read.readlines()
	for line in data[1:] :
		line = line.split("\t")
		GO = line[0]
		isNormal = line[1]
		isNormalRelDict[GO] = isNormal

with open("pvalues_from_zscores\\alldataGO_normal.tsv" , 'a+') as file_write :
	file_write.write("covidGO\tcovid-degree\tZ-score\tcovid-degree-rel\tZ-score-rel\tisNormalAbs\tisNormalRel\n")
	with open("pvalues_from_zscores\\alldataGO.tsv" , 'r') as file_read :
		data = file_read.readlines()
		for line in data[1:] :
			line = line.split("\t")
			GO = line[0]
			covidDeg = line[1]
			Zscore = line[2]
			covidDegRel = line[3]
			ZscoreRel = line[4].replace("\n", "")
			file_write.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (GO, covidDeg, Zscore, covidDegRel, ZscoreRel, isNormalAbsDict[GO], isNormalRelDict[GO]))