#!/usr/bin/python
# -*- coding: utf-8 -*-
# jupyter notebook --browser='C:/Program Files (x86)/Google/Chrome/Application/chrome.exe %s'


import timeit
start = timeit.default_timer()

import sys
import os
import csv
from operator import itemgetter
import networkx as nx
from networkx.algorithms import community
from collections import Counter
from collections import defaultdict
from collections import OrderedDict







#####################################################################################################
#########################  Compute network source to protein degree distribution ####################
##########################################  based on edges to protein ###############################
#####################################################################################################


# nodes_path = 'A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid-2\\Virus_host_interactomes_thresh25\\thresh=0.25\\covid19\\nodes.csv'
# edges_path = 'A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid-2\\Virus_host_interactomes_thresh25\\thresh=0.25\\covid19\\edges.csv'

# nodes_path = 'C:\\Downloads\\Projects\\covid_networks\\Virus_host_interactomes_thresh25\\thresh=0.25\\covid19\\nodes.csv'
# edges_path = 'C:\\Downloads\\Projects\\covid_networks\\Virus_host_interactomes_thresh25\\thresh=0.25\\covid19\\edges.csv'

# rootdir = 'C:\\Downloads\\Projects\\covid_networks\\Virus_host_interactomes_thresh25\\thresh=0.25'
rootdir = 'A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid-2\\Virus_host_interactomes_thresh25\\thresh=0.25'

# for subdir, dirs, files in os.walk(rootdir):
#     for file in files:
#         print(os.path.join(subdir, file))

viruses_folders_list = []
for subdir, dirs, files in os.walk(rootdir):
	for name in dirs:
	    viruses_folders_list.append(os.path.join(rootdir, name))
# print(viruses_folders_list)

covid_interactors = ["WFS1","QSOX2","FAM162A","EXOSC3","SRP54","SLC1A3","GIGYF2","PUSL1","TK2","TAPT1","MFGE8","TOMM70A","NOL10","NUP98","TARS2","PPIL3","USP54","PPT1","POR","HDAC2","CEP135","RAB5C","SAAL1","POLA2","VPS11","SRP19","ACSL3","HSBP1","DCAF7","SDHB","TCTN3","MRPS2","GCC2","NUP62","CENPF","RTN4","ZYG11B","MPHOSPH10","PRIM1","NGLY1","ARF6","RRP9","RAB7A","RDX","ATP5I","BCS1L","GOLGA2","SLC25A21","HECTD1","TCF12","RAE1","ECSIT","PRIM2","MDN1","LOX","FKBP15","ERP44","FAM134C","CWC27","ATP5L","REEP6","RAB2A","UGGT2","MTCH1","YIF1A","INHBE","SLC30A9","NINL","GCC1","SLC30A6","ETFA","NUP88","PDE4DIP","TUBGCP3","SEPSECS","PCM1","DCTPP1","NDUFB9","SCARB1","NUP210","RHOT2","RAB8A","ZNF318","ITGB1","CUL2","RAB18","GORASP1","AATF","SIL1","MAT2B","EXOSC5","UPF1","KDELC2","TRMT1","SIGMAR1","CHPF2","NUTF2","GFER","INTS4","GRPEL1","NUP54","ABCC1","G3BP1","STOM","CISD3","ERLEC1","SMOC1","DDX10","AKAP8L","STOML2","ERGIC1","RHOA","SLU7","HS6ST2","CEP350","TBK1","NGDN","NUPL1","FBN1","PIGO","GOLGB1","ALG11","PKP2","NLRX1","AASS","TIMM10B","TRIM59","AAR2","TOR1AIP1","DNMT1","HEATR3","AGPS","PLEKHF2","ANO6","COLGALT1","PMPCA","HMOX1","SLC30A7","ACADM","TUBGCP2","ARL6IP6","TM2D3","MIB1","SLC9A3R1","FBXL12","VIMP","CSDE1","SLC44A2","DNAJC11","CYB5B","ACAD9","CIT","EXOSC8","SPG20","CHCHD1","USP13","CEP68","CLIP4","TIMM8B","GLA","TMEM39B","PTGES2","EDEM3","CCDC86","FAM98A","NARS2","DPH5","MYCBP2","ALG8","GGCX","CEP112","OS9","NAT14","G3BP2","ATP6V1A","BAG5","RAB1A","GGH","RAB14","JAKMIP1","PSMD8","PRKACA","ADAMTS1","NIN","LMAN2","ALG5","BRD4","F2RL1","SCCPDH","COL6A1","PRKAR2B","RIPK1","ATP1B1","FKBP7","GOLGA3","CHMP2A","HS2ST1","SRP72","DNAJC19","ATP6AP1","MIPOL1","SIRT5","NPTX1","REEP5","ZC3H7A","AP2M1","TIMM9","AKAP8","UBXN8","MARC1","VPS39","TMED5","FASTKD5","FBN2","NEK9","PITRM1","FBLN5","THTPA","GNG5","LARP4B","AP3B1","PRKAR2A","PCNT","LARP7","LARP1","TBCA","MARK3","AKAP9","MARK2","DPY19L1","DCAKD","UBAP2","CSNK2A2","ADCK4","GNB1","SBNO1","CHPF","EMC1","MEPCE","MAP7D1","STC2","GPX1","CRTC3","GHITM","ZDHHC5","PABPC1","SLC6A15","TLE1","GOLGA7","WHSC1","NPC2","OCIAD1","CEP250","HOOK1","AES","SLC27A2","CDK5RAP2","FAR2","IMPDH2","RAP1GDS1","AP2A2","ERC1","FGFR1","ZC3H18","CYB5R3","MRPS25","GRIPAP1","ZNF503","TIMM10","ERO1LB","RAB10","PTBP2","BCKDK","ATP13A3","SUN2","IDE","ATE1","TLE3","NDFIP2","RBM41","PLAT","NDUFAF2","POFUT1","NDUFAF1","CLCC1","GPAA1","PMPCB","SNIP1","FOXRED2","GTF2F2","NUP214","SCAP","FAM8A1","MOV10","EIF4E2","UBAP2L","RNF41","NEU1","ERMP1","PABPC4","BRD2","POLA1","CNTRL","MOGS","KIAA1033","TBKBP1","KIAA0368","COMT","CSNK2B","IL17RA","TMEM97","PCSK6","SDF2","MARK1","BZW2","TYSND1","RALA","PVR","KDELC1","TOR1A","GDF15","PLOD2","PLD3","PRRC2B","PLEKHA5","LETMD1","FYCO1","PDZD11","HYOU1","EIF4H"]
viruses_with_shared_interactors = []
for path in viruses_folders_list :
	virus = path.split("\\")[-1]
	print(virus)
	nodes_path = '%s\\nodes.csv' % path
	edges_path = '%s\\edges.csv' % path
	# make viral_proteins dict
	descr_dict_nodes, sapiens_dict, virus_dict = {},{},{}
	with open(nodes_path, 'r') as nodecsv: # Open the file                       
		nodereader = csv.reader(nodecsv) # Read the csv  
		for line in nodereader :
			node_id = line[0]
			node_name = line[1]
			node_type = line[2].replace('\n', '')
			descr_dict_nodes[node_id] = node_name
			if (node_type == '1.0') :
				sapiens_dict[node_id] = node_name
			elif (node_type == '0.0') :
				virus_dict[node_id] = node_name
			else :
				print(line, "Warning : unidentified nodes !\n")

	with open(edges_path, 'r') as edgecsv: # Open the file
		edgereader = csv.reader(edgecsv) # Read the csv     
		edges = [tuple(e) for e in edgereader] # Retrieve the data
		# print(edges)

	G_interactome = nx.Graph()
	G_interactome.add_nodes_from(descr_dict_nodes.keys())
	G_interactome.add_edges_from(edges)

	print(nx.info(G_interactome))
	# Number of nodes: 19972
	# Number of edges: 737999
	# Average degree:  73.9034

	#### fetch viral direct interactors (1st order)
	# direct_interactors1 = []
	# direct_interactors2 = []
	
	for u,v,c in G_interactome.edges(data=True) :
		if (u in covid_interactors and v in virus_dict) :
			print(virus, u, v)
			viruses_with_shared_interactors.append(virus)
		elif (v in covid_interactors and u in virus_dict) :
			print(virus, u, v)
			direct_interactors.append(virus)

viruses_with_shared_interactors = list(set(viruses_with_shared_interactors))
print(viruses_with_shared_interactors)









stop = timeit.default_timer()
print(stop - start)  