check_random_nets_stats.py,

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
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
from collections import Counter
from collections import defaultdict
from collections import OrderedDict
import statistics 
import re
import gzip



# nodefiles_numbers = []
# edgefiles_numbers = []
# for filename in os.listdir("A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\random_networks\\Mock_networks_21Apr"):
# 	if 'nodes_' in filename :
# 		filename = filename.split(".csv")[0]
# 		nodefiles_numbers.append(int(filename[19:]))
# 	if 'edges_' in filename :
# 		filename = filename.split(".csv")[0]
# 		edgefiles_numbers.append(int(filename[19:]))

# # log_write.write("len(nodefiles_numbers) = %s\n" % len(nodefiles_numbers)) #2051
# # log_write.write("len(edgefiles_numbers) = %s\n" % len(edgefiles_numbers)) #2051
# # log_write.write("difference between edgefiles_numbers and nodefiles_numbers = %s\n" % len(set(nodefiles_numbers).difference(set(edgefiles_numbers)))) #0

# nodefiles_numbers = sorted(nodefiles_numbers)
# edgefiles_numbers = sorted(edgefiles_numbers)

# networks_to_avoid = [3, 11, 22, 23, 33, 34, 38, 57, 72, 77, 123, 127, 148, 156, 157, 183, 205, 209, 210, 217, 232, 252, 268, 281, 291, 312, 341, 353, 362, 364, 413, 421, 433, 443, 451, 494, 525, 532, 569, 592, 629, 658, 662, 674, 676, 701, 705, 712, 736, 757, 788, 796, 825, 827, 838, 873, 876, 908, 910, 916, 917, 918, 938, 958, 965, 1088, 1125, 1128, 1145, 1146, 1165, 1171, 1183, 1185, 1192, 1195, 1211, 1244, 1254, 1257, 1300, 1353, 1354, 1364, 1372, 1394, 1395, 1401, 1410, 1427, 1454, 1477, 1497, 1517, 1532, 1549, 1551, 1575, 1578, 1611, 1624, 1629, 1637, 1650, 1655, 1697, 1700, 1713, 1736, 1745, 1762, 1790, 1800, 1815, 1824, 1827, 1866, 1872, 1886, 1887, 2053, 2072, 2106, 2115, 2128, 2160, 2217, 2343, 2383, 2388, 2413, 2458, 2465]
# nodefiles_numbers = [item for item in nodefiles_numbers if item not in networks_to_avoid]

# protein_nodes_total = []
# ppi_total = []

# for iteration in nodefiles_numbers[1:] :
# 	print("\n########################## ITERATION %s  ###########################\n" % (iteration))
# 	# log_write.write("\n########################## ITERATION %s ###########################\n" % iteration)

# 	with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\random_networks\\Mock_networks_21Apr\\COVID19_GDDS_nodes_%s.csv' % str(iteration), 'r') as nodecsv: # Open the file                   
# 		nodereader = csv.reader(nodecsv) # Read the csv  
# 		nodes = [n for n in nodereader][1:]                     
# 		node_names = [n[0] for n in nodes] # Get a list of only the node names 

# 	with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\random_networks\\Mock_networks_21Apr\\COVID19_GDDS_edges_%s.csv' % str(iteration), 'r') as edgecsv: # Open the file
# 		edgereader = csv.reader(edgecsv) # Read the csv     
# 		edges = [tuple(e) for e in edgereader][1:] # Retrieve the data

# 		# log_write.write("len(node_names) : %s\n" % len(node_names)) #16007 --> 13448
# 		# log_write.write("len(edges) : %s\n" % len(edges)) #53917 --> 37835

# 	G_random = nx.Graph()
# 	G_random.add_nodes_from(node_names)
# 	G_random.add_edges_from(edges)

# 	descr_dict_random = {}
# 	description_random_set=set()
# 	for node in nodes: 
# 		descr_dict_random[node[0]] = node[2]
# 		description_random_set.add(node[2])

# 	protein_nodes = [node for node in descr_dict_random if descr_dict_random[node] == 'Human PPI (target)']
# 	protein_nodes_total.append(len(protein_nodes))

# 	# print(nx.info(G_random))
# 	# log_write.write("\nRANDOM GRAPH %s\n" % iteration)
# 	# log_write.write("%s\n" % nx.info(G_random))

# 	#####################################################################################################
# 	################################# GET EDGES CONNECTING sources TO proteins ################
# 	#####################################################################################################

# 	ppi = 0
# 	for u,v,c in G_random.edges(data=True) :	
# 		if (descr_dict_random[u] == 'Human PPI (target)' and descr_dict_random[v] == 'Human PPI (target)') :
# 			ppi = ppi + 1
	
# 	ppi_total.append(ppi)

# print("protein_nodes_total", protein_nodes_total)
# print("ppi_total", ppi_total)


# ppi_total = [354, 367, 389, 420, 350, 343, 384, 343, 379, 467, 295, 363, 405, 369, 347, 457, 319, 275, 365, 327, 370, 437, 307, 428, 361, 419, 414, 342, 492, 229, 423, 337, 312, 460, 442, 407, 416, 414, 368, 372, 366, 497, 300, 431, 352, 426, 348, 367, 527, 392, 427, 453, 336, 445, 391, 377, 373, 695, 323, 417, 465, 348, 337, 383, 437, 336, 426, 347, 354, 419, 347, 399, 431, 426, 357, 422, 437, 367, 296, 339, 327, 358, 317, 438, 365, 402, 383, 492, 332, 353, 445, 364, 294, 454, 426, 448, 416, 309, 393, 360, 405, 444, 466, 489, 515, 351, 307, 265, 405, 304, 287, 254, 308, 350, 338, 336, 270, 330, 403, 310, 385, 308, 418, 335, 387, 483, 431, 379, 371, 425, 381, 392, 417, 373, 346, 315, 391, 281, 389, 252, 461, 361, 371, 356, 293, 386, 404, 390, 381, 331, 321, 296, 393, 336, 369, 345, 477, 514, 421, 290, 454, 339, 401, 338, 361, 307, 459, 389, 309, 404, 288, 476, 347, 300, 374, 339, 333, 417, 402, 340, 386, 291, 349, 416, 325, 390, 317, 496, 290, 440, 385, 413, 380, 319, 442, 260, 353, 364, 376, 352, 314, 447, 348, 379, 385, 435, 410, 321, 381, 419, 333, 467, 442, 444, 296, 383, 424, 453, 356, 407, 358, 353, 337, 361, 349, 370, 332, 364, 450, 422, 393, 353, 366, 457, 415, 400, 423, 355, 297, 305, 407, 484, 343, 345, 315, 301, 465, 357, 327, 356, 309, 379, 458, 289, 426, 512, 337, 360, 368, 335, 401, 420, 372, 349, 443, 422, 257, 329, 396, 483, 384, 285, 411, 358, 396, 367, 345, 355, 441, 570, 271, 388, 345, 407, 418, 391, 368, 355, 434, 368, 306, 404, 386, 357, 326, 423, 423, 404, 402, 328, 448, 459, 436, 324, 343, 375, 320, 599, 423, 326, 340, 399, 496, 391, 447, 505, 449, 402, 476, 323, 539, 444, 482, 488, 359, 457, 414, 295, 556, 310, 396, 378, 364, 391, 459, 313, 371, 375, 436, 405, 410, 425, 364, 307, 399, 393, 406, 371, 448, 288, 360, 405, 341, 406, 387, 373, 449, 487, 380, 297, 354, 428, 331, 431, 351, 364, 372, 349, 457, 347, 532, 365, 337, 334, 421, 363, 346, 352, 465, 483, 408, 376, 355, 453, 311, 367, 356, 429, 361, 466, 389, 395, 329, 454, 325, 321, 232, 447, 352, 408, 440, 338, 434, 428, 364, 367, 397, 309, 282, 347, 461, 497, 463, 435, 478, 315, 320, 566, 342, 396, 418, 343, 413, 445, 389, 411, 469, 337, 479, 358, 394, 473, 357, 445, 425, 429, 306, 368, 322, 472, 446, 288, 364, 352, 355, 344, 388, 392, 373, 331, 352, 403, 506, 372, 352, 435, 378, 413, 476, 534, 353, 380, 353, 407, 362, 375, 324, 337, 305, 532, 449, 409, 378, 426, 318, 501, 422, 330, 478, 398, 424, 330, 376, 459, 306, 349, 337, 396, 364, 317, 396, 381, 395, 442, 404, 387, 340, 247, 389, 414, 393, 439, 272, 391, 344, 420, 373, 341, 308, 406, 395, 314, 326, 471, 434, 465, 362, 481, 352, 333, 428, 338, 291, 341, 253, 412, 468, 284, 305, 361, 311, 287, 415, 363, 472, 345, 353, 398, 344, 365, 338, 548, 321, 351, 293, 267, 336, 525, 253, 333, 313, 410, 420, 311, 389, 344, 378, 390, 405, 421, 307, 358, 338, 418, 435, 330, 462, 408, 450, 379, 352, 368, 288, 362, 440, 292, 366, 335, 522, 456, 409, 422, 390, 393, 356, 346, 313, 435, 312, 311, 488, 488, 309, 349, 387, 268, 460, 409, 284, 332, 469, 360, 357, 272, 376, 385, 334, 327, 312, 446, 326, 380, 460, 346, 296, 399, 282, 369, 317, 346, 400, 403, 448, 355, 413, 332, 396, 493, 465, 392, 489, 337, 447, 387, 419, 445, 605, 282, 308, 291, 436, 268, 346, 352, 474, 286, 395, 391, 313, 323, 325, 326, 416, 384, 430, 490, 275, 347, 341, 407, 369, 400, 389, 484, 334, 421, 417, 331, 527, 386, 439, 369, 299, 356, 484, 382, 355, 393, 413, 315, 397, 415, 313, 462, 448, 264, 407, 403, 520, 389, 391, 377, 391, 233, 354, 307, 364, 444, 390, 379, 372, 341, 477, 500, 396, 444, 460, 556, 415, 276, 347, 397, 274, 417, 446, 380, 244, 284, 422, 332, 380, 483, 488, 363, 346, 423, 311, 376, 381, 305, 375, 368, 428, 422, 410, 484, 424, 416, 362, 374, 346, 479, 470, 397, 396, 429, 356, 421, 270, 562, 463, 283, 447, 311, 393, 561, 568, 487, 430, 412, 416, 474, 477, 356, 342, 322, 439, 392, 304, 283, 349, 397, 372, 374, 431, 367, 429, 518, 367, 455, 495, 341, 349, 476, 382, 289, 440, 333, 445, 345, 293, 276, 279, 332, 478, 426, 363, 356, 325, 460, 427, 432, 460, 377, 373, 353, 443, 377, 428, 445, 339, 425, 371, 331, 303, 295, 399, 427, 409, 373, 337, 319, 291, 435, 445, 498, 409, 376, 341, 538, 338, 379, 478, 342, 354, 309, 579, 342, 341, 416, 374, 430, 408, 505, 373, 379, 412, 356, 401, 405, 314, 445, 296, 351, 567, 416, 393, 369, 372, 351, 415, 365, 370, 285, 443, 470, 492, 332, 390, 378, 363, 332, 404, 390, 307, 385, 362, 425, 481, 435, 275, 487, 344, 401, 450, 456, 423, 343, 359, 329, 352, 446, 305, 309, 330, 419, 450, 330, 403, 527, 385, 425, 335, 375, 422, 353, 279, 487, 550, 283, 350, 388, 430, 297, 541, 393, 350, 362, 420, 291, 322, 359, 344, 390, 367, 367, 381, 415, 506, 371, 385, 618, 475, 345, 410, 292, 467, 404, 440, 325, 439, 422, 424, 476, 576, 397, 446, 373, 337, 404, 305, 497, 525, 335, 373, 380, 437, 463, 346, 392, 278, 309, 361, 424, 397, 417, 332, 474, 292, 369, 457, 389, 342, 461, 381, 325, 306, 285, 437, 460, 313, 412, 341, 319, 396, 347, 549, 373, 449, 329, 370, 471, 423, 440, 362, 362, 378, 353, 420, 464, 273, 334, 374, 416, 387, 385, 457, 383, 445, 310, 310, 392, 369, 351, 458, 373, 362, 419, 428, 404, 350, 543, 474, 373, 424, 572, 458, 331, 460, 441, 467, 489, 430, 356, 389, 437, 413, 428, 416, 472, 416, 429, 416, 275, 413, 295, 434, 395, 325, 360, 392, 462, 447, 401, 392, 408, 486, 394, 468, 535, 318, 381, 349, 393, 394, 666, 345, 405, 443, 353, 445, 526, 467, 405, 307, 426, 366, 379, 359, 387, 377, 446, 281, 361, 405, 350, 408, 537, 423, 336, 381, 261, 418, 306, 352, 461, 390, 380, 446, 381, 304, 295, 306, 338, 324, 329, 365, 592, 436, 397, 515, 359, 437, 445, 292, 463, 395, 400, 473, 354, 424, 318, 490, 413, 336, 370, 402, 453, 399, 308, 286, 532, 427, 340, 401, 432, 299, 540, 622, 385, 305, 414, 413, 554, 344, 373, 376, 320, 483, 345, 301, 392, 349, 476, 319, 384, 363, 447, 467, 291, 289, 312, 463, 462, 342, 386, 528, 345, 379, 586, 484, 310, 338, 370, 243, 426, 554, 411, 349, 509, 303, 369, 395, 460, 452, 436, 409, 318, 399, 388, 427, 371, 334, 256, 407, 366, 319, 452, 344, 510, 449, 319, 372, 403, 420, 397, 379, 368, 406, 378, 379, 372, 240, 361, 466, 273, 458, 460, 309, 390, 361, 373, 335, 464, 441, 376, 412, 295, 351, 331, 399, 350, 539, 320, 360, 322, 446, 455, 408, 361, 384, 394, 217, 402, 431, 446, 285, 426, 418, 281, 392, 495, 376, 367, 392, 382, 449, 366, 368, 281, 346, 441, 348, 360, 335, 395, 397, 367, 318, 335, 342, 332, 350, 432, 387, 376, 420, 363, 358, 394, 340, 307, 491, 463, 543, 408, 337, 360, 297, 413, 395, 531, 442, 337, 401, 406, 337, 548, 365, 426, 439, 396, 323, 332, 340, 532, 356, 357, 342, 335, 522, 410, 457, 330, 475, 343, 421, 397, 315, 414, 365, 345, 373, 343, 311, 447, 403, 469, 367, 466, 319, 360, 434, 481, 425, 305, 391, 441, 388, 425, 406, 309, 312, 341, 333, 403, 457, 432, 347, 321, 362, 392, 412, 465, 324, 474, 279, 275, 387, 341, 406, 393, 404, 368, 348, 375, 378, 415, 351, 346, 414, 324, 320, 332, 426, 407, 425, 340, 339, 470, 408, 352, 305, 391, 305, 357, 423, 397, 350, 428, 530, 379, 457, 351, 432, 268, 421, 357, 355, 417, 369, 609, 394, 336, 230, 324, 376, 343, 398, 334, 382, 454, 376, 323, 484, 307, 426, 350, 383, 511, 336, 356, 317, 422, 337, 412, 335, 333, 485, 342, 408, 505, 355, 325, 243, 352, 299, 463, 369, 450, 369, 475, 410, 356, 512, 470, 373, 408, 388, 322, 432, 388, 403, 385, 359, 430, 444, 420, 410, 365, 362, 492, 422, 454, 449, 303, 355, 339, 469, 456, 396, 428, 379, 439, 332, 433, 410, 406, 365, 324, 426, 440, 301, 365, 375, 357, 262, 385, 314, 496, 438, 372, 355, 401, 316, 397, 313, 431, 354, 296, 411, 306, 306, 400, 557, 472, 309, 434, 446, 353, 311, 449, 313, 361, 296, 368, 384, 363, 462, 324, 416, 519, 429, 294, 394, 365, 396, 385, 379, 328, 383, 302, 341, 339, 294, 364, 289, 454, 265, 412, 324, 434, 400, 388, 506, 407, 392, 362, 435, 478, 294, 372, 331, 529, 367, 314, 391, 368, 411, 380, 459, 344, 433, 464, 465, 354, 416, 340, 357, 344, 330, 416, 253, 355, 436, 547, 400, 328, 356, 362, 319, 348, 522, 468, 298, 365, 391, 418, 537, 370, 361, 396, 402, 344, 328, 294, 371, 420, 413, 368, 430, 295, 537, 416, 390, 376, 402, 441, 372, 399, 349, 431, 329, 424, 368, 419, 486, 296, 333, 388, 360, 376, 322, 357, 461, 383, 386, 402, 480, 535, 377, 366, 380, 371, 307, 297, 413, 396, 306, 412, 513, 443, 311, 309, 533, 305, 436, 348, 445, 426, 419, 451, 466, 413, 432, 349, 381, 389, 374, 405, 402, 417, 402, 427, 441, 320, 364, 347, 364, 425, 357, 334, 331, 488, 497, 379, 414, 392, 487, 356, 324, 382, 429, 449, 318, 516, 273, 349, 546, 343, 347, 346, 424, 416, 535, 369, 391, 310, 252, 356, 283, 284, 369, 342, 473, 427, 403, 427, 357, 570, 396, 305, 343, 298, 368, 422, 398, 382, 408, 352, 449, 424, 389, 377, 317, 331, 455, 538, 379, 378, 403, 337, 376, 361, 350, 423, 388, 293, 361, 380, 443, 400, 293, 295, 387, 332, 436, 416, 452, 381, 391, 325, 372, 460, 348, 315, 316, 348, 365, 332, 402, 289, 360, 299, 412, 385, 317, 389, 469, 422, 326, 337, 426, 425, 482, 245, 395, 390, 332, 394, 315, 404, 321, 374, 436, 465, 381, 438, 388, 556, 345, 354, 427, 314, 274, 388, 436, 406, 318, 359, 449, 403, 356, 301, 430, 258, 357, 376, 416, 432, 439, 443, 338, 464, 367, 355, 330, 454, 402, 452, 278, 314, 486, 356, 374, 334, 290, 359, 425, 338, 443, 435, 448, 436, 419, 399, 357, 252, 261, 479, 389, 337, 489, 358, 281, 346, 361, 315, 482, 434, 393, 413, 309, 366, 536, 389, 428, 335, 324, 360, 378, 443, 334, 313, 379, 355, 325, 403, 376, 291, 499, 360, 370, 329, 356, 401, 314, 329, 323, 375, 307, 481, 428, 397, 568, 354, 322, 393, 342, 354, 254, 399, 384, 348, 429, 339, 426, 305, 420, 364, 395, 481, 334, 399, 391, 367, 410, 385, 286, 374, 337, 355, 319, 421, 401, 340, 447]

# print(statistics.mean(ppi_total))
# print(statistics.pstdev(ppi_total))



with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\COVID19_GDDS_nodes.csv', 'r') as nodecsv: # Open the file                       
	nodereader = csv.reader(nodecsv) # Read the csv  
	nodes = [n for n in nodereader][1:]                     
	node_names = [n[0] for n in nodes] # Get a list of only the node names 

with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\COVID19_GDDS_edges.csv', 'r') as edgecsv: # Open the file
	edgereader = csv.reader(edgecsv) # Read the csv     
	edges = [tuple(e) for e in edgereader][1:] # Retrieve the data

G_covid = nx.Graph()
G_covid.add_nodes_from(node_names)
G_covid.add_edges_from(edges)

descr_dict_covid = {}
for node in nodes: 
	descr_dict_covid[node[0]] = node[2]

ppi = 0
for u,v,c in G_covid.edges(data=True) :	
	if (descr_dict_covid[u] == 'Human PPI (target)' and descr_dict_covid[v] == 'Human PPI (target)') :
		ppi = ppi + 1


print("ppi", ppi)	

stop = timeit.default_timer()
print(stop - start)


convert_pval_manlio.py,
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

convert_pval_manlio_adjust.py,

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

get_structural_degrees_distributions_across_mock_networks.py,

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
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
from collections import Counter
from collections import defaultdict
from collections import OrderedDict
import statistics 
import re
import gzip



with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\COVID19_GDDS_nodes.csv', 'r') as nodecsv: # Open the file                       
		nodereader = csv.reader(nodecsv) # Read the csv  
		nodes = [n for n in nodereader][1:]                     
		node_names = [n[0] for n in nodes] # Get a list of only the node names 

descr_dict_covid = {}
for node in nodes: 
	descr_dict_covid[node[0]] = node[2]

# covid_proteins =  [node for node in descr_dict_covid if descr_dict_covid[node] == 'Human PPI (target)']
# covid_GO = [node for node in descr_dict_covid if descr_dict_covid[node] == 'GO']
# covid_symptom = [node for node in descr_dict_covid if descr_dict_covid[node] == 'Symptom']
# covid_disease = [node for node in descr_dict_covid if descr_dict_covid[node] == 'Disease']
covid_drug = [node for node in descr_dict_covid if descr_dict_covid[node] == 'Drug']


nodefiles_numbers = []
edgefiles_numbers = []
for filename in os.listdir("A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\random_networks\\Mock_networks_21Apr"):
	if 'nodes_' in filename :
		filename = filename.split(".csv")[0]
		nodefiles_numbers.append(int(filename[19:]))
	if 'edges_' in filename :
		filename = filename.split(".csv")[0]
		edgefiles_numbers.append(int(filename[19:]))

nodefiles_numbers = sorted(nodefiles_numbers)
edgefiles_numbers = sorted(edgefiles_numbers)

networks_to_avoid = [3, 11, 22, 23, 33, 34, 38, 57, 72, 77, 123, 127, 148, 156, 157, 183, 205, 209, 210, 217, 232, 252, 268, 281, 291, 312, 341, 353, 362, 364, 413, 421, 433, 443, 451, 494, 525, 532, 569, 592, 629, 658, 662, 674, 676, 701, 705, 712, 736, 757, 788, 796, 825, 827, 838, 873, 876, 908, 910, 916, 917, 918, 938, 958, 965, 1088, 1125, 1128, 1145, 1146, 1165, 1171, 1183, 1185, 1192, 1195, 1211, 1244, 1254, 1257, 1300, 1353, 1354, 1364, 1372, 1394, 1395, 1401, 1410, 1427, 1454, 1477, 1497, 1517, 1532, 1549, 1551, 1575, 1578, 1611, 1624, 1629, 1637, 1650, 1655, 1697, 1700, 1713, 1736, 1745, 1762, 1790, 1800, 1815, 1824, 1827, 1866, 1872, 1886, 1887, 2053, 2072, 2106, 2115, 2128, 2160, 2217, 2343, 2383, 2388, 2413, 2458, 2465]
nodefiles_numbers = [item for item in nodefiles_numbers if item not in networks_to_avoid]
print(nodefiles_numbers)


### for covid proteins
# covid_proteins_structural_degree_distribution,covid_proteins_structural_strength_distribution = {},{} 
# for iteration in nodefiles_numbers[1:] :
# 	print("\n########################## ITERATION %s ###########################\n" % iteration)
# 	with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\Further_Analysis\\full_analysis_for_source_to_proteins_undirected\\Human PPI (target)-prot\\COVID19_GDDS_Human PPI (target)_degrees_to_proteins_%s.tsv' % str(iteration), 'r') as file_read :                 
# 		data = file_read.readlines()
# 		for line in data :
# 			line = line.split('\t')
# 			protein = line[0]
# 			if protein in covid_proteins :
# 				structural_degree = int(line[1])
# 				structural_strength = float(line[2].replace('\n',''))
# 				if protein not in covid_proteins_structural_degree_distribution :
# 					covid_proteins_structural_degree_distribution[protein] = []
# 				if protein not in covid_proteins_structural_strength_distribution :
# 					covid_proteins_structural_strength_distribution[protein] = []
# 				covid_proteins_structural_degree_distribution[protein].append(structural_degree)
# 				covid_proteins_structural_strength_distribution[protein].append(structural_strength)
# with open('pvalues_from_zscores\\covid_proteins_structural_degree_distribution_in_mock_networks.tsv', 'a+') as file_write :
# 	file_write.write('covid_protein\tstructural_degree_in_mock_networks\n')
# 	for protein in covid_proteins_structural_degree_distribution :
# 		file_write.write('%s\t' % protein)
# 		for degree in covid_proteins_structural_degree_distribution[protein][:-1] :
# 			file_write.write('%s\t' % degree)
# 		file_write.write('%s\n' % covid_proteins_structural_degree_distribution[protein][-1])
# with open('pvalues_from_zscores\\covid_proteins_structural_strength_distribution_in_mock_networks.tsv', 'a+') as file_write :
# 	file_write.write('covid_protein\tstructural_strength_in_mock_networks\n')
# 	for protein in covid_proteins_structural_strength_distribution :
# 		file_write.write('%s\t' % protein)
# 		for degree in covid_proteins_structural_strength_distribution[protein][:-1] :
# 			file_write.write('%s\t' % degree)
# 		file_write.write('%s\n' % covid_proteins_structural_strength_distribution[protein][-1])


### for covid GO
# covid_GO_degree_to_proteins_distribution,covid_GO_structural_strength_to_proteins_distribution = {},{} 
# for iteration in nodefiles_numbers[1:] :
# 	print("\n########################## ITERATION %s ###########################\n" % iteration)
# 	with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\Further_Analysis\\full_analysis_for_source_to_proteins_undirected\\GO-prot\\COVID19_GDDS_GO_degrees_to_proteins_%s.tsv' % str(iteration), 'r') as file_read :                 
# 		data = file_read.readlines()
# 		for line in data :
# 			line = line.split('\t')
# 			GO = line[0]
# 			if GO in covid_GO :
# 				degree_to_proteins = int(line[1])
# 				structural_strength_to_proteins = float(line[2].replace('\n',''))
# 				if GO  not in covid_GO_degree_to_proteins_distribution :
# 					covid_GO_degree_to_proteins_distribution[GO] = []
# 				if GO not in covid_GO_structural_strength_to_proteins_distribution :
# 					covid_GO_structural_strength_to_proteins_distribution[GO] = []
# 				covid_GO_degree_to_proteins_distribution[GO].append(degree_to_proteins)
# 				covid_GO_structural_strength_to_proteins_distribution[GO].append(structural_strength_to_proteins)
# with open('pvalues_from_zscores\\covid_GO_degree_to_proteins_distribution_in_mock_networks.tsv', 'a+') as file_write :
# 	file_write.write('covid_GO\tdegree_to_proteins_in_mock_networks\n')
# 	for GO in covid_GO_degree_to_proteins_distribution :
# 		file_write.write('%s\t' % GO)
# 		for degree in covid_GO_degree_to_proteins_distribution[GO][:-1] :
# 			file_write.write('%s\t' % degree)
# 		file_write.write('%s\n' % covid_GO_degree_to_proteins_distribution[GO][-1])
# with open('pvalues_from_zscores\\covid_GO_structural_strength_to_proteins_distribution_in_mock_networks.tsv', 'a+') as file_write :
# 	file_write.write('covid_GO\tstructural_strength_to_proteins_in_mock_networks\n')
# 	for GO in covid_GO_structural_strength_to_proteins_distribution :
# 		file_write.write('%s\t' % GO)
# 		for degree in covid_GO_structural_strength_to_proteins_distribution[GO][:-1] :
# 			file_write.write('%s\t' % degree)
# 		file_write.write('%s\n' % covid_GO_structural_strength_to_proteins_distribution[GO][-1])


# ####for covid symptom
# covid_symptom_degree_to_proteins_distribution,covid_symptom_structural_strength_to_proteins_distribution = {},{} 
# for iteration in nodefiles_numbers[1:] :
# 	print("\n########################## ITERATION %s ###########################\n" % iteration)
# 	with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\Further_Analysis\\full_analysis_for_source_to_proteins_undirected\\Symptom-prot\\COVID19_GDDS_symptom_degrees_to_proteins_%s.tsv' % str(iteration), 'r') as file_read :                 
# 		data = file_read.readlines()
# 		for line in data :
# 			line = line.split('\t')
# 			symptom = line[0]
# 			if symptom in covid_symptom :
# 				degree_to_proteins = int(line[1])
# 				structural_strength_to_proteins = float(line[2].replace('\n',''))
# 				if symptom  not in covid_symptom_degree_to_proteins_distribution :
# 					covid_symptom_degree_to_proteins_distribution[symptom] = []
# 				if symptom not in covid_symptom_structural_strength_to_proteins_distribution :
# 					covid_symptom_structural_strength_to_proteins_distribution[symptom] = []
# 				covid_symptom_degree_to_proteins_distribution[symptom].append(degree_to_proteins)
# 				covid_symptom_structural_strength_to_proteins_distribution[symptom].append(structural_strength_to_proteins)
# with open('pvalues_from_zscores\\covid_symptom_degree_to_proteins_distribution_in_mock_networks.tsv', 'a+') as file_write :
# 	file_write.write('covid_symptom\tdegree_to_proteins_in_mock_networks\n')
# 	for symptom in covid_symptom_degree_to_proteins_distribution :
# 		file_write.write('%s\t' % symptom)
# 		for degree in covid_symptom_degree_to_proteins_distribution[symptom][:-1] :
# 			file_write.write('%s\t' % degree)
# 		file_write.write('%s\n' % covid_symptom_degree_to_proteins_distribution[symptom][-1])
# with open('pvalues_from_zscores\\covid_symptom_structural_strength_to_proteins_distribution_in_mock_networks.tsv', 'a+') as file_write :
# 	file_write.write('covid_symptom\tstructural_strength_to_proteins_in_mock_networks\n')
# 	for symptom in covid_symptom_structural_strength_to_proteins_distribution :
# 		file_write.write('%s\t' % symptom)
# 		for degree in covid_symptom_structural_strength_to_proteins_distribution[symptom][:-1] :
# 			file_write.write('%s\t' % degree)
# 		file_write.write('%s\n' % covid_symptom_structural_strength_to_proteins_distribution[symptom][-1])


# ## for covid disease		
# covid_disease_degree_to_proteins_distribution,covid_disease_structural_strength_to_proteins_distribution = {},{} 
# for iteration in nodefiles_numbers[1:] :
# 	print("\n########################## ITERATION %s ###########################\n" % iteration)
# 	with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\Further_Analysis\\full_analysis_for_source_to_proteins_undirected\\Disease-prot\\COVID19_GDDS_Disease_degrees_to_proteins_%s.tsv' % str(iteration), 'r') as file_read :                 
# 		data = file_read.readlines()
# 		for line in data :
# 			line = line.split('\t')
# 			disease = line[0]
# 			if disease in covid_disease :
# 				degree_to_proteins = int(line[1])
# 				structural_strength_to_proteins = float(line[2].replace('\n',''))
# 				if disease  not in covid_disease_degree_to_proteins_distribution :
# 					covid_disease_degree_to_proteins_distribution[disease] = []
# 				if disease not in covid_disease_structural_strength_to_proteins_distribution :
# 					covid_disease_structural_strength_to_proteins_distribution[disease] = []
# 				covid_disease_degree_to_proteins_distribution[disease].append(degree_to_proteins)
# 				covid_disease_structural_strength_to_proteins_distribution[disease].append(structural_strength_to_proteins)
# with open('pvalues_from_zscores\\covid_disease_degree_to_proteins_distribution_in_mock_networks.tsv', 'a+') as file_write :
# 	file_write.write('covid_disease\tdegree_to_proteins_in_mock_networks\n')
# 	for disease in covid_disease_degree_to_proteins_distribution :
# 		file_write.write('%s\t' % disease)
# 		for degree in covid_disease_degree_to_proteins_distribution[disease][:-1] :
# 			file_write.write('%s\t' % degree)
# 		file_write.write('%s\n' % covid_disease_degree_to_proteins_distribution[disease][-1])
# with open('pvalues_from_zscores\\covid_disease_structural_strength_to_proteins_distribution_in_mock_networks.tsv', 'a+') as file_write :
# 	file_write.write('covid_disease\tstructural_strength_to_proteins_in_mock_networks\n')
# 	for disease in covid_disease_structural_strength_to_proteins_distribution :
# 		file_write.write('%s\t' % disease)
# 		for degree in covid_disease_structural_strength_to_proteins_distribution[disease][:-1] :
# 			file_write.write('%s\t' % degree)
# 		file_write.write('%s\n' % covid_disease_structural_strength_to_proteins_distribution[disease][-1])




## for covid drug		
covid_drug_degree_to_proteins_distribution,covid_drug_structural_strength_to_proteins_distribution = {},{} 
for iteration in nodefiles_numbers[1:] :
	print("\n########################## ITERATION %s ###########################\n" % iteration)
	with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\Further_Analysis\\full_analysis_for_source_to_proteins_undirected\\Drug-prot\\COVID19_GDDS_Drug_degrees_to_proteins_%s.tsv' % str(iteration), 'r') as file_read :                 
		data = file_read.readlines()
		for line in data :
			line = line.split('\t')
			drug = line[0]
			if drug in covid_drug :
				degree_to_proteins = int(line[1])
				structural_strength_to_proteins = float(line[2].replace('\n',''))
				if drug  not in covid_drug_degree_to_proteins_distribution :
					covid_drug_degree_to_proteins_distribution[drug] = []
				if drug not in covid_drug_structural_strength_to_proteins_distribution :
					covid_drug_structural_strength_to_proteins_distribution[drug] = []
				covid_drug_degree_to_proteins_distribution[drug].append(degree_to_proteins)
				covid_drug_structural_strength_to_proteins_distribution[drug].append(structural_strength_to_proteins)
with open('covid_drug_degree_to_proteins_distribution_in_mock_networks.tsv', 'a+') as file_write :
	file_write.write('covid_drug\tdegree_to_proteins_in_mock_networks\n')
	for drug in covid_drug_degree_to_proteins_distribution :
		file_write.write('%s\t' % drug)
		for degree in covid_drug_degree_to_proteins_distribution[drug][:-1] :
			file_write.write('%s\t' % degree)
		file_write.write('%s\n' % covid_drug_degree_to_proteins_distribution[drug][-1])
with open('covid_drug_structural_strength_to_proteins_distribution_in_mock_networks.tsv', 'a+') as file_write :
	file_write.write('covid_drug\tstructural_strength_to_proteins_in_mock_networks\n')
	for drug in covid_drug_structural_strength_to_proteins_distribution :
		file_write.write('%s\t' % drug)
		for degree in covid_drug_structural_strength_to_proteins_distribution[drug][:-1] :
			file_write.write('%s\t' % degree)
		file_write.write('%s\n' % covid_drug_structural_strength_to_proteins_distribution[drug][-1])







# # dirname = "full_analysis_for_source_to_proteins_undirected"
# # os.mkdir(dirname)
# source_list = ['Symptom', 'Disease','Drug', 'Human PPI (target)']
# # source_list = ['GO']

# with open('Log_full_analysis_for_source_to_proteins_undirected.txt', 'a+') as full_analysis_log_write :

# 	#####################################################################################################
# 	#########################  Compute network source to protein degree distribution ####################
# 	##########################################  based on edges to protein ###############################
# 	#####################################################################################################

# 	with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\COVID19_GDDS_nodes.csv', 'r') as nodecsv: # Open the file                       
# 		nodereader = csv.reader(nodecsv) # Read the csv  
# 		nodes = [n for n in nodereader][1:]                     
# 		node_names = [n[0] for n in nodes] # Get a list of only the node names 

# 	with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\COVID19_GDDS_edges.csv', 'r') as edgecsv: # Open the file
# 		edgereader = csv.reader(edgecsv) # Read the csv     
# 		edges = [tuple(e) for e in edgereader][1:] # Retrieve the data

# 	G_covid = nx.Graph()
# 	G_covid.add_nodes_from(node_names)
# 	G_covid.add_edges_from(edges)

# 	descr_dict_covid = {}
# 	for node in nodes: 
# 		descr_dict_covid[node[0]] = node[2]

# 	print(nx.info(G_covid))
# 	# full_analysis_log_write.write(nx.info(G_covid))
# 	# Number of nodes: 16007
# 	# Number of edges: 53917
# 	# Average degree:   6.7367

# 	# print(nx.degree(G_covid)) # this is the list containing as many tuples as nodes, indicating each node's degree.

# 	average_degree = statistics.mean([ tpl[1] for tpl in nx.degree(G_covid) ])
# 	# print(average_degree) # 6.736677703504717

# 	for source_node in source_list :

# 		source_degrees = nx.degree(G_covid, nbunch=[n for n in G_covid.nodes if descr_dict_covid[n] == source_node])
# 		# print("degree of source nodes only", source_degrees)

# 		average_degree_of_sources = statistics.mean([ tpl[1] for tpl in source_degrees ])
# 		print("average degree of %s nodes" % source_node, average_degree_of_sources) # 
# 		# full_analysis_log_write.write("average degree of %s nodes : %s\n" % (source_node, average_degree_of_sources))
# 	# average degree of GO nodes 2.64123888729567
# 	# average degree of Symptom nodes 6.644413537320353
# 	# average degree of Disease nodes 4.329980842911877
# 	# average degree of Drug nodes 2.2826582500438364
# 	# average degree of Human PPI (target) nodes 115.66739606126914



# 	# # # # # #####################################################################################################
# 	# # # # # ################################# CALCULATION ON EDGES WITH FILTERING source - PROT #####################
# 	# # # # # #####################################################################################################

# 	for source_node in source_list :
# 		subfoldername = 'full_analysis_for_source_to_proteins_undirected/%s-prot' % (source_node)
# 		os.mkdir(subfoldername)
# 		print("Searching %s - Human PPI (target) links\n" % source_node)
# 		full_analysis_log_write.write("Searching %s - Human PPI (target) links\n" % source_node)
		
# 		# now we keep only the edges of source_node that interact with proteins
# 		sources_linked_to_proteins = []
# 		for u,v,c in G_covid.edges(data=True) :	
# 			if (descr_dict_covid[u] == source_node and descr_dict_covid[v] == 'Human PPI (target)') :
# 				sources_linked_to_proteins.append(u)
# 			if (descr_dict_covid[v] == source_node and descr_dict_covid[u] == 'Human PPI (target)') :
# 				sources_linked_to_proteins.append(v)	
# 		print("len(%s_linked_to_proteins), len(set(%s_linked_to_proteins))" % (source_node, source_node), len(sources_linked_to_proteins), len(set(sources_linked_to_proteins))) #
# 		full_analysis_log_write.write("len(%s_linked_to_proteins), len(set(%s_linked_to_proteins)) : %s, %s\n" % (source_node, source_node, len(sources_linked_to_proteins), len(set(sources_linked_to_proteins)))) #


# 		source_counts = sorted(sources_linked_to_proteins)

# 		source_counts_dict = defaultdict( int )
# 		for source in  source_counts:
# 		    source_counts_dict[source] += 1

# 		my_sources_nodes = [n for n in G_covid.nodes if descr_dict_covid[n] == source_node]
# 		for source in  my_sources_nodes :
# 			if source_counts_dict.get(source) == None :
# 				source_counts_dict[source]= 0

# 		average_degree_of_my_sources_nodes = statistics.mean([ source_counts_dict[source] for source in  source_counts_dict ])
# 		print("average_degree_of_my_sources_nodes", average_degree_of_my_sources_nodes)
# 		full_analysis_log_write.write("average_degree_of_my_%s_nodes : %s\n" % (source_node, average_degree_of_my_sources_nodes))

# 		max_degree_of_my_sources_nodes = max([ source_counts_dict[source] for source in  source_counts_dict ])
# 		print("max_degree_of_my_sources_nodes", max_degree_of_my_sources_nodes)

# 		ordered_source_counts_dict = OrderedDict(sorted(source_counts_dict.items()))


# 		with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\Further_Analysis\\full_analysis_for_source_to_proteins_undirected\\COVID19_GDDS_%s_degrees_to_proteins.tsv' % source_node, 'a+') as file_write:
# 			for source in ordered_source_counts_dict :
# 				file_write.write("%s\t%s\t%s\n" % (source, ordered_source_counts_dict[source], ordered_source_counts_dict[source] / sum(ordered_source_counts_dict.values())))


# 		# # # # # ##################################################################
# 		# # # # # ############ SAME ANALYSIS FOR RANDOM NETWORKS ###################
# 		# # # # # # ##################################################################
		
# 		with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\Further_Analysis\\full_analysis_for_source_to_proteins_undirected\\Log_degree_analysis_on_random_networks_for_%s_proteins.txt' % source_node, 'a+') as log_write:
# 			log_write.write("In covid net : len(sources_linked_to_proteins) and len(set(sources_linked_to_proteins)) are : %s and %s \n" % (len(sources_linked_to_proteins), len(set(sources_linked_to_proteins))) ) # # 9208 435
# 			log_write.write("In covid net : average_degree_of_my_sources_nodes = %s\n" % average_degree_of_my_sources_nodes)
# 			log_write.write("In covid net : max_degree_of_my_sources_nodes = %s\n" % max_degree_of_my_sources_nodes)

# 			nodefiles_numbers = []
# 			edgefiles_numbers = []
# 			for filename in os.listdir("A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\random_networks\\Mock_networks_21Apr"):
# 				if 'nodes_' in filename :
# 					filename = filename.split(".csv")[0]
# 					nodefiles_numbers.append(int(filename[19:]))
# 				if 'edges_' in filename :
# 					filename = filename.split(".csv")[0]
# 					edgefiles_numbers.append(int(filename[19:]))

# 			# log_write.write("len(nodefiles_numbers) = %s\n" % len(nodefiles_numbers)) #2051
# 			# log_write.write("len(edgefiles_numbers) = %s\n" % len(edgefiles_numbers)) #2051
# 			# log_write.write("difference between edgefiles_numbers and nodefiles_numbers = %s\n" % len(set(nodefiles_numbers).difference(set(edgefiles_numbers)))) #0

# 			nodefiles_numbers = sorted(nodefiles_numbers)
# 			edgefiles_numbers = sorted(edgefiles_numbers)

# 			networks_to_avoid = [3, 11, 22, 23, 33, 34, 38, 57, 72, 77, 123, 127, 148, 156, 157, 183, 205, 209, 210, 217, 232, 252, 268, 281, 291, 312, 341, 353, 362, 364, 413, 421, 433, 443, 451, 494, 525, 532, 569, 592, 629, 658, 662, 674, 676, 701, 705, 712, 736, 757, 788, 796, 825, 827, 838, 873, 876, 908, 910, 916, 917, 918, 938, 958, 965, 1088, 1125, 1128, 1145, 1146, 1165, 1171, 1183, 1185, 1192, 1195, 1211, 1244, 1254, 1257, 1300, 1353, 1354, 1364, 1372, 1394, 1395, 1401, 1410, 1427, 1454, 1477, 1497, 1517, 1532, 1549, 1551, 1575, 1578, 1611, 1624, 1629, 1637, 1650, 1655, 1697, 1700, 1713, 1736, 1745, 1762, 1790, 1800, 1815, 1824, 1827, 1866, 1872, 1886, 1887, 2053, 2072, 2106, 2115, 2128, 2160, 2217, 2343, 2383, 2388, 2413, 2458, 2465]
# 			nodefiles_numbers = [item for item in nodefiles_numbers if item not in networks_to_avoid]

# 			for iteration in nodefiles_numbers[1:] :
# 				print("\n########################## ITERATION %s -- %s ###########################\n" % (iteration, source_node))
# 				log_write.write("\n########################## ITERATION %s ###########################\n" % iteration)
			
# 				with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\random_networks\\Mock_networks_21Apr\\COVID19_GDDS_nodes_%s.csv' % str(iteration), 'r') as nodecsv: # Open the file                   
# 					nodereader = csv.reader(nodecsv) # Read the csv  
# 					nodes = [n for n in nodereader][1:]                     
# 					node_names = [n[0] for n in nodes] # Get a list of only the node names 

# 				with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\random_networks\\Mock_networks_21Apr\\COVID19_GDDS_edges_%s.csv' % str(iteration), 'r') as edgecsv: # Open the file
# 					edgereader = csv.reader(edgecsv) # Read the csv     
# 					edges = [tuple(e) for e in edgereader][1:] # Retrieve the data

# 					log_write.write("len(node_names) : %s\n" % len(node_names)) #16007 --> 13448
# 					log_write.write("len(edges) : %s\n" % len(edges)) #53917 --> 37835

# 				G_random = nx.Graph()
# 				G_random.add_nodes_from(node_names)
# 				G_random.add_edges_from(edges)

# 				descr_dict_random = {}
# 				description_random_set=set()
# 				for node in nodes: 
# 					descr_dict_random[node[0]] = node[2]
# 					description_random_set.add(node[2])

# 				# print(nx.info(G_random))
# 				log_write.write("\nRANDOM GRAPH %s\n" % iteration)
# 				log_write.write("%s\n" % nx.info(G_random))

# 				#####################################################################################################
# 				################################# GET EDGES CONNECTING sources TO proteins ################
# 				#####################################################################################################

# 				sources_linked_to_proteins = []
# 				for u,v,c in G_random.edges(data=True) :	
# 					if (descr_dict_random[u] == source_node and descr_dict_random[v] == 'Human PPI (target)') :
# 						sources_linked_to_proteins.append(u)
# 					if (descr_dict_random[v] == source_node and descr_dict_random[u] == 'Human PPI (target)') :
# 						sources_linked_to_proteins.append(v)		
				
# 				source_counts = sorted(sources_linked_to_proteins)

# 				source_counts_dict = defaultdict( int )
# 				for source in  source_counts:
# 				    source_counts_dict[source] += 1

# 				my_sources_nodes = [n for n in G_random.nodes if descr_dict_random[n] == source]
# 				for source in  my_sources_nodes :
# 					if source_counts_dict.get(source) == None :
# 						# print("%s is missing in the linked sources : Adding it with a frequency of 0.\n" % protein)
# 						source_counts_dict[source]= 0

# 				average_degree_of_my_sources_nodes = statistics.mean([ source_counts_dict[source] for source in  source_counts_dict ])
# 				log_write.write("Average degree of %s linked to Human PPI (target) : %s\n" % (source_node, average_degree_of_my_sources_nodes))

# 				ordered_source_counts_dict = OrderedDict(sorted(source_counts_dict.items()))

# 				with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\Further_Analysis\\full_analysis_for_source_to_proteins_undirected\\%s-prot\\COVID19_GDDS_%s_degrees_to_proteins_%s.tsv' % (source_node, source_node, str(iteration)), 'a+') as file_write:
# 					for source in ordered_source_counts_dict :
# 						if (sum(ordered_source_counts_dict.values()) != 0) :
# 							file_write.write("%s\t%s\t%s\n" % (source, ordered_source_counts_dict[source], ordered_source_counts_dict[source] / sum(ordered_source_counts_dict.values())))
# 						else :
# 							file_write.write("%s\t%s\t%s\n" % (source, ordered_source_counts_dict[source], ordered_source_counts_dict[source]))



# 		# # # # ##################################################################
# 		# # # # ############ COMPUTE DEGREES MEANS AND Z-SCORES ##################
# 		# # # # ######## INCLUDING ABSENCE OF source AS ZERO DEGREES ###########
# 		# # # # ##################################################################


# 		# # For covid network, compute the mean and the standard deviation of proteins degrees and normalized degrees
# 		covid_sources = []
# 		covid_source_degree_dict = {}
# 		covid_source_norm_degree_dict = {}

# 		total_source_degree_dict = defaultdict(list)
# 		total_source_normalized_degree_dict = defaultdict(list)

# 		with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\Further_Analysis\\full_analysis_for_source_to_proteins_undirected\\COVID19_GDDS_%s_degrees_to_proteins.tsv' % source_node, 'r') as data :
# 			for line in data :
# 				line = line.split('\t')
# 				source = line[0]
# 				covid_sources.append(source)
# 				source_degree = int(line[1])
# 				source_norm_degree = float(line[2].replace('\n', ''))
# 				if covid_source_degree_dict.get(source) == None :
# 					covid_source_degree_dict[source] = source_degree
# 				else :
# 					print("warning : %s duplicate ?!" % protein)
# 				covid_source_norm_degree_dict[source] = source_norm_degree

# 				if total_source_degree_dict.get(source) == None :
# 					total_source_degree_dict[source] = [source_degree]
# 				else :
# 					total_source_degree_dict[source].append(source_degree)

# 				if total_source_normalized_degree_dict.get(source) == None :
# 					total_source_normalized_degree_dict[source] = [source_norm_degree]
# 				else :
# 					total_source_normalized_degree_dict[source].append(source_norm_degree)

# 		mean_source_deg = statistics.mean(covid_source_degree_dict.values())
# 		mean_norm_source_deg = statistics.mean(covid_source_norm_degree_dict.values())
# 		sd_source_deg = statistics.pstdev(covid_source_degree_dict.values())
# 		sd_norm_source_deg = statistics.pstdev(covid_source_norm_degree_dict.values())
# 		print(mean_source_deg, sd_source_deg, mean_norm_source_deg, sd_norm_source_deg) #
# 		print(len(covid_sources),len(covid_source_degree_dict),len(covid_source_norm_degree_dict)) #
# 		# 20.14879649890591 18.9696222570052 0.002188183807439825 0.00206012405050013
# 		# 457 457 457






# 		nodefiles_numbers = []
# 		edgefiles_numbers = []
# 		for filename in os.listdir("A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\random_networks\\Mock_networks_21Apr"):
# 			if 'nodes_' in filename :
# 				filename = filename.split(".csv")[0]
# 				nodefiles_numbers.append(int(filename[19:]))
# 			if 'edges_' in filename :
# 				filename = filename.split(".csv")[0]
# 				edgefiles_numbers.append(int(filename[19:]))

# 		# log_write.write("len(nodefiles_numbers) = %s\n" % len(nodefiles_numbers)) #2051
# 		# log_write.write("len(edgefiles_numbers) = %s\n" % len(edgefiles_numbers)) #2051
# 		# log_write.write("difference between edgefiles_numbers and nodefiles_numbers = %s\n" % len(set(nodefiles_numbers).difference(set(edgefiles_numbers)))) #0
# 		print("len(nodefiles_numbers) = %s\n" % len(nodefiles_numbers)) #2051
# 		print("len(edgefiles_numbers) = %s\n" % len(edgefiles_numbers)) #2051
# 		print("difference between edgefiles_numbers and nodefiles_numbers = %s\n" % len(set(nodefiles_numbers).difference(set(edgefiles_numbers)))) #0

# 		# nodefiles_numbers = sorted(nodefiles_numbers)
# 		# edgefiles_numbers = sorted(edgefiles_numbers)

# 		networks_to_avoid = [3, 11, 22, 23, 33, 34, 38, 57, 72, 77, 123, 127, 148, 156, 157, 183, 205, 209, 210, 217, 232, 252, 268, 281, 291, 312, 341, 353, 362, 364, 413, 421, 433, 443, 451, 494, 525, 532, 569, 592, 629, 658, 662, 674, 676, 701, 705, 712, 736, 757, 788, 796, 825, 827, 838, 873, 876, 908, 910, 916, 917, 918, 938, 958, 965, 1088, 1125, 1128, 1145, 1146, 1165, 1171, 1183, 1185, 1192, 1195, 1211, 1244, 1254, 1257, 1300, 1353, 1354, 1364, 1372, 1394, 1395, 1401, 1410, 1427, 1454, 1477, 1497, 1517, 1532, 1549, 1551, 1575, 1578, 1611, 1624, 1629, 1637, 1650, 1655, 1697, 1700, 1713, 1736, 1745, 1762, 1790, 1800, 1815, 1824, 1827, 1866, 1872, 1886, 1887, 2053, 2072, 2106, 2115, 2128, 2160, 2217, 2343, 2383, 2388, 2413, 2458, 2465]
# 		nodefiles_numbers = [item for item in nodefiles_numbers if item not in networks_to_avoid]

# 		for iteration in nodefiles_numbers[1:] :
# 			# log_write.write("\n########################## ITERATION %s ###########################\n" % iteration)
# 			print("\n########################## ITERATION %s ###########################\n" % iteration)

# 			with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\Further_Analysis\\full_analysis_for_source_to_proteins_undirected\\%s-prot\\COVID19_GDDS_%s_degrees_to_proteins_%s.tsv' % (source_node, source_node, str(iteration)), 'r') as data: # Open the file  
# 				for line in data :
# 					line = line.split('\t')
# 					source = line[0]
# 					degree = float(line[1])
# 					normalized_degree = float(line[2].replace('\n', ''))

# 					if total_source_degree_dict.get(source) == None :
# 						total_source_degree_dict[source] = [degree]
# 					else :
# 						total_source_degree_dict[source].append(degree)

# 					if total_source_normalized_degree_dict.get(source) == None :
# 						total_source_normalized_degree_dict[source] = [normalized_degree]
# 					else :
# 						total_source_normalized_degree_dict[source].append(normalized_degree)

# 		print(len(total_source_degree_dict)) # 19928

# 		print(len(total_source_normalized_degree_dict)) # 19928

# 		### make the files with all calculations
# 		with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\Further_Analysis\\full_analysis_for_source_to_proteins_undirected\\COVID19_GDDS_Results_%s_to_Proteins-withZeros.tsv' % source_node, 'a+') as file_write:
# 			file_write.write("%s\t# of nets it is found in\tmin-degree\tmax-degree\tmean-degree\tmedian-degree\tsd-degree\n" % source_node)
# 			for source in  total_source_degree_dict :
# 				nets_number = len(total_source_degree_dict[source])
# 				min_deg = min(total_source_degree_dict[source])
# 				max_deg = max(total_source_degree_dict[source])
# 				mean_deg = statistics.mean(total_source_degree_dict[source])
# 				median_deg = statistics.median(total_source_degree_dict[source])
# 				sd_deg = statistics.pstdev(total_source_degree_dict[source])
# 				file_write.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (source, nets_number, min_deg, max_deg, mean_deg, median_deg, sd_deg))


# 		with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\Further_Analysis\\full_analysis_for_source_to_proteins_undirected\\COVID19_GDDS_Results_%s_to_Proteins_Relative-withZeros.tsv' % source_node, 'a+') as file_write:
# 			file_write.write("%s\t# of nets it is found in\tmin-degree-rel\tmax-degree-rel\tmean-degree-rel\tmedian-degree-rel\tsd-degree-rel\n" % source_node)
# 			for source in  total_source_normalized_degree_dict :
# 				nets_number = len(total_source_normalized_degree_dict[source])
# 				min_deg = min(total_source_normalized_degree_dict[source])
# 				max_deg = max(total_source_normalized_degree_dict[source])
# 				mean_deg = statistics.mean(total_source_normalized_degree_dict[source])
# 				median_deg = statistics.median(total_source_normalized_degree_dict[source])
# 				sd_deg = statistics.pstdev(total_source_normalized_degree_dict[source])
# 				file_write.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (source, nets_number, min_deg, max_deg, mean_deg, median_deg, sd_deg))

# 		### compute Z-score for the covid sources
# 		# For each source found in covid network, compute 2 Z-scores.

# 		with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\Further_Analysis\\full_analysis_for_source_to_proteins_undirected\\COVID19_GDDS_Results_%s_to_CovidProteins-withZeros.tsv' % source_node, 'a+') as file_write:
# 			file_write.write("covid %s\t# of nets it is found in\tmin-degree\tmax-degree\tmean-degree\tmedian-degree\tsd-degree\tcovid-degree\tZ-score\n" % source_node)
# 			for source in  covid_sources :
# 				nets_number = len(total_source_degree_dict[source])
# 				min_deg = min(total_source_degree_dict[source])
# 				max_deg = max(total_source_degree_dict[source])
# 				mean_deg = statistics.mean(total_source_degree_dict[source])
# 				median_deg = statistics.median(total_source_degree_dict[source])
# 				sd_deg = statistics.pstdev(total_source_degree_dict[source])
# 				covid_deg = covid_source_degree_dict[source]
# 				if sd_deg != 0 :
# 					z_score = ( covid_deg - mean_deg) / sd_deg
# 				else : 
# 					z_score = 'NaN'
# 				file_write.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (source, nets_number, min_deg, max_deg, mean_deg, median_deg, sd_deg, covid_deg, z_score))


# 		with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\Further_Analysis\\full_analysis_for_source_to_proteins_undirected\\COVID19_GDDS_Results_%s_to_CovidProteins_Relative-withZeros.tsv' % source_node, 'a+') as file_write:
# 			file_write.write("covid %s\t# of nets it is found in\tmin-degree-rel\tmax-degree-rel\tmean-degree-rel\tmedian-degree-rel\tsd-degree-rel\tcovid-degree-rel\tZ-score-rel\n" % source_node)
# 			for source in  covid_sources :
# 				nets_number = len(total_source_normalized_degree_dict[source])
# 				min_deg = min(total_source_normalized_degree_dict[source])
# 				max_deg = max(total_source_normalized_degree_dict[source])
# 				mean_deg = statistics.mean(total_source_normalized_degree_dict[source])
# 				median_deg = statistics.median(total_source_normalized_degree_dict[source])
# 				sd_deg = statistics.pstdev(total_source_normalized_degree_dict[source])
# 				covid_deg = covid_source_norm_degree_dict[source]
# 				if sd_deg != 0 :
# 					z_score = ( covid_deg - mean_deg) / sd_deg
# 				else : 
# 					z_score = 'NaN'
# 				file_write.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (source, nets_number, min_deg, max_deg, mean_deg, median_deg, sd_deg, covid_deg, z_score))










stop = timeit.default_timer()
print(stop - start)


plots.py,
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
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
from collections import Counter
from collections import defaultdict
from collections import OrderedDict
import statistics 
import re
import gzip
import math




# source_list = ['GO', 'Symptom', 'Disease','Drug', 'Human PPI (target)']
source_list = ['Drug', 'Symptom', 'Human PPI (target)','Disease','GO']



# zscore_dict = {}
# for source in source_list :
# 	node_type = source
# 	with open('COVID19_GDDS_Results_%s_to_CovidProteins_Relative-withZeros.tsv' % source, 'r') as file_read :
# 		data = file_read.readlines()
# 		for line in data[1:] :
# 			line = line.split('\t')
# 			node_name = line[0]
# 			zscore = float(line[-1].replace('\n', ''))
# 			if node_type not in zscore_dict :
# 				zscore_dict[node_type] = []
# 			zscore_dict[node_type].append(zscore)

# df = pd.DataFrame.from_dict(zscore_dict,  orient='index')
# print(len(df))
# # print(df)



deg_dict = {}
for source in source_list :
	node_type = source
	with open('COVID19_GDDS_%s_degrees_to_proteins.tsv' % source, 'r') as file_read :
		data = file_read.readlines()
		for line in data[1:] :
			line = line.split('\t')
			node_name = line[0]
			deg = float(line[1])
			if node_type not in deg_dict :
				deg_dict[node_type] = []
			deg_dict[node_type].append(deg)
df = pd.DataFrame.from_dict(deg_dict,  orient='index')



color_list = ["green", "blue", "orangered", "hotpink", "deepskyblue"]

fig, axes = plt.subplots(nrows=3, ncols=2, sharex=True, sharey=True, figsize=(15, 7))


source_n = 0 

for row in range(0,3) :
	for column in range(0,2) :
		if (row==0 and column==0) :
			for source in range(0,len(source_list)) :
				sns.kdeplot(df.transpose()[source_list[source]], color=color_list[source], shade=False, ax=axes[row][column], legend=False)
		else :
			sns.kdeplot(df.transpose()[source_list[source_n]], color=color_list[source_n], shade=True, legend=True, ax=axes[row][column])
			print(source_list[source_n], row, column)
			source_n = source_n + 1
# plt.savefig('all_zscores_density_final.png')
# plt.close()
fig.suptitle('Absolute degrees to proteins in covid networks - distributions',  fontsize=16)
plt.show()






### plot each plot separately
# for source in source_list :
# 	df[source].plot.kde(title='%s Zscores densities (%s-prot relative degrees)' % (source, source) )
# 	plt.savefig('%s_zscores_density.png' % source)
# 	# plt.show()
# 	plt.close()

# plot all densities in one plot
# df = df.transpose() # this is a trick to not be bothered about the columns not all having the same length.
# df.plot.kde(title='Zscores densities',colormap='Paired')
# plt.savefig('all_zscores_density_paired.png')
# plt.close()














# for j in range(0,len(color_list)) :
# 	sns.kdeplot(df.transpose()[source_list[j]], color=color_list[j], shade=True)
# plt.title('Zscores densities')
# plt.show()


# ncol = 2 # pick one dimension
# nrow = (len(df)+ ncol-1) / ncol # make sure enough subplots
# fig, ax = plt.subplots(nrows=nrow, ncols=ncol) # create the axes

# i = 0
# for source in source_list :
#   ix = np.unravel_index(i, ax.shape) # compute an appropriate index (1d or 2d)
#   df_transposed = df.transpose()
#   df_transposed[source].plot.kde(title='%s Zscores densities (%s-prot relative degrees)' % (source, source), ax=ax[ix] )
#   # ax[ix].plot(...)   # or direct axis object method plot (/scatter/bar/...)
#   plt.savefig('stacked_zscores_density.png')
#   i = i+1

# print(df.index.tolist())
# print(df)
# prop_cycle = plt.rcParams['axes.prop_cycle']
# colors = prop_cycle.by_key()['color']



# df1 = df.loc['GO']
# df1_trans = df1.transpose()
# df1_trans.plot.kde(ax=axes[0], legend=True, color='Paired', shade=True)
# df2 = df.loc['Symptom']
# df2_trans = df2.transpose()
# df2_trans.plot.kde(ax=axes[1], legend=True, color='Paired', shade=True)
# df3 = df.loc['Disease']
# df3_trans = df3.transpose()
# df3_trans.plot.kde(ax=axes[2], legend=True, color='Paired')
# df4 = df.loc['Drug']
# df4_trans = df4.transpose()
# df4_trans.plot.kde(ax=axes[3], legend=True, color='Paired')
# df5 = df.loc['Human PPI (target)']
# df5_trans = df5.transpose()
# df5_trans.plot.kde(ax=axes[4], legend=True, color='Paired')
# plt.show()



# ax = plt.gca()
# df = pd.DataFrame.from_dict(zscore_dict,  orient='index')
# df = df.transpose()
# print(df)
# for source in source_list:
# 	df_plot = df[source]
# 	# print(df_plot)
# 	df_plot.plot.kde(title='%s Zscores densities (%s-prot relative degrees)' % (source, source), ax=ax)
# 	# df_plot.plot.line(ax=ax)
# 	plt.savefig('zscores_density_test.png')




# df.plot.kde()
# df.groupby(['GO','Drug']).mean().unstack().plot()
# plt.show()

# COVID19_GDDS_proteins_degrees_to_Human PPI (target)

# target_list = ['GO', 'Symptom', 'Disease','Drug', 'Human PPI (target)']

# # degrees in covid-net
# deg_dict = {}
# for target in target_list :

# 	with open('full_analysis_for_proteins_undirected\\COVID19_GDDS_proteins_degrees_to_%s.tsv' % target, 'r') as file_read :
# 		data = file_read.readlines()
# 		for line in data :
# 			line = line.split('\t')
# 			protein = line[0]
# 			degree = float(line[1])
# 			degree_rel = float(line[2].replace('\n',''))
# 			if protein not in deg_dict :
# 				deg_dict[protein] = {}
# 			deg_dict[protein][target] = degree
 
# for protein in deg_dict :
# 	print(protein, deg_dict[protein])		

# with open('pulled_results_COVID19_GDDS_proteins_degrees_to_targets_undirected.tsv', 'a+') as file_write :
# 	#heading
# 	file_write.write('protein\t')
# 	for target in target_list[:-1] :
# 		file_write.write('%s\t' % target)
# 	file_write.write('%s\n' % target_list[-1])

# 	for protein in deg_dict :
# 		file_write.write('%s\t' % protein)
# 		for target in target_list[:-1] :
# 			file_write.write('%s\t' % deg_dict[protein][target])
# 		file_write.write('%s\n' % deg_dict[protein][target_list[-1]])









# target_list = ['GO', 'Symptom', 'Disease','Drug', 'Human PPI (target)']

# # degrees in covid-net
# deg_dict = {}
# for target in target_list :

# 	with open('full_analysis_for_proteins_undirected\\COVID19_GDDS_proteins_degrees_to_%s.tsv' % target, 'r') as file_read :
# 		data = file_read.readlines()
# 		for line in data :
# 			line = line.split('\t')
# 			protein = line[0]
# 			degree = float(line[1])
# 			degree_rel = float(line[2].replace('\n',''))
# 			if protein not in deg_dict :
# 				deg_dict[protein] = {}
# 			deg_dict[protein][target] = degree_rel
 
# for protein in deg_dict :
# 	print(protein, deg_dict[protein])		

# with open('pulled_results_COVID19_GDDS_proteins_relative_degrees_to_targets_undirected.tsv', 'a+') as file_write :
# 	#heading
# 	file_write.write('protein\t')
# 	for target in target_list[:-1] :
# 		file_write.write('%s\t' % target)
# 	file_write.write('%s\n' % target_list[-1])

# 	for protein in deg_dict :
# 		file_write.write('%s\t' % protein)
# 		for target in target_list[:-1] :
# 			file_write.write('%s\t' % deg_dict[protein][target])
# 		file_write.write('%s\n' % deg_dict[protein][target_list[-1]])









# target_list = ['GO', 'Symptom', 'Disease','Drug', 'Human PPI (target)']

# # degrees in covid-net
# score_dict = {}
# for target in target_list :

# 	with open('full_analysis_for_proteins_undirected\\COVID19_GDDS_Results_ForCovidProteins_to_%s_Relative-withZeros.tsv' % target, 'r') as file_read :
# 		data = file_read.readlines()
# 		for line in data[1:] :
# 			line = line.split('\t')
# 			protein = line[0]
# 			zscore = float(line[-1].replace('\n',''))
# 			if protein not in score_dict :
# 				score_dict[protein] = {}
# 			score_dict[protein][target] = zscore
 
# for protein in score_dict :
# 	print(protein, score_dict[protein])		

# with open('pulled_results_COVID19_GDDS_proteins_relative_zscores_to_targets_undirected.tsv', 'a+') as file_write :
# 	#heading
# 	file_write.write('protein\t')
# 	for target in target_list[:-1] :
# 		file_write.write('%s\t' % target)
# 	file_write.write('%s\n' % target_list[-1])

# 	for protein in score_dict :
# 		file_write.write('%s\t' % protein)
# 		for target in target_list[:-1] :
# 			file_write.write('%s\t' % score_dict[protein][target])
# 		file_write.write('%s\n' % score_dict[protein][target_list[-1]])









# from scipy.stats import spearmanr
# def corrfunc(x,y, ax=None, **kws):
#     """Plot the correlation coefficient in the top left hand corner of a plot."""
#     r, _ = spearmanr(x, y)
#     ax = ax or plt.gca()
#     # Unicode for lowercase rho (ρ)
#     rho = '\u03C1'
#     ax.annotate(f'{rho}_spearman = {r:.2f}', xy=(.1, .9), xycoords=ax.transAxes)

# # data = pd.read_csv('pulled_results_COVID19_GDDS_proteins_relative_degrees_to_targets_undirected.tsv', usecols=['Drug', 'Human PPI (target)'], sep='\t')
# data = pd.read_csv('pulled_results_COVID19_GDDS_proteins_relative_degrees_to_targets_undirected.tsv', sep='\t')

# # Correlation Matrix Heatmap
# # f, ax = plt.subplots(figsize=(10, 6))
# # corr = data.corr(method='spearman')
# # hm = sns.heatmap(round(corr,2), annot=True, ax=ax, cmap="coolwarm",fmt='.2f',
# #                  linewidths=.05)
# # f.subplots_adjust(top=0.93)
# # t= f.suptitle('Protein Abs Degrees Undirected Correlation Heatmap', fontsize=14)
# # plt.show()


# g = sns.pairplot(data, corner=True, diag_kind="kde", kind="reg", dropna=True)
# # g.map_lower(corrfunc)
# g.fig.suptitle('Protein Relative Degrees Undirected Pairplot', x=0.2, y=1, ha='left')
# plt.show()




# # Pair-wise Scatter Plots
# cols = ['GO', 'Symptom', 'Disease','Drug', 'Human PPI (target)']
# pp = sns.pairplot(data[cols], height=1.8, aspect=1.8,
#                   plot_kws=dict(edgecolor="k", linewidth=0.5),
#                   diag_kind="kde", diag_kws=dict(shade=True))

# fig = pp.fig 
# fig.subplots_adjust(top=0.93, wspace=0.3)
# t = fig.suptitle('Protein Relative Degree Zscores Pairwise Plots', fontsize=14)
# fig.show()


# sns.set(style="ticks", color_codes=True)
# iris = sns.load_dataset("iris")
# print(iris)
# g = sns.pairplot(iris)
# g = sns.pairplot(iris, hue="species")
# g = sns.pairplot(data, corner=True, diag_kind="kde", kind="reg")
# g.fig.suptitle('Protein Abs Degrees Pairplot', x=0.2, y=1, ha='left')
# plt.show()

# sns.set(style="ticks", color_codes=True)
# # g = sns.pairplot(iris)
# # g = sns.pairplot(iris, hue="species")
# # g = sns.pairplot(iris, corner=True, diag_kind="kde", kind="reg")
# g = sns.PairGrid(data)
# g = g.map(plt.scatter)
# plt.show()














# targets = [1,2,3,4,5]
# b = {}
# data = ['prot1','prot2','prot3']
# for target in targets :
# 	for prot in data :
# 		deg = 'deg_%s_%s' % (prot,target)
# 		degrel = 'degrel_%s_%s' % (prot,target)
# 		b[prot] = {target : (deg, degrel)}

# for i in b :
# 	print(i, b[i])







# with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\COVID19_GDDS_nodes.csv', 'r') as nodecsv: # Open the file                       
# 	nodereader = csv.reader(nodecsv) # Read the csv  
# 	nodes = [n for n in nodereader][1:]                     
# 	node_names = [n[0] for n in nodes] # Get a list of only the node names 

# with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\COVID19_GDDS_edges.csv', 'r') as edgecsv: # Open the file
# 	edgereader = csv.reader(edgecsv) # Read the csv     
# 	edges = [tuple(e) for e in edgereader][1:] # Retrieve the data

# G = nx.Graph()
# G.add_edges_from(edges)


# descr_dict = {}
# for node in nodes: 
# 	descr_dict[node[0]] = node[2]

# prots_linked_to_targets = []
# for u,v,c in G.edges(data=True) :	
# 	if (descr_dict[u] == 'Human PPI (target)' and descr_dict[v] == 'Drug') :
# 		prots_linked_to_targets.append(u)
# if 'ACPP' in prots_linked_to_targets :
# 	print('found ACPP')






stop = timeit.default_timer()
print(stop - start)

  
protein-GO-degrees-directed-rankings-corrected.py,
#!/usr/bin/python
# -*- coding: utf-8 -*-
# jupyter notebook --browser='C:/Program Files (x86)/Google/Chrome/Application/chrome.exe %s'


import timeit
start = timeit.default_timer()

import os
import csv
from operator import itemgetter
import networkx as nx
from networkx.algorithms import community
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
from collections import Counter
from collections import defaultdict
from collections import OrderedDict
import statistics 
import re
import gzip





#####################################################################################################
#########################  Compute network Protein degree distribution #################################
##########################################  based on edges to GO #####################################
#####################################################################################################

with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\COVID19_GDDS_nodes.csv', 'r') as nodecsv: # Open the file                       
	nodereader = csv.reader(nodecsv) # Read the csv  
	nodes = [n for n in nodereader][1:]                     
	node_names = [n[0] for n in nodes] # Get a list of only the node names 

with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\COVID19_GDDS_edges.csv', 'r') as edgecsv: # Open the file
	edgereader = csv.reader(edgecsv) # Read the csv     
	edges = [tuple(e) for e in edgereader][1:] # Retrieve the data

G = nx.Graph()
G.add_nodes_from(node_names)
G.add_edges_from(edges)

descr_dict = {}
for node in nodes: 
	descr_dict[node[0]] = node[2]

print(nx.info(G))
# Number of nodes: 16007
# Number of edges: 53917
# Average degree:   6.7367

# print(nx.degree(G)) # this is the list containing as many tuples as nodes, indicating each node's degree.

average_degree = statistics.mean([ tpl[1] for tpl in nx.degree(G) ])
# print(average_degree) # 6.736677703504717

proteins_degrees = nx.degree(G, nbunch=[n for n in G.nodes if descr_dict[n] == 'Human PPI (target)'])
# print("degree of protein nodes only", proteins_degrees)

average_degree_of_proteins = statistics.mean([ tpl[1] for tpl in proteins_degrees ])
print("average degree of proteins nodes", average_degree_of_proteins) # 

# Number of nodes: 16007
# Number of edges: 53917
# Average degree:   6.7367
# average degree of proteins nodes 115.66739606126914


# # # # #####################################################################################################
# # # # ################################# CALCULATION ON EDGES WITH FILTERING PROT-GO #####################
# # # # #####################################################################################################

# now we keep only the edges of proteins that interact with other proteins
prots_linked_to_XXX = []
for u,v,c in G.edges(data=True) :	
	if (descr_dict[u] == 'Human PPI (target)' and descr_dict[v] == 'GO') :
		prots_linked_to_XXX.append(u)
print(len(prots_linked_to_XXX), len(set(prots_linked_to_XXX))) # # 9208 435

# print(sorted(prots_linked_to_XXX))

protein_counts = sorted(prots_linked_to_XXX)

protein_counts_dict = defaultdict( int )
for protein in protein_counts:
    protein_counts_dict[protein] += 1

my_proteins_nodes = [n for n in G.nodes if descr_dict[n] == 'Human PPI (target)']
for protein in my_proteins_nodes :
	if protein_counts_dict.get(protein) == None :
		# print("%s is missing in the linked proteins : Adding it with a frequency of 0.\n" % protein)
		protein_counts_dict[protein]= 0

average_degree_of_my_proteins_nodes = statistics.mean([ protein_counts_dict[protein] for protein in protein_counts_dict ])
print("average_degree_of_my_proteins_nodes", average_degree_of_my_proteins_nodes)

max_degree_of_my_proteins_nodes = max([ protein_counts_dict[protein] for protein in protein_counts_dict ])
print("max_degree_of_my_proteins_nodes", max_degree_of_my_proteins_nodes)

# average_degree_of_my_proteins_nodes 20.14879649890591
# max_degree_of_my_proteins_nodes 141


ordered_protein_counts_dict = OrderedDict(sorted(protein_counts_dict.items()))


# with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\Further_Analysis\\Proteins\\COVID19_GDDS_proteins_degrees_to_GO.tsv', 'a+') as file_write:
# 	for protein in ordered_protein_counts_dict :
# 		file_write.write("%s\t%s\t%s\n" % (protein, ordered_protein_counts_dict[protein], ordered_protein_counts_dict[protein] / sum(ordered_protein_counts_dict.values())))


# # # # ##################################################################
# # # # ############ SAME ANALYSIS FOR RANDOM NETWORKS ###################
# # # # # ##################################################################

# with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\Further_Analysis\\Proteins\\Log_degree_analysis_on_random_networks_for_proteins_GO.txt', 'a+') as log_write:

# 	nodefiles_numbers = []
# 	edgefiles_numbers = []
# 	for filename in os.listdir("A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\random_networks\\Mock_networks_21Apr"):
# 		if 'nodes_' in filename :
# 			filename = filename.split(".csv")[0]
# 			nodefiles_numbers.append(int(filename[19:]))
# 		if 'edges_' in filename :
# 			filename = filename.split(".csv")[0]
# 			edgefiles_numbers.append(int(filename[19:]))

# 	log_write.write("len(nodefiles_numbers) = %s\n" % len(nodefiles_numbers)) #2051
# 	log_write.write("len(edgefiles_numbers) = %s\n" % len(edgefiles_numbers)) #2051
# 	log_write.write("difference between edgefiles_numbers and nodefiles_numbers = %s\n" % len(set(nodefiles_numbers).difference(set(edgefiles_numbers)))) #0

# 	nodefiles_numbers = sorted(nodefiles_numbers)
# 	edgefiles_numbers = sorted(edgefiles_numbers)

# 	networks_to_avoid = [3, 11, 22, 23, 33, 34, 38, 57, 72, 77, 123, 127, 148, 156, 157, 183, 205, 209, 210, 217, 232, 252, 268, 281, 291, 312, 341, 353, 362, 364, 413, 421, 433, 443, 451, 494, 525, 532, 569, 592, 629, 658, 662, 674, 676, 701, 705, 712, 736, 757, 788, 796, 825, 827, 838, 873, 876, 908, 910, 916, 917, 918, 938, 958, 965, 1088, 1125, 1128, 1145, 1146, 1165, 1171, 1183, 1185, 1192, 1195, 1211, 1244, 1254, 1257, 1300, 1353, 1354, 1364, 1372, 1394, 1395, 1401, 1410, 1427, 1454, 1477, 1497, 1517, 1532, 1549, 1551, 1575, 1578, 1611, 1624, 1629, 1637, 1650, 1655, 1697, 1700, 1713, 1736, 1745, 1762, 1790, 1800, 1815, 1824, 1827, 1866, 1872, 1886, 1887, 2053, 2072, 2106, 2115, 2128, 2160, 2217, 2343, 2383, 2388, 2413, 2458, 2465]
# 	nodefiles_numbers = [item for item in nodefiles_numbers if item not in networks_to_avoid]

# 	for iteration in nodefiles_numbers[1:] :
# 		print("\n########################## ITERATION %s ###########################\n" % iteration)
# 		log_write.write("\n########################## ITERATION %s ###########################\n" % iteration)
	
# 		with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\random_networks\\Mock_networks_21Apr\\COVID19_GDDS_nodes_%s.csv' % str(iteration), 'r') as nodecsv: # Open the file                   
# 			nodereader = csv.reader(nodecsv) # Read the csv  
# 			nodes = [n for n in nodereader][1:]                     
# 			node_names = [n[0] for n in nodes] # Get a list of only the node names 

# 		with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\random_networks\\Mock_networks_21Apr\\COVID19_GDDS_edges_%s.csv' % str(iteration), 'r') as edgecsv: # Open the file
# 			edgereader = csv.reader(edgecsv) # Read the csv     
# 			edges = [tuple(e) for e in edgereader][1:] # Retrieve the data

# 			log_write.write("len(node_names) : %s\n" % len(node_names)) #16007 --> 13448
# 			log_write.write("len(edges) : %s\n" % len(edges)) #53917 --> 37835

# 		G = nx.Graph()
# 		G.add_nodes_from(node_names)
# 		G.add_edges_from(edges)

# 		descr_dict = {}
# 		description_set=set()
# 		for node in nodes: 
# 			descr_dict[node[0]] = node[2]
# 			description_set.add(node[2])

# 		# print(nx.info(G))
# 		log_write.write("\nRANDOM GRAPH %s\n" % iteration)
# 		log_write.write("%s\n" % nx.info(G))

# 		#####################################################################################################
# 		################################# GET EDGES CONNECTING proteinS TO GO ################
# 		#####################################################################################################

# 		prots_linked_to_XXX = []
# 		for u,v,c in G.edges(data=True) :	
# 			if (descr_dict[u] == 'Human PPI (target)' and descr_dict[v] == 'GO') :
# 				prots_linked_to_XXX.append(u)
		
# 		protein_counts = sorted(prots_linked_to_XXX)

# 		protein_counts_dict = defaultdict( int )
# 		for protein in protein_counts:
# 		    protein_counts_dict[protein] += 1

# 		my_proteins_nodes = [n for n in G.nodes if descr_dict[n] == 'Human PPI (target)']
# 		for protein in my_proteins_nodes :
# 			if protein_counts_dict.get(protein) == None :
# 				# print("%s is missing in the linked proteins : Adding it with a frequency of 0.\n" % protein)
# 				protein_counts_dict[protein]= 0

# 		average_degree_of_my_proteins_nodes = statistics.mean([ protein_counts_dict[protein] for protein in protein_counts_dict ])
# 		log_write.write("Average degree of proteins linked to GO(s) : %s\n" % average_degree_of_my_proteins_nodes)

# 		ordered_protein_counts_dict = OrderedDict(sorted(protein_counts_dict.items()))

# 		with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\Further_Analysis\\Proteins\\prot-GO\\COVID19_GDDS_proteins_degrees_to_GO_%s.tsv' % str(iteration), 'a+') as file_write:
# 			for protein in ordered_protein_counts_dict :
# 				if (sum(ordered_protein_counts_dict.values()) != 0) :
# 					file_write.write("%s\t%s\t%s\n" % (protein, ordered_protein_counts_dict[protein], ordered_protein_counts_dict[protein] / sum(ordered_protein_counts_dict.values())))
# 				else :
# 					file_write.write("%s\t%s\t%s\n" % (protein, ordered_protein_counts_dict[protein], ordered_protein_counts_dict[protein]))



# # # # ##################################################################
# # # # ############ COMPUTE DEGREES MEANS AND Z-SCORES ##################
# # # # ######## INCLUDING ABSENCE OF GO AS ZERO DEGREES ###########
# # # # ##################################################################


# # For covid network, compute the mean and the standard deviation of proteins degrees and normalized degrees
covid_proteins = []
covid_protein_degree_dict = {}
covid_protein_norm_degree_dict = {}

total_protein_degree_dict = defaultdict(list)
total_protein_normalized_degree_dict = defaultdict(list)

with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\Further_Analysis\\Proteins\\COVID19_GDDS_proteins_degrees_to_GO.tsv', 'r') as data :
	for line in data :
		line = line.split('\t')
		protein = line[0]
		covid_proteins.append(protein)
		protein_degree = int(line[1])
		protein_norm_degree = float(line[2].replace('\n', ''))
		if covid_protein_degree_dict.get(protein) == None :
			covid_protein_degree_dict[protein] = protein_degree
		else :
			print("warning : %s duplicate ?!" % protein)
		covid_protein_norm_degree_dict[protein] = protein_norm_degree

		if total_protein_degree_dict.get(protein) == None :
			total_protein_degree_dict[protein] = [protein_degree]
		else :
			total_protein_degree_dict[protein].append(protein_degree)

		if total_protein_normalized_degree_dict.get(protein) == None :
			total_protein_normalized_degree_dict[protein] = [protein_norm_degree]
		else :
			total_protein_normalized_degree_dict[protein].append(protein_norm_degree)

mean_protein_deg = statistics.mean(covid_protein_degree_dict.values())
mean_norm_protein_deg = statistics.mean(covid_protein_norm_degree_dict.values())
sd_protein_deg = statistics.pstdev(covid_protein_degree_dict.values())
sd_norm_protein_deg = statistics.pstdev(covid_protein_norm_degree_dict.values())
print(mean_protein_deg, sd_protein_deg, mean_norm_protein_deg, sd_norm_protein_deg) #
print(len(covid_proteins),len(covid_protein_degree_dict),len(covid_protein_norm_degree_dict)) #
# 20.14879649890591 18.9696222570052 0.002188183807439825 0.00206012405050013
# 457 457 457






nodefiles_numbers = []
edgefiles_numbers = []
for filename in os.listdir("A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\random_networks\\Mock_networks_21Apr"):
	if 'nodes_' in filename :
		filename = filename.split(".csv")[0]
		nodefiles_numbers.append(int(filename[19:]))
	if 'edges_' in filename :
		filename = filename.split(".csv")[0]
		edgefiles_numbers.append(int(filename[19:]))

# log_write.write("len(nodefiles_numbers) = %s\n" % len(nodefiles_numbers)) #2051
# log_write.write("len(edgefiles_numbers) = %s\n" % len(edgefiles_numbers)) #2051
# log_write.write("difference between edgefiles_numbers and nodefiles_numbers = %s\n" % len(set(nodefiles_numbers).difference(set(edgefiles_numbers)))) #0
print("len(nodefiles_numbers) = %s\n" % len(nodefiles_numbers)) #2051
print("len(edgefiles_numbers) = %s\n" % len(edgefiles_numbers)) #2051
print("difference between edgefiles_numbers and nodefiles_numbers = %s\n" % len(set(nodefiles_numbers).difference(set(edgefiles_numbers)))) #0

# nodefiles_numbers = sorted(nodefiles_numbers)
# edgefiles_numbers = sorted(edgefiles_numbers)

networks_to_avoid = [3, 11, 22, 23, 33, 34, 38, 57, 72, 77, 123, 127, 148, 156, 157, 183, 205, 209, 210, 217, 232, 252, 268, 281, 291, 312, 341, 353, 362, 364, 413, 421, 433, 443, 451, 494, 525, 532, 569, 592, 629, 658, 662, 674, 676, 701, 705, 712, 736, 757, 788, 796, 825, 827, 838, 873, 876, 908, 910, 916, 917, 918, 938, 958, 965, 1088, 1125, 1128, 1145, 1146, 1165, 1171, 1183, 1185, 1192, 1195, 1211, 1244, 1254, 1257, 1300, 1353, 1354, 1364, 1372, 1394, 1395, 1401, 1410, 1427, 1454, 1477, 1497, 1517, 1532, 1549, 1551, 1575, 1578, 1611, 1624, 1629, 1637, 1650, 1655, 1697, 1700, 1713, 1736, 1745, 1762, 1790, 1800, 1815, 1824, 1827, 1866, 1872, 1886, 1887, 2053, 2072, 2106, 2115, 2128, 2160, 2217, 2343, 2383, 2388, 2413, 2458, 2465]
nodefiles_numbers = [item for item in nodefiles_numbers if item not in networks_to_avoid]

for iteration in nodefiles_numbers[1:] :
	# log_write.write("\n########################## ITERATION %s ###########################\n" % iteration)
	print("\n########################## ITERATION %s ###########################\n" % iteration)

	with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\Further_Analysis\\Proteins\\prot-GO\\COVID19_GDDS_proteins_degrees_to_GO_%s.tsv' % str(iteration), 'r') as data: # Open the file  
		for line in data :
			line = line.split('\t')
			protein = line[0]
			degree = float(line[1])
			normalized_degree = float(line[2].replace('\n', ''))

			if total_protein_degree_dict.get(protein) == None :
				total_protein_degree_dict[protein] = [degree]
			else :
				total_protein_degree_dict[protein].append(degree)

			if total_protein_normalized_degree_dict.get(protein) == None :
				total_protein_normalized_degree_dict[protein] = [normalized_degree]
			else :
				total_protein_normalized_degree_dict[protein].append(normalized_degree)

print(len(total_protein_degree_dict)) # 19928

print(len(total_protein_normalized_degree_dict)) # 19928

### make the files with all calculations
with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\Further_Analysis\\Proteins\\COVID19_GDDS_Results_Proteins_to_GO-withZeros.tsv', 'a+') as file_write:
	file_write.write("protein\t# of nets it is found in\tmin-degree\tmax-degree\tmean-degree\tmedian-degree\tsd-degree\n")
	for protein in total_protein_degree_dict :
		nets_number = len(total_protein_degree_dict[protein])
		min_deg = min(total_protein_degree_dict[protein])
		max_deg = max(total_protein_degree_dict[protein])
		mean_deg = statistics.mean(total_protein_degree_dict[protein])
		median_deg = statistics.median(total_protein_degree_dict[protein])
		sd_deg = statistics.pstdev(total_protein_degree_dict[protein])
		file_write.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (protein, nets_number, min_deg, max_deg, mean_deg, median_deg, sd_deg))


with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\Further_Analysis\\Proteins\\COVID19_GDDS_Results_Proteins_to_GO_Relative-withZeros.tsv', 'a+') as file_write:
	file_write.write("protein\t# of nets it is found in\tmin-degree-rel\tmax-degree-rel\tmean-degree-rel\tmedian-degree-rel\tsd-degree-rel\n")
	for protein in total_protein_normalized_degree_dict :
		nets_number = len(total_protein_normalized_degree_dict[protein])
		min_deg = min(total_protein_normalized_degree_dict[protein])
		max_deg = max(total_protein_normalized_degree_dict[protein])
		mean_deg = statistics.mean(total_protein_normalized_degree_dict[protein])
		median_deg = statistics.median(total_protein_normalized_degree_dict[protein])
		sd_deg = statistics.pstdev(total_protein_normalized_degree_dict[protein])
		file_write.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (protein, nets_number, min_deg, max_deg, mean_deg, median_deg, sd_deg))

### compute Z-score for the covid proteins
# For each protein found in covid network, compute 2 Z-scores.

with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\Further_Analysis\\Proteins\\COVID19_GDDS_Results_ForCovidProteins_to_GO-withZeros.tsv', 'a+') as file_write:
	file_write.write("covid protein\t# of nets it is found in\tmin-degree\tmax-degree\tmean-degree\tmedian-degree\tsd-degree\tcovid-degree\tZ-score\n")
	for protein in covid_proteins :
		nets_number = len(total_protein_degree_dict[protein])
		min_deg = min(total_protein_degree_dict[protein])
		max_deg = max(total_protein_degree_dict[protein])
		mean_deg = statistics.mean(total_protein_degree_dict[protein])
		median_deg = statistics.median(total_protein_degree_dict[protein])
		sd_deg = statistics.pstdev(total_protein_degree_dict[protein])
		covid_deg = covid_protein_degree_dict[protein]
		if sd_deg != 0 :
			z_score = ( covid_deg - mean_deg) / sd_deg
		else : 
			z_score = 'NaN'
		file_write.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (protein, nets_number, min_deg, max_deg, mean_deg, median_deg, sd_deg, covid_deg, z_score))


with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\Further_Analysis\\Proteins\\COVID19_GDDS_Results_ForCovidProteins_to_GO_Relative-withZeros.tsv', 'a+') as file_write:
	file_write.write("covid protein\t# of nets it is found in\tmin-degree-rel\tmax-degree-rel\tmean-degree-rel\tmedian-degree-rel\tsd-degree-rel\tcovid-degree-rel\tZ-score-rel\n")
	for protein in covid_proteins :
		nets_number = len(total_protein_normalized_degree_dict[protein])
		min_deg = min(total_protein_normalized_degree_dict[protein])
		max_deg = max(total_protein_normalized_degree_dict[protein])
		mean_deg = statistics.mean(total_protein_normalized_degree_dict[protein])
		median_deg = statistics.median(total_protein_normalized_degree_dict[protein])
		sd_deg = statistics.pstdev(total_protein_normalized_degree_dict[protein])
		covid_deg = covid_protein_norm_degree_dict[protein]
		if sd_deg != 0 :
			z_score = ( covid_deg - mean_deg) / sd_deg
		else : 
			z_score = 'NaN'
		file_write.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (protein, nets_number, min_deg, max_deg, mean_deg, median_deg, sd_deg, covid_deg, z_score))










stop = timeit.default_timer()
print(stop - start)
  
protein-disease-degrees-directed-rankings-corrected.py,
#!/usr/bin/python
# -*- coding: utf-8 -*-
# jupyter notebook --browser='C:/Program Files (x86)/Google/Chrome/Application/chrome.exe %s'


import timeit
start = timeit.default_timer()

import os
import csv
from operator import itemgetter
import networkx as nx
from networkx.algorithms import community
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
from collections import Counter
from collections import defaultdict
from collections import OrderedDict
import statistics 
import re
import gzip





#####################################################################################################
#########################  Compute network Protein degree distribution #################################
##########################################  based on edges to diseases #####################################
#####################################################################################################

with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\COVID19_GDDS_nodes.csv', 'r') as nodecsv: # Open the file                       
	nodereader = csv.reader(nodecsv) # Read the csv  
	nodes = [n for n in nodereader][1:]                     
	node_names = [n[0] for n in nodes] # Get a list of only the node names 

with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\COVID19_GDDS_edges.csv', 'r') as edgecsv: # Open the file
	edgereader = csv.reader(edgecsv) # Read the csv     
	edges = [tuple(e) for e in edgereader][1:] # Retrieve the data

G = nx.Graph()
G.add_nodes_from(node_names)
G.add_edges_from(edges)

descr_dict = {}
for node in nodes: 
	descr_dict[node[0]] = node[2]

print(nx.info(G))
# Number of nodes: 16007
# Number of edges: 53917
# Average degree:   6.7367

# print(nx.degree(G)) # this is the list containing as many tuples as nodes, indicating each node's degree.

average_degree = statistics.mean([ tpl[1] for tpl in nx.degree(G) ])
# print(average_degree) # 6.736677703504717

proteins_degrees = nx.degree(G, nbunch=[n for n in G.nodes if descr_dict[n] == 'Human PPI (target)'])
# print("degree of protein nodes only", proteins_degrees)

average_degree_of_proteins = statistics.mean([ tpl[1] for tpl in proteins_degrees ])
print("average degree of proteins nodes", average_degree_of_proteins) # 

# Number of nodes: 16007
# Number of edges: 53917
# Average degree:   6.7367
# average degree of proteins nodes 115.66739606126914


# # # # #####################################################################################################
# # # # ################################# CALCULATION ON EDGES WITH FILTERING PROT-disease #####################
# # # # #####################################################################################################

# now we keep only the edges of proteins that interact with other proteins
prots_linked_to_XXX = []
for u,v,c in G.edges(data=True) :	
	if (descr_dict[u] == 'Human PPI (target)' and descr_dict[v] == 'Disease') :
		prots_linked_to_XXX.append(u)
print(len(prots_linked_to_XXX), len(set(prots_linked_to_XXX))) # 16124 387

# ({'Drug': 5703, 'Disease': 4176, 'GO': 3487, 'Symptom': 2157, 'Human PPI (target)': 457, 'Viral Gene': 27})

# print(sorted(prots_linked_to_XXX))

protein_counts = sorted(prots_linked_to_XXX)

protein_counts_dict = defaultdict( int )
for protein in protein_counts:
    protein_counts_dict[protein] += 1

my_proteins_nodes = [n for n in G.nodes if descr_dict[n] == 'Human PPI (target)']
for protein in my_proteins_nodes :
	if protein_counts_dict.get(protein) == None :
		# print("%s is missing in the linked proteins : Adding it with a frequency of 0.\n" % protein)
		protein_counts_dict[protein]= 0

average_degree_of_my_proteins_nodes = statistics.mean([ protein_counts_dict[protein] for protein in protein_counts_dict ])
print("average_degree_of_my_proteins_nodes", average_degree_of_my_proteins_nodes)

max_degree_of_my_proteins_nodes = max([ protein_counts_dict[protein] for protein in protein_counts_dict ])
print("max_degree_of_my_proteins_nodes", max_degree_of_my_proteins_nodes)

# average_degree_of_my_proteins_nodes 35.282275711159734
# max_degree_of_my_proteins_nodes 844




ordered_protein_counts_dict = OrderedDict(sorted(protein_counts_dict.items()))


# with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\Further_Analysis\\Proteins\\COVID19_GDDS_proteins_degrees_to_disease.tsv', 'a+') as file_write:
# 	for protein in ordered_protein_counts_dict :
# 		file_write.write("%s\t%s\t%s\n" % (protein, ordered_protein_counts_dict[protein], ordered_protein_counts_dict[protein] / sum(ordered_protein_counts_dict.values())))


# # # ##################################################################
# # # ############ SAME ANALYSIS FOR RANDOM NETWORKS ###################
# # # # ##################################################################

# with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\Further_Analysis\\Proteins\\Log_degree_analysis_on_random_networks_for_proteins_diseases.txt', 'a+') as log_write:

# 	nodefiles_numbers = []
# 	edgefiles_numbers = []
# 	for filename in os.listdir("A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\random_networks\\Mock_networks_21Apr"):
# 		if 'nodes_' in filename :
# 			filename = filename.split(".csv")[0]
# 			nodefiles_numbers.append(int(filename[19:]))
# 		if 'edges_' in filename :
# 			filename = filename.split(".csv")[0]
# 			edgefiles_numbers.append(int(filename[19:]))

# 	log_write.write("len(nodefiles_numbers) = %s\n" % len(nodefiles_numbers)) #2051
# 	log_write.write("len(edgefiles_numbers) = %s\n" % len(edgefiles_numbers)) #2051
# 	log_write.write("difference between edgefiles_numbers and nodefiles_numbers = %s\n" % len(set(nodefiles_numbers).difference(set(edgefiles_numbers)))) #0

# 	nodefiles_numbers = sorted(nodefiles_numbers)
# 	edgefiles_numbers = sorted(edgefiles_numbers)

# 	networks_to_avoid = [3, 11, 22, 23, 33, 34, 38, 57, 72, 77, 123, 127, 148, 156, 157, 183, 205, 209, 210, 217, 232, 252, 268, 281, 291, 312, 341, 353, 362, 364, 413, 421, 433, 443, 451, 494, 525, 532, 569, 592, 629, 658, 662, 674, 676, 701, 705, 712, 736, 757, 788, 796, 825, 827, 838, 873, 876, 908, 910, 916, 917, 918, 938, 958, 965, 1088, 1125, 1128, 1145, 1146, 1165, 1171, 1183, 1185, 1192, 1195, 1211, 1244, 1254, 1257, 1300, 1353, 1354, 1364, 1372, 1394, 1395, 1401, 1410, 1427, 1454, 1477, 1497, 1517, 1532, 1549, 1551, 1575, 1578, 1611, 1624, 1629, 1637, 1650, 1655, 1697, 1700, 1713, 1736, 1745, 1762, 1790, 1800, 1815, 1824, 1827, 1866, 1872, 1886, 1887, 2053, 2072, 2106, 2115, 2128, 2160, 2217, 2343, 2383, 2388, 2413, 2458, 2465]
# 	nodefiles_numbers = [item for item in nodefiles_numbers if item not in networks_to_avoid]

# 	for iteration in nodefiles_numbers[1:] :
# 		print("\n########################## ITERATION %s ###########################\n" % iteration)
# 		log_write.write("\n########################## ITERATION %s ###########################\n" % iteration)
	
# 		with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\random_networks\\Mock_networks_21Apr\\COVID19_GDDS_nodes_%s.csv' % str(iteration), 'r') as nodecsv: # Open the file                   
# 			nodereader = csv.reader(nodecsv) # Read the csv  
# 			nodes = [n for n in nodereader][1:]                     
# 			node_names = [n[0] for n in nodes] # Get a list of only the node names 

# 		with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\random_networks\\Mock_networks_21Apr\\COVID19_GDDS_edges_%s.csv' % str(iteration), 'r') as edgecsv: # Open the file
# 			edgereader = csv.reader(edgecsv) # Read the csv     
# 			edges = [tuple(e) for e in edgereader][1:] # Retrieve the data

# 			log_write.write("len(node_names) : %s\n" % len(node_names)) #16007 --> 13448
# 			log_write.write("len(edges) : %s\n" % len(edges)) #53917 --> 37835

# 		G = nx.Graph()
# 		G.add_nodes_from(node_names)
# 		G.add_edges_from(edges)

# 		descr_dict = {}
# 		description_set=set()
# 		for node in nodes: 
# 			descr_dict[node[0]] = node[2]
# 			description_set.add(node[2])

# 		# print(description_set) #{'Symptom', 'Human PPI (target)', 'GO', 'Disease', 'Drug'}
# 		# recounted = Counter(list(descr_dict.values()))
# 		# print("%s\n" % recounted)

# 		# print(nx.info(G))
# 		log_write.write("\nRANDOM GRAPH %s\n" % iteration)
# 		log_write.write("%s\n" % nx.info(G))

# 		#####################################################################################################
# 		################################# GET EDGES CONNECTING proteinS TO FIRST ORDER PROTEINS ################
# 		#####################################################################################################

# 		prots_linked_to_XXX = []
# 		for u,v,c in G.edges(data=True) :	
# 			if (descr_dict[u] == 'Human PPI (target)' and descr_dict[v] == 'Disease') :
# 				prots_linked_to_XXX.append(u)
		
# 		protein_counts = sorted(prots_linked_to_XXX)

# 		protein_counts_dict = defaultdict( int )
# 		for protein in protein_counts:
# 		    protein_counts_dict[protein] += 1

# 		my_proteins_nodes = [n for n in G.nodes if descr_dict[n] == 'Human PPI (target)']
# 		for protein in my_proteins_nodes :
# 			if protein_counts_dict.get(protein) == None :
# 				# print("%s is missing in the linked proteins : Adding it with a frequency of 0.\n" % protein)
# 				protein_counts_dict[protein]= 0

# 		average_degree_of_my_proteins_nodes = statistics.mean([ protein_counts_dict[protein] for protein in protein_counts_dict ])
# 		log_write.write("Average degree of proteins linked to disease(s) : %s\n" % average_degree_of_my_proteins_nodes)

# 		ordered_protein_counts_dict = OrderedDict(sorted(protein_counts_dict.items()))

# 		with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\Further_Analysis\\Proteins\\prot-dis\\COVID19_GDDS_proteins_degrees_to_diseases_%s.tsv' % str(iteration), 'a+') as file_write:
# 			for protein in ordered_protein_counts_dict :
# 				if (sum(ordered_protein_counts_dict.values()) != 0) :
# 					file_write.write("%s\t%s\t%s\n" % (protein, ordered_protein_counts_dict[protein], ordered_protein_counts_dict[protein] / sum(ordered_protein_counts_dict.values())))
# 				else :
# 					file_write.write("%s\t%s\t%s\n" % (protein, ordered_protein_counts_dict[protein], ordered_protein_counts_dict[protein]))



# # # # ##################################################################
# # # # ############ COMPUTE DEGREES MEANS AND Z-SCORES ##################
# # # # ######## INCLUDING ABSENCE OF diseases AS ZERO DEGREES ###########
# # # # ##################################################################


# # For covid network, compute the mean and the standard deviation of proteins degrees and normalized degrees
covid_proteins = []
covid_protein_degree_dict = {}
covid_protein_norm_degree_dict = {}

total_protein_degree_dict = defaultdict(list)
total_protein_normalized_degree_dict = defaultdict(list)

with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\Further_Analysis\\Proteins\\COVID19_GDDS_proteins_degrees_to_disease.tsv', 'r') as data :
	for line in data :
		line = line.split('\t')
		protein = line[0]
		covid_proteins.append(protein)
		protein_degree = int(line[1])
		protein_norm_degree = float(line[2].replace('\n', ''))
		if covid_protein_degree_dict.get(protein) == None :
			covid_protein_degree_dict[protein] = protein_degree
		else :
			print("warning : %s duplicate ?!" % protein)
		covid_protein_norm_degree_dict[protein] = protein_norm_degree

		if total_protein_degree_dict.get(protein) == None :
			total_protein_degree_dict[protein] = [protein_degree]
		else :
			total_protein_degree_dict[protein].append(protein_degree)

		if total_protein_normalized_degree_dict.get(protein) == None :
			total_protein_normalized_degree_dict[protein] = [protein_norm_degree]
		else :
			total_protein_normalized_degree_dict[protein].append(protein_norm_degree)

mean_protein_deg = statistics.mean(covid_protein_degree_dict.values())
mean_norm_protein_deg = statistics.mean(covid_protein_norm_degree_dict.values())
sd_protein_deg = statistics.pstdev(covid_protein_degree_dict.values())
sd_norm_protein_deg = statistics.pstdev(covid_protein_norm_degree_dict.values())
print(mean_protein_deg, sd_protein_deg, mean_norm_protein_deg, sd_norm_protein_deg) #
print(len(covid_proteins),len(covid_protein_degree_dict),len(covid_protein_norm_degree_dict)) #
# 35.282275711159734 70.47315218003904 0.002188183807439825 0.004370699093279524
# 457 457 457




nodefiles_numbers = []
edgefiles_numbers = []
for filename in os.listdir("A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\random_networks\\Mock_networks_21Apr"):
	if 'nodes_' in filename :
		filename = filename.split(".csv")[0]
		nodefiles_numbers.append(int(filename[19:]))
	if 'edges_' in filename :
		filename = filename.split(".csv")[0]
		edgefiles_numbers.append(int(filename[19:]))

# log_write.write("len(nodefiles_numbers) = %s\n" % len(nodefiles_numbers)) #2051
# log_write.write("len(edgefiles_numbers) = %s\n" % len(edgefiles_numbers)) #2051
# log_write.write("difference between edgefiles_numbers and nodefiles_numbers = %s\n" % len(set(nodefiles_numbers).difference(set(edgefiles_numbers)))) #0
print("len(nodefiles_numbers) = %s\n" % len(nodefiles_numbers)) #2051
print("len(edgefiles_numbers) = %s\n" % len(edgefiles_numbers)) #2051
print("difference between edgefiles_numbers and nodefiles_numbers = %s\n" % len(set(nodefiles_numbers).difference(set(edgefiles_numbers)))) #0

# nodefiles_numbers = sorted(nodefiles_numbers)
# edgefiles_numbers = sorted(edgefiles_numbers)

networks_to_avoid = [3, 11, 22, 23, 33, 34, 38, 57, 72, 77, 123, 127, 148, 156, 157, 183, 205, 209, 210, 217, 232, 252, 268, 281, 291, 312, 341, 353, 362, 364, 413, 421, 433, 443, 451, 494, 525, 532, 569, 592, 629, 658, 662, 674, 676, 701, 705, 712, 736, 757, 788, 796, 825, 827, 838, 873, 876, 908, 910, 916, 917, 918, 938, 958, 965, 1088, 1125, 1128, 1145, 1146, 1165, 1171, 1183, 1185, 1192, 1195, 1211, 1244, 1254, 1257, 1300, 1353, 1354, 1364, 1372, 1394, 1395, 1401, 1410, 1427, 1454, 1477, 1497, 1517, 1532, 1549, 1551, 1575, 1578, 1611, 1624, 1629, 1637, 1650, 1655, 1697, 1700, 1713, 1736, 1745, 1762, 1790, 1800, 1815, 1824, 1827, 1866, 1872, 1886, 1887, 2053, 2072, 2106, 2115, 2128, 2160, 2217, 2343, 2383, 2388, 2413, 2458, 2465]
nodefiles_numbers = [item for item in nodefiles_numbers if item not in networks_to_avoid]

# for iteration in nodefiles_numbers[1:] :
# 	# log_write.write("\n########################## ITERATION %s ###########################\n" % iteration)
# 	print("\n########################## ITERATION %s ###########################\n" % iteration)

# 	with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\Further_Analysis\\Proteins\\prot-dis\\COVID19_GDDS_proteins_degrees_to_diseases_%s.tsv' % str(iteration), 'r') as data: # Open the file  
# 		for line in data :
# 			line = line.split('\t')
# 			protein = line[0]
# 			degree = float(line[1])
# 			normalized_degree = float(line[2].replace('\n', ''))

# 			if total_protein_degree_dict.get(protein) == None :
# 				total_protein_degree_dict[protein] = [degree]
# 			else :
# 				total_protein_degree_dict[protein].append(degree)

# 			if total_protein_normalized_degree_dict.get(protein) == None :
# 				total_protein_normalized_degree_dict[protein] = [normalized_degree]
# 			else :
# 				total_protein_normalized_degree_dict[protein].append(normalized_degree)

# print(len(total_protein_degree_dict)) # 19928

# print(len(total_protein_normalized_degree_dict)) # 19928

# ### make the files with all calculations
# with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\Further_Analysis\\Proteins\\COVID19_GDDS_Results_Proteins_to_diseases-withZeros.tsv', 'a+') as file_write:
# 	file_write.write("protein\t# of nets it is found in\tmin-degree\tmax-degree\tmean-degree\tmedian-degree\tsd-degree\n")
# 	for protein in total_protein_degree_dict :
# 		nets_number = len(total_protein_degree_dict[protein])
# 		min_deg = min(total_protein_degree_dict[protein])
# 		max_deg = max(total_protein_degree_dict[protein])
# 		mean_deg = statistics.mean(total_protein_degree_dict[protein])
# 		median_deg = statistics.median(total_protein_degree_dict[protein])
# 		sd_deg = statistics.pstdev(total_protein_degree_dict[protein])
# 		file_write.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (protein, nets_number, min_deg, max_deg, mean_deg, median_deg, sd_deg))


# with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\Further_Analysis\\Proteins\\COVID19_GDDS_Results_Proteins_to_diseases_Relative-withZeros.tsv', 'a+') as file_write:
# 	file_write.write("protein\t# of nets it is found in\tmin-degree-rel\tmax-degree-rel\tmean-degree-rel\tmedian-degree-rel\tsd-degree-rel\n")
# 	for protein in total_protein_normalized_degree_dict :
# 		nets_number = len(total_protein_normalized_degree_dict[protein])
# 		min_deg = min(total_protein_normalized_degree_dict[protein])
# 		max_deg = max(total_protein_normalized_degree_dict[protein])
# 		mean_deg = statistics.mean(total_protein_normalized_degree_dict[protein])
# 		median_deg = statistics.median(total_protein_normalized_degree_dict[protein])
# 		sd_deg = statistics.pstdev(total_protein_normalized_degree_dict[protein])
# 		file_write.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (protein, nets_number, min_deg, max_deg, mean_deg, median_deg, sd_deg))

# ### compute Z-score for the covid proteins
# # For each protein found in covid network, compute 2 Z-scores.

# with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\Further_Analysis\\Proteins\\COVID19_GDDS_Results_ForCovidProteins_to_diseases-withZeros.tsv', 'a+') as file_write:
# 	file_write.write("covid protein\t# of nets it is found in\tmin-degree\tmax-degree\tmean-degree\tmedian-degree\tsd-degree\tcovid-degree\tZ-score\n")
# 	for protein in covid_proteins :
# 		nets_number = len(total_protein_degree_dict[protein])
# 		min_deg = min(total_protein_degree_dict[protein])
# 		max_deg = max(total_protein_degree_dict[protein])
# 		mean_deg = statistics.mean(total_protein_degree_dict[protein])
# 		median_deg = statistics.median(total_protein_degree_dict[protein])
# 		sd_deg = statistics.pstdev(total_protein_degree_dict[protein])
# 		covid_deg = covid_protein_degree_dict[protein]
# 		if sd_deg != 0 :
# 			z_score = ( covid_deg - mean_deg) / sd_deg
# 		else : 
# 			z_score = 'NaN'
# 		file_write.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (protein, nets_number, min_deg, max_deg, mean_deg, median_deg, sd_deg, covid_deg, z_score))


# with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\Further_Analysis\\Proteins\\COVID19_GDDS_Results_ForCovidProteins_to_diseases_Relative-withZeros.tsv', 'a+') as file_write:
# 	file_write.write("covid protein\t# of nets it is found in\tmin-degree-rel\tmax-degree-rel\tmean-degree-rel\tmedian-degree-rel\tsd-degree-rel\tcovid-degree-rel\tZ-score-rel\n")
# 	for protein in covid_proteins :
# 		nets_number = len(total_protein_normalized_degree_dict[protein])
# 		min_deg = min(total_protein_normalized_degree_dict[protein])
# 		max_deg = max(total_protein_normalized_degree_dict[protein])
# 		mean_deg = statistics.mean(total_protein_normalized_degree_dict[protein])
# 		median_deg = statistics.median(total_protein_normalized_degree_dict[protein])
# 		sd_deg = statistics.pstdev(total_protein_normalized_degree_dict[protein])
# 		covid_deg = covid_protein_norm_degree_dict[protein]
# 		if sd_deg != 0 :
# 			z_score = ( covid_deg - mean_deg) / sd_deg
# 		else : 
# 			z_score = 'NaN'
# 		file_write.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (protein, nets_number, min_deg, max_deg, mean_deg, median_deg, sd_deg, covid_deg, z_score))










stop = timeit.default_timer()
print(stop - start)
  
protein-drug-degrees-directed-rankings.py,
#!/usr/bin/python
# -*- coding: utf-8 -*-
# jupyter notebook --browser='C:/Program Files (x86)/Google/Chrome/Application/chrome.exe %s'


import timeit
start = timeit.default_timer()

import os
import csv
from operator import itemgetter
import networkx as nx
from networkx.algorithms import community
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
from collections import Counter
from collections import defaultdict
from collections import OrderedDict
import statistics 
import re
import gzip





#####################################################################################################
#########################  Compute network Protein degree distribution #################################
##########################################  based on edges to proteins #####################################
#####################################################################################################

with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\COVID19_GDDS_nodes.csv', 'r') as nodecsv: # Open the file                       
	nodereader = csv.reader(nodecsv) # Read the csv  
	nodes = [n for n in nodereader][1:]                     
	node_names = [n[0] for n in nodes] # Get a list of only the node names 

with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\COVID19_GDDS_edges.csv', 'r') as edgecsv: # Open the file
	edgereader = csv.reader(edgecsv) # Read the csv     
	edges = [tuple(e) for e in edgereader][1:] # Retrieve the data

G = nx.Graph()
G.add_nodes_from(node_names)
G.add_edges_from(edges)

descr_dict = {}
for node in nodes: 
	descr_dict[node[0]] = node[2]

print(nx.info(G))
# Number of nodes: 16007
# Number of edges: 53917
# Average degree:   6.7367

# print(nx.degree(G)) # this is the list containing as many tuples as nodes, indicating each node's degree.

average_degree = statistics.mean([ tpl[1] for tpl in nx.degree(G) ])
# print(average_degree) # 6.736677703504717

proteins_degrees = nx.degree(G, nbunch=[n for n in G.nodes if descr_dict[n] == 'Human PPI (target)'])
# print("degree of protein nodes only", proteins_degrees)

average_degree_of_proteins = statistics.mean([ tpl[1] for tpl in proteins_degrees ])
print("average degree of proteins nodes", average_degree_of_proteins) # 

# Number of nodes: 16007
# Number of edges: 53917
# Average degree:   6.7367
# average degree of proteins nodes 115.66739606126914


# # # # #####################################################################################################
# # # # ################################# CALCULATION ON EDGES WITH FILTERING PROT-DRUG #####################
# # # # #####################################################################################################

# now we keep only the edges of proteins that interact with other proteins
prots_linked_to_drugs = []
for u,v,c in G.edges(data=True) :	
	if (descr_dict[u] == 'Human PPI (target)' and descr_dict[v] == 'Drug') :
		prots_linked_to_drugs.append(u)
print(len(prots_linked_to_drugs), len(set(prots_linked_to_drugs))) # 13018 169

# ({'Drug': 5703, 'Disease': 4176, 'GO': 3487, 'Symptom': 2157, 'Human PPI (target)': 457, 'Viral Gene': 27})

# print(sorted(prots_linked_to_drugs))

protein_counts = sorted(prots_linked_to_drugs)

protein_counts_dict = defaultdict( int )
for protein in protein_counts:
    protein_counts_dict[protein] += 1

my_proteins_nodes = [n for n in G.nodes if descr_dict[n] == 'Human PPI (target)']
for protein in my_proteins_nodes :
	if protein_counts_dict.get(protein) == None :
		# print("%s is missing in the linked proteins : Adding it with a frequency of 0.\n" % protein)
		protein_counts_dict[protein]= 0

average_degree_of_my_proteins_nodes = statistics.mean([ protein_counts_dict[protein] for protein in protein_counts_dict ])
print("average_degree_of_my_proteins_nodes", average_degree_of_my_proteins_nodes)
# 4.37417943107221

max_degree_of_my_proteins_nodes = max([ protein_counts_dict[protein] for protein in protein_counts_dict ])
print("max_degree_of_my_proteins_nodes", max_degree_of_my_proteins_nodes)
# average_degree_of_my_proteins_nodes 28.485776805251643
# max_degree_of_my_proteins_nodes 1105



ordered_protein_counts_dict = OrderedDict(sorted(protein_counts_dict.items()))


# with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\Further_Analysis\\Proteins\\COVID19_GDDS_proteins_degrees_to_drugs.tsv', 'a+') as file_write:
# 	for protein in ordered_protein_counts_dict :
# 		file_write.write("%s\t%s\t%s\n" % (protein, ordered_protein_counts_dict[protein], ordered_protein_counts_dict[protein] / sum(ordered_protein_counts_dict.values())))


# # # ##################################################################
# # # ############ SAME ANALYSIS FOR RANDOM NETWORKS ###################
# # # # ##################################################################

# with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\Further_Analysis\\Proteins\\Log_degree_analysis_on_random_networks_for_proteins_drugs.txt', 'a+') as log_write:

# 	nodefiles_numbers = []
# 	edgefiles_numbers = []
# 	for filename in os.listdir("A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\random_networks\\Mock_networks_21Apr"):
# 		if 'nodes_' in filename :
# 			filename = filename.split(".csv")[0]
# 			nodefiles_numbers.append(int(filename[19:]))
# 		if 'edges_' in filename :
# 			filename = filename.split(".csv")[0]
# 			edgefiles_numbers.append(int(filename[19:]))

# 	log_write.write("len(nodefiles_numbers) = %s\n" % len(nodefiles_numbers)) #2051
# 	log_write.write("len(edgefiles_numbers) = %s\n" % len(edgefiles_numbers)) #2051
# 	log_write.write("difference between edgefiles_numbers and nodefiles_numbers = %s\n" % len(set(nodefiles_numbers).difference(set(edgefiles_numbers)))) #0

# 	nodefiles_numbers = sorted(nodefiles_numbers)
# 	edgefiles_numbers = sorted(edgefiles_numbers)

# 	for iteration in nodefiles_numbers[1:] :
# 		print("\n########################## ITERATION %s ###########################\n" % iteration)
# 		log_write.write("\n########################## ITERATION %s ###########################\n" % iteration)
	
# 		with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\random_networks\\Mock_networks_21Apr\\COVID19_GDDS_nodes_%s.csv' % str(iteration), 'r') as nodecsv: # Open the file                   
# 			nodereader = csv.reader(nodecsv) # Read the csv  
# 			nodes = [n for n in nodereader][1:]                     
# 			node_names = [n[0] for n in nodes] # Get a list of only the node names 

# 		with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\random_networks\\Mock_networks_21Apr\\COVID19_GDDS_edges_%s.csv' % str(iteration), 'r') as edgecsv: # Open the file
# 			edgereader = csv.reader(edgecsv) # Read the csv     
# 			edges = [tuple(e) for e in edgereader][1:] # Retrieve the data

# 			log_write.write("len(node_names) : %s\n" % len(node_names)) #16007 --> 13448
# 			log_write.write("len(edges) : %s\n" % len(edges)) #53917 --> 37835

# 		G = nx.Graph()
# 		G.add_nodes_from(node_names)
# 		G.add_edges_from(edges)

# 		descr_dict = {}
# 		description_set=set()
# 		for node in nodes: 
# 			descr_dict[node[0]] = node[2]
# 			description_set.add(node[2])

# 		# print(description_set) #{'Symptom', 'Human PPI (target)', 'GO', 'Disease', 'Drug'}
# 		# recounted = Counter(list(descr_dict.values()))
# 		# print("%s\n" % recounted)

# 		# print(nx.info(G))
# 		log_write.write("\nRANDOM GRAPH %s\n" % iteration)
# 		log_write.write("%s\n" % nx.info(G))

# 		#####################################################################################################
# 		################################# GET EDGES CONNECTING proteinS TO FIRST ORDER PROTEINS ################
# 		#####################################################################################################

# 		prots_linked_to_drugs = []
# 		for u,v,c in G.edges(data=True) :	
# 			if (descr_dict[u] == 'Human PPI (target)' and descr_dict[v] == 'Drug') :
# 				prots_linked_to_drugs.append(u)
		
# 		protein_counts = sorted(prots_linked_to_drugs)

# 		protein_counts_dict = defaultdict( int )
# 		for protein in protein_counts:
# 		    protein_counts_dict[protein] += 1

# 		my_proteins_nodes = [n for n in G.nodes if descr_dict[n] == 'Human PPI (target)']
# 		for protein in my_proteins_nodes :
# 			if protein_counts_dict.get(protein) == None :
# 				# print("%s is missing in the linked proteins : Adding it with a frequency of 0.\n" % protein)
# 				protein_counts_dict[protein]= 0

# 		average_degree_of_my_proteins_nodes = statistics.mean([ protein_counts_dict[protein] for protein in protein_counts_dict ])
# 		log_write.write("Average degree of proteins linked to drug(s) : %s\n" % average_degree_of_my_proteins_nodes)

# 		ordered_protein_counts_dict = OrderedDict(sorted(protein_counts_dict.items()))

# 		with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\Further_Analysis\\Proteins\\prot_drug\\COVID19_GDDS_proteins_degrees_to_drugs_%s.tsv' % str(iteration), 'a+') as file_write:
# 			for protein in ordered_protein_counts_dict :
# 				if (sum(ordered_protein_counts_dict.values()) != 0) :
# 					file_write.write("%s\t%s\t%s\n" % (protein, ordered_protein_counts_dict[protein], ordered_protein_counts_dict[protein] / sum(ordered_protein_counts_dict.values())))
# 				else :
# 					file_write.write("%s\t%s\t%s\n" % (protein, ordered_protein_counts_dict[protein], ordered_protein_counts_dict[protein]))



# # # # ##################################################################
# # # # ############ COMPUTE DEGREES MEANS AND Z-SCORES ##################
# # # # ######## INCLUDING ABSENCE OF drugs AS ZERO DEGREES ###########
# # # # ##################################################################


# # For covid network, compute the mean and the standard deviation of proteins degrees and normalized degrees
covid_proteins = []
covid_protein_degree_dict = {}
covid_protein_norm_degree_dict = {}

total_protein_degree_dict = defaultdict(list)
total_protein_normalized_degree_dict = defaultdict(list)

with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\Further_Analysis\\Proteins\\COVID19_GDDS_proteins_degrees_to_drugs.tsv', 'r') as data :
	for line in data :
		line = line.split('\t')
		protein = line[0]
		covid_proteins.append(protein)
		protein_degree = int(line[1])
		protein_norm_degree = float(line[2].replace('\n', ''))
		if covid_protein_degree_dict.get(protein) == None :
			covid_protein_degree_dict[protein] = protein_degree
		else :
			print("warning : %s duplicate ?!" % protein)
		covid_protein_norm_degree_dict[protein] = protein_norm_degree

		if total_protein_degree_dict.get(protein) == None :
			total_protein_degree_dict[protein] = [protein_degree]
		else :
			total_protein_degree_dict[protein].append(protein_degree)

		if total_protein_normalized_degree_dict.get(protein) == None :
			total_protein_normalized_degree_dict[protein] = [protein_norm_degree]
		else :
			total_protein_normalized_degree_dict[protein].append(protein_norm_degree)

mean_protein_deg = statistics.mean(covid_protein_degree_dict.values())
mean_norm_protein_deg = statistics.mean(covid_protein_norm_degree_dict.values())
sd_protein_deg = statistics.pstdev(covid_protein_degree_dict.values())
sd_norm_protein_deg = statistics.pstdev(covid_protein_norm_degree_dict.values())
print(mean_protein_deg, sd_protein_deg, mean_norm_protein_deg, sd_norm_protein_deg) #
print(len(covid_proteins),len(covid_protein_degree_dict),len(covid_protein_norm_degree_dict)) #
# 28.485776805251643 119.87329340458344 0.002188183807439825 0.00920827265360143
# 457 457 457



nodefiles_numbers = []
edgefiles_numbers = []
for filename in os.listdir("A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\random_networks\\Mock_networks_21Apr"):
	if 'nodes_' in filename :
		filename = filename.split(".csv")[0]
		nodefiles_numbers.append(int(filename[19:]))
	if 'edges_' in filename :
		filename = filename.split(".csv")[0]
		edgefiles_numbers.append(int(filename[19:]))

# log_write.write("len(nodefiles_numbers) = %s\n" % len(nodefiles_numbers)) #2051
# log_write.write("len(edgefiles_numbers) = %s\n" % len(edgefiles_numbers)) #2051
# log_write.write("difference between edgefiles_numbers and nodefiles_numbers = %s\n" % len(set(nodefiles_numbers).difference(set(edgefiles_numbers)))) #0
print("len(nodefiles_numbers) = %s\n" % len(nodefiles_numbers)) #2051
print("len(edgefiles_numbers) = %s\n" % len(edgefiles_numbers)) #2051
print("difference between edgefiles_numbers and nodefiles_numbers = %s\n" % len(set(nodefiles_numbers).difference(set(edgefiles_numbers)))) #0

# nodefiles_numbers = sorted(nodefiles_numbers)
# edgefiles_numbers = sorted(edgefiles_numbers)


for iteration in nodefiles_numbers[1:] :
	# log_write.write("\n########################## ITERATION %s ###########################\n" % iteration)
	print("\n########################## ITERATION %s ###########################\n" % iteration)

	with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\Further_Analysis\\Proteins\\prot_drug\\COVID19_GDDS_proteins_degrees_to_drugs_%s.tsv' % str(iteration), 'r') as data: # Open the file  
		for line in data :
			line = line.split('\t')
			protein = line[0]
			degree = float(line[1])
			normalized_degree = float(line[2].replace('\n', ''))

			if total_protein_degree_dict.get(protein) == None :
				total_protein_degree_dict[protein] = [degree]
			else :
				total_protein_degree_dict[protein].append(degree)

			if total_protein_normalized_degree_dict.get(protein) == None :
				total_protein_normalized_degree_dict[protein] = [normalized_degree]
			else :
				total_protein_normalized_degree_dict[protein].append(normalized_degree)

print(len(total_protein_degree_dict)) # 19928

print(len(total_protein_normalized_degree_dict)) # 19928

### make the files with all calculations
with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\Further_Analysis\\Proteins\\COVID19_GDDS_Results_Proteins_to_drugs-withZeros.tsv', 'a+') as file_write:
	file_write.write("protein\t# of nets it is found in\tmin-degree\tmax-degree\tmean-degree\tmedian-degree\tsd-degree\n")
	for protein in total_protein_degree_dict :
		nets_number = len(total_protein_degree_dict[protein])
		min_deg = min(total_protein_degree_dict[protein])
		max_deg = max(total_protein_degree_dict[protein])
		mean_deg = statistics.mean(total_protein_degree_dict[protein])
		median_deg = statistics.median(total_protein_degree_dict[protein])
		sd_deg = statistics.pstdev(total_protein_degree_dict[protein])
		file_write.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (protein, nets_number, min_deg, max_deg, mean_deg, median_deg, sd_deg))


with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\Further_Analysis\\Proteins\\COVID19_GDDS_Results_Proteins_to_drugs_Relative-withZeros.tsv', 'a+') as file_write:
	file_write.write("protein\t# of nets it is found in\tmin-degree-rel\tmax-degree-rel\tmean-degree-rel\tmedian-degree-rel\tsd-degree-rel\n")
	for protein in total_protein_normalized_degree_dict :
		nets_number = len(total_protein_normalized_degree_dict[protein])
		min_deg = min(total_protein_normalized_degree_dict[protein])
		max_deg = max(total_protein_normalized_degree_dict[protein])
		mean_deg = statistics.mean(total_protein_normalized_degree_dict[protein])
		median_deg = statistics.median(total_protein_normalized_degree_dict[protein])
		sd_deg = statistics.pstdev(total_protein_normalized_degree_dict[protein])
		file_write.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (protein, nets_number, min_deg, max_deg, mean_deg, median_deg, sd_deg))

### compute Z-score for the covid proteins
# For each protein found in covid network, compute 2 Z-scores.

with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\Further_Analysis\\Proteins\\COVID19_GDDS_Results_ForCovidProteins_to_drugs-withZeros.tsv', 'a+') as file_write:
	file_write.write("covid protein\t# of nets it is found in\tmin-degree\tmax-degree\tmean-degree\tmedian-degree\tsd-degree\tcovid-degree\tZ-score\n")
	for protein in covid_proteins :
		nets_number = len(total_protein_degree_dict[protein])
		min_deg = min(total_protein_degree_dict[protein])
		max_deg = max(total_protein_degree_dict[protein])
		mean_deg = statistics.mean(total_protein_degree_dict[protein])
		median_deg = statistics.median(total_protein_degree_dict[protein])
		sd_deg = statistics.pstdev(total_protein_degree_dict[protein])
		covid_deg = covid_protein_degree_dict[protein]
		if sd_deg != 0 :
			z_score = ( covid_deg - mean_deg) / sd_deg
		else : 
			z_score = 'NaN'
		file_write.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (protein, nets_number, min_deg, max_deg, mean_deg, median_deg, sd_deg, covid_deg, z_score))


with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\Further_Analysis\\Proteins\\COVID19_GDDS_Results_ForCovidProteins_to_drugs_Relative-withZeros.tsv', 'a+') as file_write:
	file_write.write("covid protein\t# of nets it is found in\tmin-degree-rel\tmax-degree-rel\tmean-degree-rel\tmedian-degree-rel\tsd-degree-rel\tcovid-degree-rel\tZ-score-rel\n")
	for protein in covid_proteins :
		nets_number = len(total_protein_normalized_degree_dict[protein])
		min_deg = min(total_protein_normalized_degree_dict[protein])
		max_deg = max(total_protein_normalized_degree_dict[protein])
		mean_deg = statistics.mean(total_protein_normalized_degree_dict[protein])
		median_deg = statistics.median(total_protein_normalized_degree_dict[protein])
		sd_deg = statistics.pstdev(total_protein_normalized_degree_dict[protein])
		covid_deg = covid_protein_norm_degree_dict[protein]
		if sd_deg != 0 :
			z_score = ( covid_deg - mean_deg) / sd_deg
		else : 
			z_score = 'NaN'
		file_write.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (protein, nets_number, min_deg, max_deg, mean_deg, median_deg, sd_deg, covid_deg, z_score))










stop = timeit.default_timer()
print(stop - start)
 
protein-drug-degrees-undirected-rankings.py,
#!/usr/bin/python
# -*- coding: utf-8 -*-
# jupyter notebook --browser='C:/Program Files (x86)/Google/Chrome/Application/chrome.exe %s'


import timeit
start = timeit.default_timer()

import os
import csv
from operator import itemgetter
import networkx as nx
from networkx.algorithms import community
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
from collections import Counter
from collections import defaultdict
from collections import OrderedDict
import statistics 
import re
import gzip





#####################################################################################################
#########################  Compute network Protein degree distribution #################################
##########################################  based on edges to proteins #####################################
#####################################################################################################

with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\COVID19_GDDS_nodes.csv', 'r') as nodecsv: # Open the file                       
	nodereader = csv.reader(nodecsv) # Read the csv  
	nodes = [n for n in nodereader][1:]                     
	node_names = [n[0] for n in nodes] # Get a list of only the node names 

with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\COVID19_GDDS_edges.csv', 'r') as edgecsv: # Open the file
	edgereader = csv.reader(edgecsv) # Read the csv     
	edges = [tuple(e) for e in edgereader][1:] # Retrieve the data

G = nx.Graph()
G.add_nodes_from(node_names)
G.add_edges_from(edges)

descr_dict = {}
for node in nodes: 
	descr_dict[node[0]] = node[2]

print(nx.info(G))
# Number of nodes: 16007
# Number of edges: 53917
# Average degree:   6.7367

# print(nx.degree(G)) # this is the list containing as many tuples as nodes, indicating each node's degree.

average_degree = statistics.mean([ tpl[1] for tpl in nx.degree(G) ])
# print(average_degree) # 6.736677703504717

proteins_degrees = nx.degree(G, nbunch=[n for n in G.nodes if descr_dict[n] == 'Human PPI (target)'])
# print("degree of protein nodes only", proteins_degrees)

average_degree_of_proteins = statistics.mean([ tpl[1] for tpl in proteins_degrees ])
print("average degree of proteins nodes", average_degree_of_proteins) # 

# Number of nodes: 16007
# Number of edges: 53917
# Average degree:   6.7367
# average degree of proteins nodes 115.66739606126914


# # # # #####################################################################################################
# # # # ################################# CALCULATION ON EDGES WITH FILTERING PROT-DRUG #####################
# # # # #####################################################################################################

# now we keep only the edges of proteins that interact with other proteins
prots_linked_to_drugs = []
for u,v,c in G.edges(data=True) :	
	if (descr_dict[u] == 'Human PPI (target)' and descr_dict[v] == 'Drug') :
		prots_linked_to_drugs.append(u)
	if (descr_dict[v] == 'Human PPI (target)' and descr_dict[u] == 'Drug') :
		prots_linked_to_drugs.append(v)		
print(len(prots_linked_to_drugs), len(set(prots_linked_to_drugs))) # 13018 169

# ({'Drug': 5703, 'Disease': 4176, 'GO': 3487, 'Symptom': 2157, 'Human PPI (target)': 457, 'Viral Gene': 27})

# print(sorted(prots_linked_to_drugs))

protein_counts = sorted(prots_linked_to_drugs)

protein_counts_dict = defaultdict( int )
for protein in protein_counts:
    protein_counts_dict[protein] += 1

my_proteins_nodes = [n for n in G.nodes if descr_dict[n] == 'Human PPI (target)']
for protein in my_proteins_nodes :
	if protein_counts_dict.get(protein) == None :
		# print("%s is missing in the linked proteins : Adding it with a frequency of 0.\n" % protein)
		protein_counts_dict[protein]= 0

average_degree_of_my_proteins_nodes = statistics.mean([ protein_counts_dict[protein] for protein in protein_counts_dict ])
print("average_degree_of_my_proteins_nodes", average_degree_of_my_proteins_nodes)
# 4.37417943107221

max_degree_of_my_proteins_nodes = max([ protein_counts_dict[protein] for protein in protein_counts_dict ])
print("max_degree_of_my_proteins_nodes", max_degree_of_my_proteins_nodes)
# average_degree_of_my_proteins_nodes 28.485776805251643
# max_degree_of_my_proteins_nodes 1105



ordered_protein_counts_dict = OrderedDict(sorted(protein_counts_dict.items()))


# with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\Further_Analysis\\Proteins\\COVID19_GDDS_proteins_degrees_to_drugs_undirected.tsv', 'a+') as file_write:
# 	for protein in ordered_protein_counts_dict :
# 		file_write.write("%s\t%s\t%s\n" % (protein, ordered_protein_counts_dict[protein], ordered_protein_counts_dict[protein] / sum(ordered_protein_counts_dict.values())))


# # # ##################################################################
# # # ############ SAME ANALYSIS FOR RANDOM NETWORKS ###################
# # # # ##################################################################

# with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\Further_Analysis\\Proteins\\Log_degree_analysis_on_random_networks_for_proteins_drugs_undirected_corrected.txt', 'a+') as log_write:

# 	nodefiles_numbers = []
# 	edgefiles_numbers = []
# 	for filename in os.listdir("A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\random_networks\\Mock_networks_21Apr"):
# 		if 'nodes_' in filename :
# 			filename = filename.split(".csv")[0]
# 			nodefiles_numbers.append(int(filename[19:]))
# 		if 'edges_' in filename :
# 			filename = filename.split(".csv")[0]
# 			edgefiles_numbers.append(int(filename[19:]))

# 	log_write.write("len(nodefiles_numbers) = %s\n" % len(nodefiles_numbers)) #2051
# 	log_write.write("len(edgefiles_numbers) = %s\n" % len(edgefiles_numbers)) #2051
# 	log_write.write("difference between edgefiles_numbers and nodefiles_numbers = %s\n" % len(set(nodefiles_numbers).difference(set(edgefiles_numbers)))) #0

# 	nodefiles_numbers = sorted(nodefiles_numbers)
# 	edgefiles_numbers = sorted(edgefiles_numbers)

# 	networks_to_avoid = [3, 11, 22, 23, 33, 34, 38, 57, 72, 77, 123, 127, 148, 156, 157, 183, 205, 209, 210, 217, 232, 252, 268, 281, 291, 312, 341, 353, 362, 364, 413, 421, 433, 443, 451, 494, 525, 532, 569, 592, 629, 658, 662, 674, 676, 701, 705, 712, 736, 757, 788, 796, 825, 827, 838, 873, 876, 908, 910, 916, 917, 918, 938, 958, 965, 1088, 1125, 1128, 1145, 1146, 1165, 1171, 1183, 1185, 1192, 1195, 1211, 1244, 1254, 1257, 1300, 1353, 1354, 1364, 1372, 1394, 1395, 1401, 1410, 1427, 1454, 1477, 1497, 1517, 1532, 1549, 1551, 1575, 1578, 1611, 1624, 1629, 1637, 1650, 1655, 1697, 1700, 1713, 1736, 1745, 1762, 1790, 1800, 1815, 1824, 1827, 1866, 1872, 1886, 1887, 2053, 2072, 2106, 2115, 2128, 2160, 2217, 2343, 2383, 2388, 2413, 2458, 2465]
# 	nodefiles_numbers = [item for item in nodefiles_numbers if item not in networks_to_avoid]

# 	for iteration in nodefiles_numbers[1:] :
# 		print("\n########################## ITERATION %s ###########################\n" % iteration)
# 		log_write.write("\n########################## ITERATION %s ###########################\n" % iteration)
	
# 		with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\random_networks\\Mock_networks_21Apr\\COVID19_GDDS_nodes_%s.csv' % str(iteration), 'r') as nodecsv: # Open the file                   
# 			nodereader = csv.reader(nodecsv) # Read the csv  
# 			nodes = [n for n in nodereader][1:]                     
# 			node_names = [n[0] for n in nodes] # Get a list of only the node names 

# 		with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\random_networks\\Mock_networks_21Apr\\COVID19_GDDS_edges_%s.csv' % str(iteration), 'r') as edgecsv: # Open the file
# 			edgereader = csv.reader(edgecsv) # Read the csv     
# 			edges = [tuple(e) for e in edgereader][1:] # Retrieve the data

# 			log_write.write("len(node_names) : %s\n" % len(node_names)) #16007 --> 13448
# 			log_write.write("len(edges) : %s\n" % len(edges)) #53917 --> 37835

# 		G = nx.Graph()
# 		G.add_nodes_from(node_names)
# 		G.add_edges_from(edges)

# 		descr_dict = {}
# 		description_set=set()
# 		for node in nodes: 
# 			descr_dict[node[0]] = node[2]
# 			description_set.add(node[2])

# 		# print(description_set) #{'Symptom', 'Human PPI (target)', 'GO', 'Disease', 'Drug'}
# 		# recounted = Counter(list(descr_dict.values()))
# 		# print("%s\n" % recounted)

# 		# print(nx.info(G))
# 		log_write.write("\nRANDOM GRAPH %s\n" % iteration)
# 		log_write.write("%s\n" % nx.info(G))

# 		#####################################################################################################
# 		################################# GET EDGES CONNECTING proteinS TO FIRST ORDER PROTEINS ################
# 		#####################################################################################################

# 		prots_linked_to_drugs = []
# 		for u,v,c in G.edges(data=True) :	
# 			if (descr_dict[u] == 'Human PPI (target)' and descr_dict[v] == 'Drug') :
# 				prots_linked_to_drugs.append(u)
# 			if (descr_dict[v] == 'Human PPI (target)' and descr_dict[u] == 'Drug') :
# 				prots_linked_to_drugs.append(v)		
		
# 		protein_counts = sorted(prots_linked_to_drugs)

# 		protein_counts_dict = defaultdict( int )
# 		for protein in protein_counts:
# 		    protein_counts_dict[protein] += 1

# 		my_proteins_nodes = [n for n in G.nodes if descr_dict[n] == 'Human PPI (target)']
# 		for protein in my_proteins_nodes :
# 			if protein_counts_dict.get(protein) == None :
# 				# print("%s is missing in the linked proteins : Adding it with a frequency of 0.\n" % protein)
# 				protein_counts_dict[protein]= 0

# 		average_degree_of_my_proteins_nodes = statistics.mean([ protein_counts_dict[protein] for protein in protein_counts_dict ])
# 		log_write.write("Average degree of proteins linked to drug(s) : %s\n" % average_degree_of_my_proteins_nodes)

# 		ordered_protein_counts_dict = OrderedDict(sorted(protein_counts_dict.items()))

# 		with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\Further_Analysis\\Proteins\\prot_drug_undirected_corrected\\COVID19_GDDS_proteins_degrees_to_drugs_%s.tsv' % str(iteration), 'a+') as file_write:
# 			for protein in ordered_protein_counts_dict :
# 				if (sum(ordered_protein_counts_dict.values()) != 0) :
# 					file_write.write("%s\t%s\t%s\n" % (protein, ordered_protein_counts_dict[protein], ordered_protein_counts_dict[protein] / sum(ordered_protein_counts_dict.values())))
# 				else :
# 					file_write.write("%s\t%s\t%s\n" % (protein, ordered_protein_counts_dict[protein], ordered_protein_counts_dict[protein]))



# # # # ##################################################################
# # # # ############ COMPUTE DEGREES MEANS AND Z-SCORES ##################
# # # # ######## INCLUDING ABSENCE OF drugs AS ZERO DEGREES ###########
# # # # ##################################################################


# # For covid network, compute the mean and the standard deviation of proteins degrees and normalized degrees
covid_proteins = []
covid_protein_degree_dict = {}
covid_protein_norm_degree_dict = {}

total_protein_degree_dict = defaultdict(list)
total_protein_normalized_degree_dict = defaultdict(list)

with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\Further_Analysis\\Proteins\\COVID19_GDDS_proteins_degrees_to_drugs_undirected.tsv', 'r') as data :
	for line in data :
		line = line.split('\t')
		protein = line[0]
		covid_proteins.append(protein)
		protein_degree = int(line[1])
		protein_norm_degree = float(line[2].replace('\n', ''))
		if covid_protein_degree_dict.get(protein) == None :
			covid_protein_degree_dict[protein] = protein_degree
		else :
			print("warning : %s duplicate ?!" % protein)
		covid_protein_norm_degree_dict[protein] = protein_norm_degree

		if total_protein_degree_dict.get(protein) == None :
			total_protein_degree_dict[protein] = [protein_degree]
		else :
			total_protein_degree_dict[protein].append(protein_degree)

		if total_protein_normalized_degree_dict.get(protein) == None :
			total_protein_normalized_degree_dict[protein] = [protein_norm_degree]
		else :
			total_protein_normalized_degree_dict[protein].append(protein_norm_degree)

mean_protein_deg = statistics.mean(covid_protein_degree_dict.values())
mean_norm_protein_deg = statistics.mean(covid_protein_norm_degree_dict.values())
sd_protein_deg = statistics.pstdev(covid_protein_degree_dict.values())
sd_norm_protein_deg = statistics.pstdev(covid_protein_norm_degree_dict.values())
print(mean_protein_deg, sd_protein_deg, mean_norm_protein_deg, sd_norm_protein_deg) #
print(len(covid_proteins),len(covid_protein_degree_dict),len(covid_protein_norm_degree_dict)) #
# 28.485776805251643 119.87329340458344 0.002188183807439825 0.00920827265360143
# 457 457 457




nodefiles_numbers = []
edgefiles_numbers = []
for filename in os.listdir("A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\random_networks\\Mock_networks_21Apr"):
	if 'nodes_' in filename :
		filename = filename.split(".csv")[0]
		nodefiles_numbers.append(int(filename[19:]))
	if 'edges_' in filename :
		filename = filename.split(".csv")[0]
		edgefiles_numbers.append(int(filename[19:]))

# log_write.write("len(nodefiles_numbers) = %s\n" % len(nodefiles_numbers)) #2051
# log_write.write("len(edgefiles_numbers) = %s\n" % len(edgefiles_numbers)) #2051
# log_write.write("difference between edgefiles_numbers and nodefiles_numbers = %s\n" % len(set(nodefiles_numbers).difference(set(edgefiles_numbers)))) #0
print("len(nodefiles_numbers) = %s\n" % len(nodefiles_numbers)) #2051
print("len(edgefiles_numbers) = %s\n" % len(edgefiles_numbers)) #2051
print("difference between edgefiles_numbers and nodefiles_numbers = %s\n" % len(set(nodefiles_numbers).difference(set(edgefiles_numbers)))) #0

# nodefiles_numbers = sorted(nodefiles_numbers)
# edgefiles_numbers = sorted(edgefiles_numbers)
networks_to_avoid = [3, 11, 22, 23, 33, 34, 38, 57, 72, 77, 123, 127, 148, 156, 157, 183, 205, 209, 210, 217, 232, 252, 268, 281, 291, 312, 341, 353, 362, 364, 413, 421, 433, 443, 451, 494, 525, 532, 569, 592, 629, 658, 662, 674, 676, 701, 705, 712, 736, 757, 788, 796, 825, 827, 838, 873, 876, 908, 910, 916, 917, 918, 938, 958, 965, 1088, 1125, 1128, 1145, 1146, 1165, 1171, 1183, 1185, 1192, 1195, 1211, 1244, 1254, 1257, 1300, 1353, 1354, 1364, 1372, 1394, 1395, 1401, 1410, 1427, 1454, 1477, 1497, 1517, 1532, 1549, 1551, 1575, 1578, 1611, 1624, 1629, 1637, 1650, 1655, 1697, 1700, 1713, 1736, 1745, 1762, 1790, 1800, 1815, 1824, 1827, 1866, 1872, 1886, 1887, 2053, 2072, 2106, 2115, 2128, 2160, 2217, 2343, 2383, 2388, 2413, 2458, 2465]
nodefiles_numbers = [item for item in nodefiles_numbers if item not in networks_to_avoid]
print("len(nodefiles_numbers) = %s\n" % len(nodefiles_numbers)) #1918



for iteration in nodefiles_numbers[1:] :
	# log_write.write("\n########################## ITERATION %s ###########################\n" % iteration)
	# print("\n########################## ITERATION %s ###########################\n" % iteration)

	with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\Further_Analysis\\Proteins\\prot_drug_undirected_corrected\\COVID19_GDDS_proteins_degrees_to_drugs_%s.tsv' % str(iteration), 'r') as data: # Open the file  
		for line in data :
			line = line.split('\t')
			protein = line[0]
			degree = float(line[1])
			normalized_degree = float(line[2].replace('\n', ''))

			if total_protein_degree_dict.get(protein) == None :
				total_protein_degree_dict[protein] = [degree]
			else :
				total_protein_degree_dict[protein].append(degree)

			if total_protein_normalized_degree_dict.get(protein) == None :
				total_protein_normalized_degree_dict[protein] = [normalized_degree]
			else :
				total_protein_normalized_degree_dict[protein].append(normalized_degree)

print(len(total_protein_degree_dict)) # 19928

print(len(total_protein_normalized_degree_dict)) # 19928

### make the files with all calculations
with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\Further_Analysis\\Proteins\\COVID19_GDDS_Results_Proteins_to_drugs-withZeros_undirected_corrected.tsv', 'a+') as file_write:
	file_write.write("protein\t# of nets it is found in\tmin-degree\tmax-degree\tmean-degree\tmedian-degree\tsd-degree\n")
	for protein in total_protein_degree_dict :
		nets_number = len(total_protein_degree_dict[protein])
		min_deg = min(total_protein_degree_dict[protein])
		max_deg = max(total_protein_degree_dict[protein])
		mean_deg = statistics.mean(total_protein_degree_dict[protein])
		median_deg = statistics.median(total_protein_degree_dict[protein])
		sd_deg = statistics.pstdev(total_protein_degree_dict[protein])
		file_write.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (protein, nets_number, min_deg, max_deg, mean_deg, median_deg, sd_deg))


with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\Further_Analysis\\Proteins\\COVID19_GDDS_Results_Proteins_to_drugs_Relative-withZeros_undirected_corrected.tsv', 'a+') as file_write:
	file_write.write("protein\t# of nets it is found in\tmin-degree-rel\tmax-degree-rel\tmean-degree-rel\tmedian-degree-rel\tsd-degree-rel\n")
	for protein in total_protein_normalized_degree_dict :
		nets_number = len(total_protein_normalized_degree_dict[protein])
		min_deg = min(total_protein_normalized_degree_dict[protein])
		max_deg = max(total_protein_normalized_degree_dict[protein])
		mean_deg = statistics.mean(total_protein_normalized_degree_dict[protein])
		median_deg = statistics.median(total_protein_normalized_degree_dict[protein])
		sd_deg = statistics.pstdev(total_protein_normalized_degree_dict[protein])
		file_write.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (protein, nets_number, min_deg, max_deg, mean_deg, median_deg, sd_deg))

### compute Z-score for the covid proteins
# For each protein found in covid network, compute 2 Z-scores.

with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\Further_Analysis\\Proteins\\COVID19_GDDS_Results_ForCovidProteins_to_drugs-withZeros_undirected_corrected.tsv', 'a+') as file_write:
	file_write.write("covid protein\t# of nets it is found in\tmin-degree\tmax-degree\tmean-degree\tmedian-degree\tsd-degree\tcovid-degree\tZ-score\n")
	for protein in covid_proteins :
		nets_number = len(total_protein_degree_dict[protein])
		min_deg = min(total_protein_degree_dict[protein])
		max_deg = max(total_protein_degree_dict[protein])
		mean_deg = statistics.mean(total_protein_degree_dict[protein])
		median_deg = statistics.median(total_protein_degree_dict[protein])
		sd_deg = statistics.pstdev(total_protein_degree_dict[protein])
		covid_deg = covid_protein_degree_dict[protein]
		if sd_deg != 0 :
			z_score = ( covid_deg - mean_deg) / sd_deg
		else : 
			z_score = 'NaN'
		file_write.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (protein, nets_number, min_deg, max_deg, mean_deg, median_deg, sd_deg, covid_deg, z_score))


with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\Further_Analysis\\Proteins\\COVID19_GDDS_Results_ForCovidProteins_to_drugs_Relative-withZeros_undirected_corrected.tsv', 'a+') as file_write:
	file_write.write("covid protein\t# of nets it is found in\tmin-degree-rel\tmax-degree-rel\tmean-degree-rel\tmedian-degree-rel\tsd-degree-rel\tcovid-degree-rel\tZ-score-rel\n")
	for protein in covid_proteins :
		nets_number = len(total_protein_normalized_degree_dict[protein])
		min_deg = min(total_protein_normalized_degree_dict[protein])
		max_deg = max(total_protein_normalized_degree_dict[protein])
		mean_deg = statistics.mean(total_protein_normalized_degree_dict[protein])
		median_deg = statistics.median(total_protein_normalized_degree_dict[protein])
		sd_deg = statistics.pstdev(total_protein_normalized_degree_dict[protein])
		covid_deg = covid_protein_norm_degree_dict[protein]
		if sd_deg != 0 :
			z_score = ( covid_deg - mean_deg) / sd_deg
		else : 
			z_score = 'NaN'
		file_write.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (protein, nets_number, min_deg, max_deg, mean_deg, median_deg, sd_deg, covid_deg, z_score))










stop = timeit.default_timer()
print(stop - start)
 
protein-protein-degrees-rankings.py,
#!/usr/bin/python
# -*- coding: utf-8 -*-
# jupyter notebook --browser='C:/Program Files (x86)/Google/Chrome/Application/chrome.exe %s'


import timeit
start = timeit.default_timer()

import os
import csv
from operator import itemgetter
import networkx as nx
from networkx.algorithms import community
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
from collections import Counter
from collections import defaultdict
from collections import OrderedDict
import statistics 
import re
import gzip





#####################################################################################################
#########################  Compute network Protein degree distribution #################################
##########################################  based on edges to proteins #####################################
#####################################################################################################

with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\COVID19_GDDS_nodes.csv', 'r') as nodecsv: # Open the file                       
	nodereader = csv.reader(nodecsv) # Read the csv  
	nodes = [n for n in nodereader][1:]                     
	node_names = [n[0] for n in nodes] # Get a list of only the node names 

with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\COVID19_GDDS_edges.csv', 'r') as edgecsv: # Open the file
	edgereader = csv.reader(edgecsv) # Read the csv     
	edges = [tuple(e) for e in edgereader][1:] # Retrieve the data

G = nx.Graph()
G.add_nodes_from(node_names)
G.add_edges_from(edges)

descr_dict = {}
for node in nodes: 
	descr_dict[node[0]] = node[2]

print(nx.info(G))
# Number of nodes: 16007
# Number of edges: 53917
# Average degree:   6.7367

# print(nx.degree(G)) # this is the list containing as many tuples as nodes, indicating each node's degree.

average_degree = statistics.mean([ tpl[1] for tpl in nx.degree(G) ])
# print(average_degree) # 6.736677703504717

proteins_degrees = nx.degree(G, nbunch=[n for n in G.nodes if descr_dict[n] == 'Human PPI (target)'])
# print("degree of protein nodes only", proteins_degrees)

average_degree_of_proteins = statistics.mean([ tpl[1] for tpl in proteins_degrees ])
print("average degree of proteins nodes", average_degree_of_proteins) # 

# Number of nodes: 16007
# Number of edges: 53917
# Average degree:   6.7367
# average degree of proteins nodes 115.66739606126914



######### temp check if edges are symmetric ##########

# set_of_edges = []
# for u,v,c in G.edges(data=True) :	
# 	if (descr_dict[u] == 'Human PPI (target)' and descr_dict[v] == 'Human PPI (target)') :
# 		set_of_edges.append((u,v))
# print(len(set_of_edges), len(set(set_of_edges))) #1999 1999
# for (a,b) in set_of_edges :
# 	if (b,a) in set_of_edges :
# 		print("(%s, %s) and (%s, %s) are both found in edges" % (a, b, b, a)) # none
#########################
# set_of_edges = [('coucou', 'hello'),('hello', 'coucou'),(1,2),(2,10)]
# for (a,b) in set_of_edges :
# 	if (b,a) in set_of_edges :
# 		print("(%s, %s) and (%s, %s) are both found in edges" % (a, b, b, a)) # none



# ######### temp check if CEP68 edges ##########

# set_of_edges = []
# for u,v,c in G.edges(data=True) :	
# 	if (u == 'MYC' and descr_dict[v] == 'Human PPI (target)') :
# 		set_of_edges.append((u,v))
# 	if (v == 'MYC' and descr_dict[u] == 'Human PPI (target)') :
# 		set_of_edges.append((u,v))
# print(len(set_of_edges), len(set(set_of_edges))) #1999 1999


# for i in set_of_edges :
# 	print(i)

#########################







# # # # #####################################################################################################
# # # # ################################# CALCULATION ON EDGES WITH FILTERING PROT-PROT #####################
# # # # #####################################################################################################

# # now we keep only the edges of proteins that interact with other proteins
# prots_linked_to_other_prots = []
# for u,v,c in G.edges(data=True) :	
# 	if (descr_dict[u] == 'Human PPI (target)' and descr_dict[v] == 'Human PPI (target)') :
# 		prots_linked_to_other_prots.append(u)
# 		prots_linked_to_other_prots.append(v)
# print(len(prots_linked_to_other_prots), len(set(prots_linked_to_other_prots))) # 3998 425

# # ({'Drug': 5703, 'Disease': 4176, 'GO': 3487, 'Symptom': 2157, 'Human PPI (target)': 457, 'Viral Gene': 27})

# # print(sorted(prots_linked_to_other_prots))

# protein_counts = sorted(prots_linked_to_other_prots)

# protein_counts_dict = defaultdict( int )
# for protein in protein_counts:
#     protein_counts_dict[protein] += 1

# my_proteins_nodes = [n for n in G.nodes if descr_dict[n] == 'Human PPI (target)']
# for protein in my_proteins_nodes :
# 	if protein_counts_dict.get(protein) == None :
# 		# print("%s is missing in the linked proteins : Adding it with a frequency of 0.\n" % protein)
# 		protein_counts_dict[protein]= 0

# average_degree_of_my_proteins_nodes = statistics.mean([ protein_counts_dict[protein] for protein in protein_counts_dict ])
# print("average_degree_of_my_proteins_nodes", average_degree_of_my_proteins_nodes)
# # 4.37417943107221

# max_degree_of_my_proteins_nodes = max([ protein_counts_dict[protein] for protein in protein_counts_dict ])
# print("max_degree_of_my_proteins_nodes", max_degree_of_my_proteins_nodes) #55
# # average_degree_of_my_proteins_nodes 8.74835886214442
# # max_degree_of_my_proteins_nodes 108


# ordered_protein_counts_dict = OrderedDict(sorted(protein_counts_dict.items()))


# with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\Further_Analysis\\Proteins\\COVID19_GDDS_proteins_degrees_to_other_proteins.tsv', 'a+') as file_write:
# 	for protein in ordered_protein_counts_dict :
# 		file_write.write("%s\t%s\t%s\n" % (protein, ordered_protein_counts_dict[protein], ordered_protein_counts_dict[protein] / sum(ordered_protein_counts_dict.values())))


# # # ##################################################################
# # # ############ SAME ANALYSIS FOR RANDOM NETWORKS ###################
# # # # ##################################################################

# with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\Further_Analysis\\Proteins\\Log_degree_analysis_on_random_networks_for_proteins_proteins.txt', 'a+') as log_write:

# 	nodefiles_numbers = []
# 	edgefiles_numbers = []
# 	for filename in os.listdir("A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\random_networks\\Mock_networks_21Apr"):
# 		if 'nodes_' in filename :
# 			filename = filename.split(".csv")[0]
# 			nodefiles_numbers.append(int(filename[19:]))
# 		if 'edges_' in filename :
# 			filename = filename.split(".csv")[0]
# 			edgefiles_numbers.append(int(filename[19:]))

# 	log_write.write("len(nodefiles_numbers) = %s\n" % len(nodefiles_numbers)) #2051
# 	log_write.write("len(edgefiles_numbers) = %s\n" % len(edgefiles_numbers)) #2051
# 	log_write.write("difference between edgefiles_numbers and nodefiles_numbers = %s\n" % len(set(nodefiles_numbers).difference(set(edgefiles_numbers)))) #0

# 	nodefiles_numbers = sorted(nodefiles_numbers)
# 	edgefiles_numbers = sorted(edgefiles_numbers)

# 	for iteration in nodefiles_numbers[1:] :
# 		print("\n########################## ITERATION %s ###########################\n" % iteration)
# 		log_write.write("\n########################## ITERATION %s ###########################\n" % iteration)
	
# 		with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\random_networks\\Mock_networks_21Apr\\COVID19_GDDS_nodes_%s.csv' % str(iteration), 'r') as nodecsv: # Open the file                   
# 			nodereader = csv.reader(nodecsv) # Read the csv  
# 			nodes = [n for n in nodereader][1:]                     
# 			node_names = [n[0] for n in nodes] # Get a list of only the node names 

# 		with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\random_networks\\Mock_networks_21Apr\\COVID19_GDDS_edges_%s.csv' % str(iteration), 'r') as edgecsv: # Open the file
# 			edgereader = csv.reader(edgecsv) # Read the csv     
# 			edges = [tuple(e) for e in edgereader][1:] # Retrieve the data

# 			log_write.write("len(node_names) : %s\n" % len(node_names)) #16007 --> 13448
# 			log_write.write("len(edges) : %s\n" % len(edges)) #53917 --> 37835

# 		G = nx.Graph()
# 		G.add_nodes_from(node_names)
# 		G.add_edges_from(edges)

# 		descr_dict = {}
# 		description_set=set()
# 		for node in nodes: 
# 			descr_dict[node[0]] = node[2]
# 			description_set.add(node[2])

# 		# print(description_set) #{'Symptom', 'Human PPI (target)', 'GO', 'Disease', 'Drug'}
# 		# recounted = Counter(list(descr_dict.values()))
# 		# print("%s\n" % recounted)

# 		# print(nx.info(G))
# 		log_write.write("\nRANDOM GRAPH %s\n" % iteration)
# 		log_write.write("%s\n" % nx.info(G))

# 		#####################################################################################################
# 		################################# GET EDGES CONNECTING proteinS TO FIRST ORDER PROTEINS ################
# 		#####################################################################################################

# 		prots_linked_to_other_prots = []
# 		for u,v,c in G.edges(data=True) :	
# 			if (descr_dict[u] == 'Human PPI (target)' and descr_dict[v] == 'Human PPI (target)') :
# 				prots_linked_to_other_prots.append(u)
# 				prots_linked_to_other_prots.append(v)
		
# 		protein_counts = sorted(prots_linked_to_other_prots)

# 		protein_counts_dict = defaultdict( int )
# 		for protein in protein_counts:
# 		    protein_counts_dict[protein] += 1

# 		my_proteins_nodes = [n for n in G.nodes if descr_dict[n] == 'Human PPI (target)']
# 		for protein in my_proteins_nodes :
# 			if protein_counts_dict.get(protein) == None :
# 				# print("%s is missing in the linked proteins : Adding it with a frequency of 0.\n" % protein)
# 				protein_counts_dict[protein]= 0

# 		average_degree_of_my_proteins_nodes = statistics.mean([ protein_counts_dict[protein] for protein in protein_counts_dict ])
# 		log_write.write("Average degree of proteins linked to another protein : %s\n" % average_degree_of_my_proteins_nodes)

# 		ordered_protein_counts_dict = OrderedDict(sorted(protein_counts_dict.items()))

# 		with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\Further_Analysis\\Proteins\\prot_prot\\COVID19_GDDS_proteins_degrees_to_other_proteins_%s.tsv' % str(iteration), 'a+') as file_write:
# 			for protein in ordered_protein_counts_dict :
# 				if (sum(ordered_protein_counts_dict.values()) != 0) :
# 					file_write.write("%s\t%s\t%s\n" % (protein, ordered_protein_counts_dict[protein], ordered_protein_counts_dict[protein] / sum(ordered_protein_counts_dict.values())))
# 				else :
# 					file_write.write("%s\t%s\t%s\n" % (protein, ordered_protein_counts_dict[protein], ordered_protein_counts_dict[protein]))



# # # # ##################################################################
# # # # ############ COMPUTE DEGREES MEANS AND Z-SCORES ##################
# # # # ######## INCLUDING ABSENCE OF proteinS AS ZERO DEGREES ###########
# # # # ##################################################################


# # # For covid network, compute the mean and the standard deviation of proteins degrees and normalized degrees
# covid_proteins = []
# covid_protein_degree_dict = {}
# covid_protein_norm_degree_dict = {}

# total_protein_degree_dict = defaultdict(list)
# total_protein_normalized_degree_dict = defaultdict(list)

# with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\Further_Analysis\\Proteins\\COVID19_GDDS_proteins_degrees_to_other_proteins.tsv', 'r') as data :
# 	for line in data :
# 		line = line.split('\t')
# 		protein = line[0]
# 		covid_proteins.append(protein)
# 		protein_degree = int(line[1])
# 		protein_norm_degree = float(line[2].replace('\n', ''))
# 		if covid_protein_degree_dict.get(protein) == None :
# 			covid_protein_degree_dict[protein] = protein_degree
# 		else :
# 			print("warning : %s duplicate ?!" % protein)
# 		covid_protein_norm_degree_dict[protein] = protein_norm_degree

# 		if total_protein_degree_dict.get(protein) == None :
# 			total_protein_degree_dict[protein] = [protein_degree]
# 		else :
# 			total_protein_degree_dict[protein].append(protein_degree)

# 		if total_protein_normalized_degree_dict.get(protein) == None :
# 			total_protein_normalized_degree_dict[protein] = [protein_norm_degree]
# 		else :
# 			total_protein_normalized_degree_dict[protein].append(protein_norm_degree)

# mean_protein_deg = statistics.mean(covid_protein_degree_dict.values())
# mean_norm_protein_deg = statistics.mean(covid_protein_norm_degree_dict.values())
# sd_protein_deg = statistics.pstdev(covid_protein_degree_dict.values())
# sd_norm_protein_deg = statistics.pstdev(covid_protein_norm_degree_dict.values())
# print(mean_protein_deg, sd_protein_deg, mean_norm_protein_deg, sd_norm_protein_deg) #
# print(len(covid_proteins),len(covid_protein_degree_dict),len(covid_protein_norm_degree_dict)) #
# # 4.37417943107221 6.090315035936245 0.002188183807439825 0.003046680858397321
# # 457 457 457


# nodefiles_numbers = []
# edgefiles_numbers = []
# for filename in os.listdir("A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\random_networks\\Mock_networks_21Apr"):
# 	if 'nodes_' in filename :
# 		filename = filename.split(".csv")[0]
# 		nodefiles_numbers.append(int(filename[19:]))
# 	if 'edges_' in filename :
# 		filename = filename.split(".csv")[0]
# 		edgefiles_numbers.append(int(filename[19:]))

# # log_write.write("len(nodefiles_numbers) = %s\n" % len(nodefiles_numbers)) #2051
# # log_write.write("len(edgefiles_numbers) = %s\n" % len(edgefiles_numbers)) #2051
# # log_write.write("difference between edgefiles_numbers and nodefiles_numbers = %s\n" % len(set(nodefiles_numbers).difference(set(edgefiles_numbers)))) #0
# print("len(nodefiles_numbers) = %s\n" % len(nodefiles_numbers)) #2051
# print("len(edgefiles_numbers) = %s\n" % len(edgefiles_numbers)) #2051
# print("difference between edgefiles_numbers and nodefiles_numbers = %s\n" % len(set(nodefiles_numbers).difference(set(edgefiles_numbers)))) #0

# # nodefiles_numbers = sorted(nodefiles_numbers)
# # edgefiles_numbers = sorted(edgefiles_numbers)


# for iteration in nodefiles_numbers[1:] :
# 	# log_write.write("\n########################## ITERATION %s ###########################\n" % iteration)
# 	print("\n########################## ITERATION %s ###########################\n" % iteration)

# 	with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\Further_Analysis\\Proteins\\prot_prot\\COVID19_GDDS_proteins_degrees_to_other_proteins_%s.tsv' % str(iteration), 'r') as data: # Open the file  
# 		for line in data :
# 			line = line.split('\t')
# 			protein = line[0]
# 			degree = float(line[1])
# 			normalized_degree = float(line[2].replace('\n', ''))

# 			if total_protein_degree_dict.get(protein) == None :
# 				total_protein_degree_dict[protein] = [degree]
# 			else :
# 				total_protein_degree_dict[protein].append(degree)

# 			if total_protein_normalized_degree_dict.get(protein) == None :
# 				total_protein_normalized_degree_dict[protein] = [normalized_degree]
# 			else :
# 				total_protein_normalized_degree_dict[protein].append(normalized_degree)

# print(len(total_protein_degree_dict)) # 19928

# print(len(total_protein_normalized_degree_dict)) # 19928

# ### make the files with all calculations
# with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\Further_Analysis\\Proteins\\COVID19_GDDS_Results_Proteins_to_other_proteins-withZeros.tsv', 'a+') as file_write:
# 	file_write.write("protein\t# of nets it is found in\tmin-degree\tmax-degree\tmean-degree\tmedian-degree\tsd-degree\n")
# 	for protein in total_protein_degree_dict :
# 		nets_number = len(total_protein_degree_dict[protein])
# 		min_deg = min(total_protein_degree_dict[protein])
# 		max_deg = max(total_protein_degree_dict[protein])
# 		mean_deg = statistics.mean(total_protein_degree_dict[protein])
# 		median_deg = statistics.median(total_protein_degree_dict[protein])
# 		sd_deg = statistics.pstdev(total_protein_degree_dict[protein])
# 		file_write.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (protein, nets_number, min_deg, max_deg, mean_deg, median_deg, sd_deg))


# with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\Further_Analysis\\Proteins\\COVID19_GDDS_Results_Proteins_to_other_proteins_Relative-withZeros.tsv', 'a+') as file_write:
# 	file_write.write("protein\t# of nets it is found in\tmin-degree-rel\tmax-degree-rel\tmean-degree-rel\tmedian-degree-rel\tsd-degree-rel\n")
# 	for protein in total_protein_normalized_degree_dict :
# 		nets_number = len(total_protein_normalized_degree_dict[protein])
# 		min_deg = min(total_protein_normalized_degree_dict[protein])
# 		max_deg = max(total_protein_normalized_degree_dict[protein])
# 		mean_deg = statistics.mean(total_protein_normalized_degree_dict[protein])
# 		median_deg = statistics.median(total_protein_normalized_degree_dict[protein])
# 		sd_deg = statistics.pstdev(total_protein_normalized_degree_dict[protein])
# 		file_write.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (protein, nets_number, min_deg, max_deg, mean_deg, median_deg, sd_deg))

# ### compute Z-score for the covid proteins
# # For each protein found in covid network, compute 2 Z-scores.

# with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\Further_Analysis\\Proteins\\COVID19_GDDS_Results_ForCovidProteins_to_other_proteins-withZeros.tsv', 'a+') as file_write:
# 	file_write.write("covid protein\t# of nets it is found in\tmin-degree\tmax-degree\tmean-degree\tmedian-degree\tsd-degree\tcovid-degree\tZ-score\n")
# 	for protein in covid_proteins :
# 		nets_number = len(total_protein_degree_dict[protein])
# 		min_deg = min(total_protein_degree_dict[protein])
# 		max_deg = max(total_protein_degree_dict[protein])
# 		mean_deg = statistics.mean(total_protein_degree_dict[protein])
# 		median_deg = statistics.median(total_protein_degree_dict[protein])
# 		sd_deg = statistics.pstdev(total_protein_degree_dict[protein])
# 		covid_deg = covid_protein_degree_dict[protein]
# 		if sd_deg != 0 :
# 			z_score = ( covid_deg - mean_deg) / sd_deg
# 		else : 
# 			z_score = 'NaN'
# 		file_write.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (protein, nets_number, min_deg, max_deg, mean_deg, median_deg, sd_deg, covid_deg, z_score))


# with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\Further_Analysis\\Proteins\\COVID19_GDDS_Results_ForCovidProteins_to_other_proteins_Relative-withZeros.tsv', 'a+') as file_write:
# 	file_write.write("covid protein\t# of nets it is found in\tmin-degree-rel\tmax-degree-rel\tmean-degree-rel\tmedian-degree-rel\tsd-degree-rel\tcovid-degree-rel\tZ-score-rel\n")
# 	for protein in covid_proteins :
# 		nets_number = len(total_protein_normalized_degree_dict[protein])
# 		min_deg = min(total_protein_normalized_degree_dict[protein])
# 		max_deg = max(total_protein_normalized_degree_dict[protein])
# 		mean_deg = statistics.mean(total_protein_normalized_degree_dict[protein])
# 		median_deg = statistics.median(total_protein_normalized_degree_dict[protein])
# 		sd_deg = statistics.pstdev(total_protein_normalized_degree_dict[protein])
# 		covid_deg = covid_protein_norm_degree_dict[protein]
# 		if sd_deg != 0 :
# 			z_score = ( covid_deg - mean_deg) / sd_deg
# 		else : 
# 			z_score = 'NaN'
# 		file_write.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (protein, nets_number, min_deg, max_deg, mean_deg, median_deg, sd_deg, covid_deg, z_score))










stop = timeit.default_timer()
print(stop - start)
 
protein-symptom-degrees-directed-rankings-corrected.py,
#!/usr/bin/python
# -*- coding: utf-8 -*-
# jupyter notebook --browser='C:/Program Files (x86)/Google/Chrome/Application/chrome.exe %s'


import timeit
start = timeit.default_timer()

import os
import csv
from operator import itemgetter
import networkx as nx
from networkx.algorithms import community
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
from collections import Counter
from collections import defaultdict
from collections import OrderedDict
import statistics 
import re
import gzip





#####################################################################################################
#########################  Compute network Protein degree distribution #################################
##########################################  based on edges to symptoms #####################################
#####################################################################################################

with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\COVID19_GDDS_nodes.csv', 'r') as nodecsv: # Open the file                       
	nodereader = csv.reader(nodecsv) # Read the csv  
	nodes = [n for n in nodereader][1:]                     
	node_names = [n[0] for n in nodes] # Get a list of only the node names 

with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\COVID19_GDDS_edges.csv', 'r') as edgecsv: # Open the file
	edgereader = csv.reader(edgecsv) # Read the csv     
	edges = [tuple(e) for e in edgereader][1:] # Retrieve the data

G = nx.Graph()
G.add_nodes_from(node_names)
G.add_edges_from(edges)

descr_dict = {}
for node in nodes: 
	descr_dict[node[0]] = node[2]

print(nx.info(G))
# Number of nodes: 16007
# Number of edges: 53917
# Average degree:   6.7367

# print(nx.degree(G)) # this is the list containing as many tuples as nodes, indicating each node's degree.

average_degree = statistics.mean([ tpl[1] for tpl in nx.degree(G) ])
# print(average_degree) # 6.736677703504717

proteins_degrees = nx.degree(G, nbunch=[n for n in G.nodes if descr_dict[n] == 'Human PPI (target)'])
# print("degree of protein nodes only", proteins_degrees)

average_degree_of_proteins = statistics.mean([ tpl[1] for tpl in proteins_degrees ])
print("average degree of proteins nodes", average_degree_of_proteins) # 

# Number of nodes: 16007
# Number of edges: 53917
# Average degree:   6.7367
# average degree of proteins nodes 115.66739606126914


# # # # #####################################################################################################
# # # # ################################# CALCULATION ON EDGES WITH FILTERING PROT-sym #####################
# # # # #####################################################################################################

# now we keep only the edges of proteins that interact with other proteins
prots_linked_to_XXX = []
for u,v,c in G.edges(data=True) :	
	if (descr_dict[u] == 'Human PPI (target)' and descr_dict[v] == 'Symptom') :
		prots_linked_to_XXX.append(u)
print(len(prots_linked_to_XXX), len(set(prots_linked_to_XXX))) # # 10178 376


# ({'Drug': 5703, 'Disease': 4176, 'GO': 3487, 'Symptom': 2157, 'Human PPI (target)': 457, 'Viral Gene': 27})

# print(sorted(prots_linked_to_XXX))

protein_counts = sorted(prots_linked_to_XXX)

protein_counts_dict = defaultdict( int )
for protein in protein_counts:
    protein_counts_dict[protein] += 1

my_proteins_nodes = [n for n in G.nodes if descr_dict[n] == 'Human PPI (target)']
for protein in my_proteins_nodes :
	if protein_counts_dict.get(protein) == None :
		# print("%s is missing in the linked proteins : Adding it with a frequency of 0.\n" % protein)
		protein_counts_dict[protein]= 0

average_degree_of_my_proteins_nodes = statistics.mean([ protein_counts_dict[protein] for protein in protein_counts_dict ])
print("average_degree_of_my_proteins_nodes", average_degree_of_my_proteins_nodes)

max_degree_of_my_proteins_nodes = max([ protein_counts_dict[protein] for protein in protein_counts_dict ])
print("max_degree_of_my_proteins_nodes", max_degree_of_my_proteins_nodes)

# average_degree_of_my_proteins_nodes 22.271334792122538
# max_degree_of_my_proteins_nodes 347





ordered_protein_counts_dict = OrderedDict(sorted(protein_counts_dict.items()))


# with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\Further_Analysis\\Proteins\\COVID19_GDDS_proteins_degrees_to_symptoms.tsv', 'a+') as file_write:
# 	for protein in ordered_protein_counts_dict :
# 		file_write.write("%s\t%s\t%s\n" % (protein, ordered_protein_counts_dict[protein], ordered_protein_counts_dict[protein] / sum(ordered_protein_counts_dict.values())))


# # # ##################################################################
# # # ############ SAME ANALYSIS FOR RANDOM NETWORKS ###################
# # # # ##################################################################

# with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\Further_Analysis\\Proteins\\Log_degree_analysis_on_random_networks_for_proteins_symptoms.txt', 'a+') as log_write:

# 	nodefiles_numbers = []
# 	edgefiles_numbers = []
# 	for filename in os.listdir("A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\random_networks\\Mock_networks_21Apr"):
# 		if 'nodes_' in filename :
# 			filename = filename.split(".csv")[0]
# 			nodefiles_numbers.append(int(filename[19:]))
# 		if 'edges_' in filename :
# 			filename = filename.split(".csv")[0]
# 			edgefiles_numbers.append(int(filename[19:]))

# 	log_write.write("len(nodefiles_numbers) = %s\n" % len(nodefiles_numbers)) #2051
# 	log_write.write("len(edgefiles_numbers) = %s\n" % len(edgefiles_numbers)) #2051
# 	log_write.write("difference between edgefiles_numbers and nodefiles_numbers = %s\n" % len(set(nodefiles_numbers).difference(set(edgefiles_numbers)))) #0

# 	nodefiles_numbers = sorted(nodefiles_numbers)
# 	edgefiles_numbers = sorted(edgefiles_numbers)

# 	networks_to_avoid = [3, 11, 22, 23, 33, 34, 38, 57, 72, 77, 123, 127, 148, 156, 157, 183, 205, 209, 210, 217, 232, 252, 268, 281, 291, 312, 341, 353, 362, 364, 413, 421, 433, 443, 451, 494, 525, 532, 569, 592, 629, 658, 662, 674, 676, 701, 705, 712, 736, 757, 788, 796, 825, 827, 838, 873, 876, 908, 910, 916, 917, 918, 938, 958, 965, 1088, 1125, 1128, 1145, 1146, 1165, 1171, 1183, 1185, 1192, 1195, 1211, 1244, 1254, 1257, 1300, 1353, 1354, 1364, 1372, 1394, 1395, 1401, 1410, 1427, 1454, 1477, 1497, 1517, 1532, 1549, 1551, 1575, 1578, 1611, 1624, 1629, 1637, 1650, 1655, 1697, 1700, 1713, 1736, 1745, 1762, 1790, 1800, 1815, 1824, 1827, 1866, 1872, 1886, 1887, 2053, 2072, 2106, 2115, 2128, 2160, 2217, 2343, 2383, 2388, 2413, 2458, 2465]
# 	nodefiles_numbers = [item for item in nodefiles_numbers if item not in networks_to_avoid]

# 	for iteration in nodefiles_numbers[1:] :
# 		print("\n########################## ITERATION %s ###########################\n" % iteration)
# 		log_write.write("\n########################## ITERATION %s ###########################\n" % iteration)
	
# 		with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\random_networks\\Mock_networks_21Apr\\COVID19_GDDS_nodes_%s.csv' % str(iteration), 'r') as nodecsv: # Open the file                   
# 			nodereader = csv.reader(nodecsv) # Read the csv  
# 			nodes = [n for n in nodereader][1:]                     
# 			node_names = [n[0] for n in nodes] # Get a list of only the node names 

# 		with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\random_networks\\Mock_networks_21Apr\\COVID19_GDDS_edges_%s.csv' % str(iteration), 'r') as edgecsv: # Open the file
# 			edgereader = csv.reader(edgecsv) # Read the csv     
# 			edges = [tuple(e) for e in edgereader][1:] # Retrieve the data

# 			log_write.write("len(node_names) : %s\n" % len(node_names)) #16007 --> 13448
# 			log_write.write("len(edges) : %s\n" % len(edges)) #53917 --> 37835

# 		G = nx.Graph()
# 		G.add_nodes_from(node_names)
# 		G.add_edges_from(edges)

# 		descr_dict = {}
# 		description_set=set()
# 		for node in nodes: 
# 			descr_dict[node[0]] = node[2]
# 			description_set.add(node[2])

# 		# print(description_set) #{'Symptom', 'Human PPI (target)', 'GO', 'Disease', 'Drug'}
# 		# recounted = Counter(list(descr_dict.values()))
# 		# print("%s\n" % recounted)

# 		# print(nx.info(G))
# 		log_write.write("\nRANDOM GRAPH %s\n" % iteration)
# 		log_write.write("%s\n" % nx.info(G))

# 		#####################################################################################################
# 		################################# GET EDGES CONNECTING proteinS TO FIRST ORDER PROTEINS ################
# 		#####################################################################################################

# 		prots_linked_to_XXX = []
# 		for u,v,c in G.edges(data=True) :	
# 			if (descr_dict[u] == 'Human PPI (target)' and descr_dict[v] == 'Symptom') :
# 				prots_linked_to_XXX.append(u)
		
# 		protein_counts = sorted(prots_linked_to_XXX)

# 		protein_counts_dict = defaultdict( int )
# 		for protein in protein_counts:
# 		    protein_counts_dict[protein] += 1

# 		my_proteins_nodes = [n for n in G.nodes if descr_dict[n] == 'Human PPI (target)']
# 		for protein in my_proteins_nodes :
# 			if protein_counts_dict.get(protein) == None :
# 				# print("%s is missing in the linked proteins : Adding it with a frequency of 0.\n" % protein)
# 				protein_counts_dict[protein]= 0

# 		average_degree_of_my_proteins_nodes = statistics.mean([ protein_counts_dict[protein] for protein in protein_counts_dict ])
# 		log_write.write("Average degree of proteins linked to symptom(s) : %s\n" % average_degree_of_my_proteins_nodes)

# 		ordered_protein_counts_dict = OrderedDict(sorted(protein_counts_dict.items()))

# 		with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\Further_Analysis\\Proteins\\prot-sym\\COVID19_GDDS_proteins_degrees_to_symptoms_%s.tsv' % str(iteration), 'a+') as file_write:
# 			for protein in ordered_protein_counts_dict :
# 				if (sum(ordered_protein_counts_dict.values()) != 0) :
# 					file_write.write("%s\t%s\t%s\n" % (protein, ordered_protein_counts_dict[protein], ordered_protein_counts_dict[protein] / sum(ordered_protein_counts_dict.values())))
# 				else :
# 					file_write.write("%s\t%s\t%s\n" % (protein, ordered_protein_counts_dict[protein], ordered_protein_counts_dict[protein]))



# # # # ##################################################################
# # # # ############ COMPUTE DEGREES MEANS AND Z-SCORES ##################
# # # # ######## INCLUDING ABSENCE OF symptoms AS ZERO DEGREES ###########
# # # # ##################################################################


# # For covid network, compute the mean and the standard deviation of proteins degrees and normalized degrees
covid_proteins = []
covid_protein_degree_dict = {}
covid_protein_norm_degree_dict = {}

total_protein_degree_dict = defaultdict(list)
total_protein_normalized_degree_dict = defaultdict(list)

with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\Further_Analysis\\Proteins\\COVID19_GDDS_proteins_degrees_to_symptoms.tsv', 'r') as data :
	for line in data :
		line = line.split('\t')
		protein = line[0]
		covid_proteins.append(protein)
		protein_degree = int(line[1])
		protein_norm_degree = float(line[2].replace('\n', ''))
		if covid_protein_degree_dict.get(protein) == None :
			covid_protein_degree_dict[protein] = protein_degree
		else :
			print("warning : %s duplicate ?!" % protein)
		covid_protein_norm_degree_dict[protein] = protein_norm_degree

		if total_protein_degree_dict.get(protein) == None :
			total_protein_degree_dict[protein] = [protein_degree]
		else :
			total_protein_degree_dict[protein].append(protein_degree)

		if total_protein_normalized_degree_dict.get(protein) == None :
			total_protein_normalized_degree_dict[protein] = [protein_norm_degree]
		else :
			total_protein_normalized_degree_dict[protein].append(protein_norm_degree)

mean_protein_deg = statistics.mean(covid_protein_degree_dict.values())
mean_norm_protein_deg = statistics.mean(covid_protein_norm_degree_dict.values())
sd_protein_deg = statistics.pstdev(covid_protein_degree_dict.values())
sd_norm_protein_deg = statistics.pstdev(covid_protein_norm_degree_dict.values())
print(mean_protein_deg, sd_protein_deg, mean_norm_protein_deg, sd_norm_protein_deg) #
print(len(covid_proteins),len(covid_protein_degree_dict),len(covid_protein_norm_degree_dict)) #
# 22.271334792122538 39.8090726298042 0.002188183807439825 0.003911286365671468
# 457 457 457





nodefiles_numbers = []
edgefiles_numbers = []
for filename in os.listdir("A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\random_networks\\Mock_networks_21Apr"):
	if 'nodes_' in filename :
		filename = filename.split(".csv")[0]
		nodefiles_numbers.append(int(filename[19:]))
	if 'edges_' in filename :
		filename = filename.split(".csv")[0]
		edgefiles_numbers.append(int(filename[19:]))

# log_write.write("len(nodefiles_numbers) = %s\n" % len(nodefiles_numbers)) #2051
# log_write.write("len(edgefiles_numbers) = %s\n" % len(edgefiles_numbers)) #2051
# log_write.write("difference between edgefiles_numbers and nodefiles_numbers = %s\n" % len(set(nodefiles_numbers).difference(set(edgefiles_numbers)))) #0
print("len(nodefiles_numbers) = %s\n" % len(nodefiles_numbers)) #2051
print("len(edgefiles_numbers) = %s\n" % len(edgefiles_numbers)) #2051
print("difference between edgefiles_numbers and nodefiles_numbers = %s\n" % len(set(nodefiles_numbers).difference(set(edgefiles_numbers)))) #0

# nodefiles_numbers = sorted(nodefiles_numbers)
# edgefiles_numbers = sorted(edgefiles_numbers)

networks_to_avoid = [3, 11, 22, 23, 33, 34, 38, 57, 72, 77, 123, 127, 148, 156, 157, 183, 205, 209, 210, 217, 232, 252, 268, 281, 291, 312, 341, 353, 362, 364, 413, 421, 433, 443, 451, 494, 525, 532, 569, 592, 629, 658, 662, 674, 676, 701, 705, 712, 736, 757, 788, 796, 825, 827, 838, 873, 876, 908, 910, 916, 917, 918, 938, 958, 965, 1088, 1125, 1128, 1145, 1146, 1165, 1171, 1183, 1185, 1192, 1195, 1211, 1244, 1254, 1257, 1300, 1353, 1354, 1364, 1372, 1394, 1395, 1401, 1410, 1427, 1454, 1477, 1497, 1517, 1532, 1549, 1551, 1575, 1578, 1611, 1624, 1629, 1637, 1650, 1655, 1697, 1700, 1713, 1736, 1745, 1762, 1790, 1800, 1815, 1824, 1827, 1866, 1872, 1886, 1887, 2053, 2072, 2106, 2115, 2128, 2160, 2217, 2343, 2383, 2388, 2413, 2458, 2465]
nodefiles_numbers = [item for item in nodefiles_numbers if item not in networks_to_avoid]

for iteration in nodefiles_numbers[1:] :
	# log_write.write("\n########################## ITERATION %s ###########################\n" % iteration)
	print("\n########################## ITERATION %s ###########################\n" % iteration)

	with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\Further_Analysis\\Proteins\\prot-sym\\COVID19_GDDS_proteins_degrees_to_symptoms_%s.tsv' % str(iteration), 'r') as data: # Open the file  
		for line in data :
			line = line.split('\t')
			protein = line[0]
			degree = float(line[1])
			normalized_degree = float(line[2].replace('\n', ''))

			if total_protein_degree_dict.get(protein) == None :
				total_protein_degree_dict[protein] = [degree]
			else :
				total_protein_degree_dict[protein].append(degree)

			if total_protein_normalized_degree_dict.get(protein) == None :
				total_protein_normalized_degree_dict[protein] = [normalized_degree]
			else :
				total_protein_normalized_degree_dict[protein].append(normalized_degree)

print(len(total_protein_degree_dict)) # 19928

print(len(total_protein_normalized_degree_dict)) # 19928

### make the files with all calculations
with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\Further_Analysis\\Proteins\\COVID19_GDDS_Results_Proteins_to_symptoms-withZeros.tsv', 'a+') as file_write:
	file_write.write("protein\t# of nets it is found in\tmin-degree\tmax-degree\tmean-degree\tmedian-degree\tsd-degree\n")
	for protein in total_protein_degree_dict :
		nets_number = len(total_protein_degree_dict[protein])
		min_deg = min(total_protein_degree_dict[protein])
		max_deg = max(total_protein_degree_dict[protein])
		mean_deg = statistics.mean(total_protein_degree_dict[protein])
		median_deg = statistics.median(total_protein_degree_dict[protein])
		sd_deg = statistics.pstdev(total_protein_degree_dict[protein])
		file_write.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (protein, nets_number, min_deg, max_deg, mean_deg, median_deg, sd_deg))


with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\Further_Analysis\\Proteins\\COVID19_GDDS_Results_Proteins_to_symptoms_Relative-withZeros.tsv', 'a+') as file_write:
	file_write.write("protein\t# of nets it is found in\tmin-degree-rel\tmax-degree-rel\tmean-degree-rel\tmedian-degree-rel\tsd-degree-rel\n")
	for protein in total_protein_normalized_degree_dict :
		nets_number = len(total_protein_normalized_degree_dict[protein])
		min_deg = min(total_protein_normalized_degree_dict[protein])
		max_deg = max(total_protein_normalized_degree_dict[protein])
		mean_deg = statistics.mean(total_protein_normalized_degree_dict[protein])
		median_deg = statistics.median(total_protein_normalized_degree_dict[protein])
		sd_deg = statistics.pstdev(total_protein_normalized_degree_dict[protein])
		file_write.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (protein, nets_number, min_deg, max_deg, mean_deg, median_deg, sd_deg))

### compute Z-score for the covid proteins
# For each protein found in covid network, compute 2 Z-scores.

with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\Further_Analysis\\Proteins\\COVID19_GDDS_Results_ForCovidProteins_to_symptoms-withZeros.tsv', 'a+') as file_write:
	file_write.write("covid protein\t# of nets it is found in\tmin-degree\tmax-degree\tmean-degree\tmedian-degree\tsd-degree\tcovid-degree\tZ-score\n")
	for protein in covid_proteins :
		nets_number = len(total_protein_degree_dict[protein])
		min_deg = min(total_protein_degree_dict[protein])
		max_deg = max(total_protein_degree_dict[protein])
		mean_deg = statistics.mean(total_protein_degree_dict[protein])
		median_deg = statistics.median(total_protein_degree_dict[protein])
		sd_deg = statistics.pstdev(total_protein_degree_dict[protein])
		covid_deg = covid_protein_degree_dict[protein]
		if sd_deg != 0 :
			z_score = ( covid_deg - mean_deg) / sd_deg
		else : 
			z_score = 'NaN'
		file_write.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (protein, nets_number, min_deg, max_deg, mean_deg, median_deg, sd_deg, covid_deg, z_score))


with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\Further_Analysis\\Proteins\\COVID19_GDDS_Results_ForCovidProteins_to_symptoms_Relative-withZeros.tsv', 'a+') as file_write:
	file_write.write("covid protein\t# of nets it is found in\tmin-degree-rel\tmax-degree-rel\tmean-degree-rel\tmedian-degree-rel\tsd-degree-rel\tcovid-degree-rel\tZ-score-rel\n")
	for protein in covid_proteins :
		nets_number = len(total_protein_normalized_degree_dict[protein])
		min_deg = min(total_protein_normalized_degree_dict[protein])
		max_deg = max(total_protein_normalized_degree_dict[protein])
		mean_deg = statistics.mean(total_protein_normalized_degree_dict[protein])
		median_deg = statistics.median(total_protein_normalized_degree_dict[protein])
		sd_deg = statistics.pstdev(total_protein_normalized_degree_dict[protein])
		covid_deg = covid_protein_norm_degree_dict[protein]
		if sd_deg != 0 :
			z_score = ( covid_deg - mean_deg) / sd_deg
		else : 
			z_score = 'NaN'
		file_write.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (protein, nets_number, min_deg, max_deg, mean_deg, median_deg, sd_deg, covid_deg, z_score))










stop = timeit.default_timer()
print(stop - start)
  
protein-target-degrees-UNdirected-rankings-corrected-generic-script.py,
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
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
from collections import Counter
from collections import defaultdict
from collections import OrderedDict
import statistics 
import re
import gzip


# target = sys.argv[1]
# print(target)
# print(type(target))

dirname = "full_analysis_for_proteins_undirected"
os.mkdir(dirname)

# with open log

#####################################################################################################
#########################  Compute network Protein degree distribution #################################
##########################################  based on edges to protein #####################################
#####################################################################################################

with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\COVID19_GDDS_nodes.csv', 'r') as nodecsv: # Open the file                       
	nodereader = csv.reader(nodecsv) # Read the csv  
	nodes = [n for n in nodereader][1:]                     
	node_names = [n[0] for n in nodes] # Get a list of only the node names 

with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\COVID19_GDDS_edges.csv', 'r') as edgecsv: # Open the file
	edgereader = csv.reader(edgecsv) # Read the csv     
	edges = [tuple(e) for e in edgereader][1:] # Retrieve the data

G_covid = nx.Graph()
G_covid.add_nodes_from(node_names)
G_covid.add_edges_from(edges)

descr_dict_covid = {}
for node in nodes: 
	descr_dict_covid[node[0]] = node[2]

print(nx.info(G_covid))
# Number of nodes: 16007
# Number of edges: 53917
# Average degree:   6.7367

# print(nx.degree(G_covid)) # this is the list containing as many tuples as nodes, indicating each node's degree.

average_degree = statistics.mean([ tpl[1] for tpl in nx.degree(G_covid) ])
# print(average_degree) # 6.736677703504717

proteins_degrees = nx.degree(G_covid, nbunch=[n for n in G_covid.nodes if descr_dict_covid[n] == 'Human PPI (target)'])
# print("degree of protein nodes only", proteins_degrees)

average_degree_of_proteins = statistics.mean([ tpl[1] for tpl in proteins_degrees ])
print("average degree of proteins nodes", average_degree_of_proteins) # 

# Number of nodes: 16007
# Number of edges: 53917
# Average degree:   6.7367
# average degree of proteins nodes 115.66739606126914


# # # # #####################################################################################################
# # # # ################################# CALCULATION ON EDGES WITH FILTERING PROT-target #####################
# # # # #####################################################################################################

target_list = ['GO', 'Symptom', 'Disease','Drug', 'Human PPI (target)']

for target in target_list :
	print("Searching Human PPI (target) - %s links\n" % target)

	subfoldername = 'full_analysis_for_proteins_undirected/prot-%s' % (target)
	os.mkdir(subfoldername)
	# now we keep only the edges of proteins that interact with targets
	prots_linked_to_targets = []
	for u,v,c in G_covid.edges(data=True) :	
		if (descr_dict_covid[u] == 'Human PPI (target)' and descr_dict_covid[v] == target) :
			prots_linked_to_targets.append(u)
		if (descr_dict_covid[v] == 'Human PPI (target)' and descr_dict_covid[u] == target) :
			prots_linked_to_targets.append(v)	
	print("len(prots_linked_to_targets), len(set(prots_linked_to_targets))", len(prots_linked_to_targets), len(set(prots_linked_to_targets))) # # 9208 435

	protein_counts = sorted(prots_linked_to_targets)

	protein_counts_dict = defaultdict( int )
	for protein in protein_counts:
	    protein_counts_dict[protein] += 1

	my_proteins_nodes = [n for n in G_covid.nodes if descr_dict_covid[n] == 'Human PPI (target)']
	for protein in my_proteins_nodes :
		if protein_counts_dict.get(protein) == None :
			protein_counts_dict[protein]= 0

	average_degree_of_my_proteins_nodes = statistics.mean([ protein_counts_dict[protein] for protein in protein_counts_dict ])
	print("average_degree_of_my_proteins_nodes", average_degree_of_my_proteins_nodes)

	max_degree_of_my_proteins_nodes = max([ protein_counts_dict[protein] for protein in protein_counts_dict ])
	print("max_degree_of_my_proteins_nodes", max_degree_of_my_proteins_nodes)

	ordered_protein_counts_dict = OrderedDict(sorted(protein_counts_dict.items()))


	with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\Further_Analysis\\Proteins\\full_analysis_for_proteins_undirected\\COVID19_GDDS_proteins_degrees_to_%s.tsv' % target, 'a+') as file_write:
		for protein in ordered_protein_counts_dict :
			file_write.write("%s\t%s\t%s\n" % (protein, ordered_protein_counts_dict[protein], ordered_protein_counts_dict[protein] / sum(ordered_protein_counts_dict.values())))


	# # # # ##################################################################
	# # # # ############ SAME ANALYSIS FOR RANDOM NETWORKS ###################
	# # # # # ##################################################################

	with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\Further_Analysis\\Proteins\\full_analysis_for_proteins_undirected\\Log_degree_analysis_on_random_networks_for_proteins_%s.txt' % target, 'a+') as log_write:
		log_write.write("In covid net : len(prots_linked_to_targets) and len(set(prots_linked_to_targets)) are : %s and %s \n" % (len(prots_linked_to_targets), len(set(prots_linked_to_targets))) ) # # 9208 435
		log_write.write("In covid net : average_degree_of_my_proteins_nodes = %s\n" % average_degree_of_my_proteins_nodes)
		log_write.write("In covid net : max_degree_of_my_proteins_nodes = %s\n" % max_degree_of_my_proteins_nodes)

		nodefiles_numbers = []
		edgefiles_numbers = []
		for filename in os.listdir("A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\random_networks\\Mock_networks_21Apr"):
			if 'nodes_' in filename :
				filename = filename.split(".csv")[0]
				nodefiles_numbers.append(int(filename[19:]))
			if 'edges_' in filename :
				filename = filename.split(".csv")[0]
				edgefiles_numbers.append(int(filename[19:]))

		# log_write.write("len(nodefiles_numbers) = %s\n" % len(nodefiles_numbers)) #2051
		# log_write.write("len(edgefiles_numbers) = %s\n" % len(edgefiles_numbers)) #2051
		# log_write.write("difference between edgefiles_numbers and nodefiles_numbers = %s\n" % len(set(nodefiles_numbers).difference(set(edgefiles_numbers)))) #0

		nodefiles_numbers = sorted(nodefiles_numbers)
		edgefiles_numbers = sorted(edgefiles_numbers)

		networks_to_avoid = [3, 11, 22, 23, 33, 34, 38, 57, 72, 77, 123, 127, 148, 156, 157, 183, 205, 209, 210, 217, 232, 252, 268, 281, 291, 312, 341, 353, 362, 364, 413, 421, 433, 443, 451, 494, 525, 532, 569, 592, 629, 658, 662, 674, 676, 701, 705, 712, 736, 757, 788, 796, 825, 827, 838, 873, 876, 908, 910, 916, 917, 918, 938, 958, 965, 1088, 1125, 1128, 1145, 1146, 1165, 1171, 1183, 1185, 1192, 1195, 1211, 1244, 1254, 1257, 1300, 1353, 1354, 1364, 1372, 1394, 1395, 1401, 1410, 1427, 1454, 1477, 1497, 1517, 1532, 1549, 1551, 1575, 1578, 1611, 1624, 1629, 1637, 1650, 1655, 1697, 1700, 1713, 1736, 1745, 1762, 1790, 1800, 1815, 1824, 1827, 1866, 1872, 1886, 1887, 2053, 2072, 2106, 2115, 2128, 2160, 2217, 2343, 2383, 2388, 2413, 2458, 2465]
		nodefiles_numbers = [item for item in nodefiles_numbers if item not in networks_to_avoid]

		for iteration in nodefiles_numbers[1:] :
			print("\n########################## ITERATION %s -- %s ###########################\n" % (iteration, target))
			log_write.write("\n########################## ITERATION %s ###########################\n" % iteration)
		
			with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\random_networks\\Mock_networks_21Apr\\COVID19_GDDS_nodes_%s.csv' % str(iteration), 'r') as nodecsv: # Open the file                   
				nodereader = csv.reader(nodecsv) # Read the csv  
				nodes = [n for n in nodereader][1:]                     
				node_names = [n[0] for n in nodes] # Get a list of only the node names 

			with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\random_networks\\Mock_networks_21Apr\\COVID19_GDDS_edges_%s.csv' % str(iteration), 'r') as edgecsv: # Open the file
				edgereader = csv.reader(edgecsv) # Read the csv     
				edges = [tuple(e) for e in edgereader][1:] # Retrieve the data

				log_write.write("len(node_names) : %s\n" % len(node_names)) #16007 --> 13448
				log_write.write("len(edges) : %s\n" % len(edges)) #53917 --> 37835

			G_random = nx.Graph()
			G_random.add_nodes_from(node_names)
			G_random.add_edges_from(edges)

			descr_dict_random = {}
			description_random_set=set()
			for node in nodes: 
				descr_dict_random[node[0]] = node[2]
				description_random_set.add(node[2])

			# print(nx.info(G_random))
			log_write.write("\nRANDOM GRAPH %s\n" % iteration)
			log_write.write("%s\n" % nx.info(G_random))

			#####################################################################################################
			################################# GET EDGES CONNECTING proteinS TO target ################
			#####################################################################################################

			prots_linked_to_targets = []
			for u,v,c in G_random.edges(data=True) :	
				if (descr_dict_random[u] == 'Human PPI (target)' and descr_dict_random[v] == target) :
					prots_linked_to_targets.append(u)
				if (descr_dict_random[v] == 'Human PPI (target)' and descr_dict_random[u] == target) :
					prots_linked_to_targets.append(v)		
			
			protein_counts = sorted(prots_linked_to_targets)

			protein_counts_dict = defaultdict( int )
			for protein in protein_counts:
			    protein_counts_dict[protein] += 1

			my_proteins_nodes = [n for n in G_random.nodes if descr_dict_random[n] == 'Human PPI (target)']
			for protein in my_proteins_nodes :
				if protein_counts_dict.get(protein) == None :
					# print("%s is missing in the linked proteins : Adding it with a frequency of 0.\n" % protein)
					protein_counts_dict[protein]= 0

			average_degree_of_my_proteins_nodes = statistics.mean([ protein_counts_dict[protein] for protein in protein_counts_dict ])
			log_write.write("Average degree of Human PPI (target) linked to %s : %s\n" % (target, average_degree_of_my_proteins_nodes))

			ordered_protein_counts_dict = OrderedDict(sorted(protein_counts_dict.items()))

			with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\Further_Analysis\\Proteins\\full_analysis_for_proteins_undirected\\prot-%s\\COVID19_GDDS_proteins_degrees_to_%s_%s.tsv' % (target, target, str(iteration)), 'a+') as file_write:
				for protein in ordered_protein_counts_dict :
					if (sum(ordered_protein_counts_dict.values()) != 0) :
						file_write.write("%s\t%s\t%s\n" % (protein, ordered_protein_counts_dict[protein], ordered_protein_counts_dict[protein] / sum(ordered_protein_counts_dict.values())))
					else :
						file_write.write("%s\t%s\t%s\n" % (protein, ordered_protein_counts_dict[protein], ordered_protein_counts_dict[protein]))



	# # # # ##################################################################
	# # # # ############ COMPUTE DEGREES MEANS AND Z-SCORES ##################
	# # # # ######## INCLUDING ABSENCE OF target AS ZERO DEGREES ###########
	# # # # ##################################################################


	# # For covid network, compute the mean and the standard deviation of proteins degrees and normalized degrees
	covid_proteins = []
	covid_protein_degree_dict = {}
	covid_protein_norm_degree_dict = {}

	total_protein_degree_dict = defaultdict(list)
	total_protein_normalized_degree_dict = defaultdict(list)

	with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\Further_Analysis\\Proteins\\full_analysis_for_proteins_undirected\\COVID19_GDDS_proteins_degrees_to_%s.tsv' % target, 'r') as data :
		for line in data :
			line = line.split('\t')
			protein = line[0]
			covid_proteins.append(protein)
			protein_degree = int(line[1])
			protein_norm_degree = float(line[2].replace('\n', ''))
			if covid_protein_degree_dict.get(protein) == None :
				covid_protein_degree_dict[protein] = protein_degree
			else :
				print("warning : %s duplicate ?!" % protein)
			covid_protein_norm_degree_dict[protein] = protein_norm_degree

			if total_protein_degree_dict.get(protein) == None :
				total_protein_degree_dict[protein] = [protein_degree]
			else :
				total_protein_degree_dict[protein].append(protein_degree)

			if total_protein_normalized_degree_dict.get(protein) == None :
				total_protein_normalized_degree_dict[protein] = [protein_norm_degree]
			else :
				total_protein_normalized_degree_dict[protein].append(protein_norm_degree)

	mean_protein_deg = statistics.mean(covid_protein_degree_dict.values())
	mean_norm_protein_deg = statistics.mean(covid_protein_norm_degree_dict.values())
	sd_protein_deg = statistics.pstdev(covid_protein_degree_dict.values())
	sd_norm_protein_deg = statistics.pstdev(covid_protein_norm_degree_dict.values())
	print(mean_protein_deg, sd_protein_deg, mean_norm_protein_deg, sd_norm_protein_deg) #
	print(len(covid_proteins),len(covid_protein_degree_dict),len(covid_protein_norm_degree_dict)) #
	# 20.14879649890591 18.9696222570052 0.002188183807439825 0.00206012405050013
	# 457 457 457






	nodefiles_numbers = []
	edgefiles_numbers = []
	for filename in os.listdir("A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\random_networks\\Mock_networks_21Apr"):
		if 'nodes_' in filename :
			filename = filename.split(".csv")[0]
			nodefiles_numbers.append(int(filename[19:]))
		if 'edges_' in filename :
			filename = filename.split(".csv")[0]
			edgefiles_numbers.append(int(filename[19:]))

	# log_write.write("len(nodefiles_numbers) = %s\n" % len(nodefiles_numbers)) #2051
	# log_write.write("len(edgefiles_numbers) = %s\n" % len(edgefiles_numbers)) #2051
	# log_write.write("difference between edgefiles_numbers and nodefiles_numbers = %s\n" % len(set(nodefiles_numbers).difference(set(edgefiles_numbers)))) #0
	print("len(nodefiles_numbers) = %s\n" % len(nodefiles_numbers)) #2051
	print("len(edgefiles_numbers) = %s\n" % len(edgefiles_numbers)) #2051
	print("difference between edgefiles_numbers and nodefiles_numbers = %s\n" % len(set(nodefiles_numbers).difference(set(edgefiles_numbers)))) #0

	# nodefiles_numbers = sorted(nodefiles_numbers)
	# edgefiles_numbers = sorted(edgefiles_numbers)

	networks_to_avoid = [3, 11, 22, 23, 33, 34, 38, 57, 72, 77, 123, 127, 148, 156, 157, 183, 205, 209, 210, 217, 232, 252, 268, 281, 291, 312, 341, 353, 362, 364, 413, 421, 433, 443, 451, 494, 525, 532, 569, 592, 629, 658, 662, 674, 676, 701, 705, 712, 736, 757, 788, 796, 825, 827, 838, 873, 876, 908, 910, 916, 917, 918, 938, 958, 965, 1088, 1125, 1128, 1145, 1146, 1165, 1171, 1183, 1185, 1192, 1195, 1211, 1244, 1254, 1257, 1300, 1353, 1354, 1364, 1372, 1394, 1395, 1401, 1410, 1427, 1454, 1477, 1497, 1517, 1532, 1549, 1551, 1575, 1578, 1611, 1624, 1629, 1637, 1650, 1655, 1697, 1700, 1713, 1736, 1745, 1762, 1790, 1800, 1815, 1824, 1827, 1866, 1872, 1886, 1887, 2053, 2072, 2106, 2115, 2128, 2160, 2217, 2343, 2383, 2388, 2413, 2458, 2465]
	nodefiles_numbers = [item for item in nodefiles_numbers if item not in networks_to_avoid]

	for iteration in nodefiles_numbers[1:] :
		# log_write.write("\n########################## ITERATION %s ###########################\n" % iteration)
		print("\n########################## ITERATION %s ###########################\n" % iteration)

		with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\Further_Analysis\\Proteins\\full_analysis_for_proteins_undirected\\prot-%s\\COVID19_GDDS_proteins_degrees_to_%s_%s.tsv' % (target, target, str(iteration)), 'r') as data: # Open the file  
			for line in data :
				line = line.split('\t')
				protein = line[0]
				degree = float(line[1])
				normalized_degree = float(line[2].replace('\n', ''))

				if total_protein_degree_dict.get(protein) == None :
					total_protein_degree_dict[protein] = [degree]
				else :
					total_protein_degree_dict[protein].append(degree)

				if total_protein_normalized_degree_dict.get(protein) == None :
					total_protein_normalized_degree_dict[protein] = [normalized_degree]
				else :
					total_protein_normalized_degree_dict[protein].append(normalized_degree)

	print(len(total_protein_degree_dict)) # 19928

	print(len(total_protein_normalized_degree_dict)) # 19928

	### make the files with all calculations
	with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\Further_Analysis\\Proteins\\full_analysis_for_proteins_undirected\\COVID19_GDDS_Results_Proteins_to_%s-withZeros.tsv' % target, 'a+') as file_write:
		file_write.write("protein\t# of nets it is found in\tmin-degree\tmax-degree\tmean-degree\tmedian-degree\tsd-degree\n")
		for protein in total_protein_degree_dict :
			nets_number = len(total_protein_degree_dict[protein])
			min_deg = min(total_protein_degree_dict[protein])
			max_deg = max(total_protein_degree_dict[protein])
			mean_deg = statistics.mean(total_protein_degree_dict[protein])
			median_deg = statistics.median(total_protein_degree_dict[protein])
			sd_deg = statistics.pstdev(total_protein_degree_dict[protein])
			file_write.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (protein, nets_number, min_deg, max_deg, mean_deg, median_deg, sd_deg))


	with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\Further_Analysis\\Proteins\\full_analysis_for_proteins_undirected\\COVID19_GDDS_Results_Proteins_to_%s_Relative-withZeros.tsv' % target, 'a+') as file_write:
		file_write.write("protein\t# of nets it is found in\tmin-degree-rel\tmax-degree-rel\tmean-degree-rel\tmedian-degree-rel\tsd-degree-rel\n")
		for protein in total_protein_normalized_degree_dict :
			nets_number = len(total_protein_normalized_degree_dict[protein])
			min_deg = min(total_protein_normalized_degree_dict[protein])
			max_deg = max(total_protein_normalized_degree_dict[protein])
			mean_deg = statistics.mean(total_protein_normalized_degree_dict[protein])
			median_deg = statistics.median(total_protein_normalized_degree_dict[protein])
			sd_deg = statistics.pstdev(total_protein_normalized_degree_dict[protein])
			file_write.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (protein, nets_number, min_deg, max_deg, mean_deg, median_deg, sd_deg))

	### compute Z-score for the covid proteins
	# For each protein found in covid network, compute 2 Z-scores.

	with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\Further_Analysis\\Proteins\\full_analysis_for_proteins_undirected\\COVID19_GDDS_Results_ForCovidProteins_to_%s-withZeros.tsv' % target, 'a+') as file_write:
		file_write.write("covid protein\t# of nets it is found in\tmin-degree\tmax-degree\tmean-degree\tmedian-degree\tsd-degree\tcovid-degree\tZ-score\n")
		for protein in covid_proteins :
			nets_number = len(total_protein_degree_dict[protein])
			min_deg = min(total_protein_degree_dict[protein])
			max_deg = max(total_protein_degree_dict[protein])
			mean_deg = statistics.mean(total_protein_degree_dict[protein])
			median_deg = statistics.median(total_protein_degree_dict[protein])
			sd_deg = statistics.pstdev(total_protein_degree_dict[protein])
			covid_deg = covid_protein_degree_dict[protein]
			if sd_deg != 0 :
				z_score = ( covid_deg - mean_deg) / sd_deg
			else : 
				z_score = 'NaN'
			file_write.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (protein, nets_number, min_deg, max_deg, mean_deg, median_deg, sd_deg, covid_deg, z_score))


	with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\Further_Analysis\\Proteins\\full_analysis_for_proteins_undirected\\COVID19_GDDS_Results_ForCovidProteins_to_%s_Relative-withZeros.tsv' % target, 'a+') as file_write:
		file_write.write("covid protein\t# of nets it is found in\tmin-degree-rel\tmax-degree-rel\tmean-degree-rel\tmedian-degree-rel\tsd-degree-rel\tcovid-degree-rel\tZ-score-rel\n")
		for protein in covid_proteins :
			nets_number = len(total_protein_normalized_degree_dict[protein])
			min_deg = min(total_protein_normalized_degree_dict[protein])
			max_deg = max(total_protein_normalized_degree_dict[protein])
			mean_deg = statistics.mean(total_protein_normalized_degree_dict[protein])
			median_deg = statistics.median(total_protein_normalized_degree_dict[protein])
			sd_deg = statistics.pstdev(total_protein_normalized_degree_dict[protein])
			covid_deg = covid_protein_norm_degree_dict[protein]
			if sd_deg != 0 :
				z_score = ( covid_deg - mean_deg) / sd_deg
			else : 
				z_score = 'NaN'
			file_write.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (protein, nets_number, min_deg, max_deg, mean_deg, median_deg, sd_deg, covid_deg, z_score))










stop = timeit.default_timer()
print(stop - start)
  
protein-target-degrees-directed-rankings-corrected-generic-script.py,
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
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
from collections import Counter
from collections import defaultdict
from collections import OrderedDict
import statistics 
import re
import gzip


# target = sys.argv[1]
# print(target)
# print(type(target))

dirname = "full_analysis_for_proteins"
os.mkdir(dirname)

# with open log

#####################################################################################################
#########################  Compute network Protein degree distribution #################################
##########################################  based on edges to protein #####################################
#####################################################################################################

with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\COVID19_GDDS_nodes.csv', 'r') as nodecsv: # Open the file                       
	nodereader = csv.reader(nodecsv) # Read the csv  
	nodes = [n for n in nodereader][1:]                     
	node_names = [n[0] for n in nodes] # Get a list of only the node names 

with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\COVID19_GDDS_edges.csv', 'r') as edgecsv: # Open the file
	edgereader = csv.reader(edgecsv) # Read the csv     
	edges = [tuple(e) for e in edgereader][1:] # Retrieve the data

G_covid = nx.Graph()
G_covid.add_nodes_from(node_names)
G_covid.add_edges_from(edges)

descr_dict_covid = {}
for node in nodes: 
	descr_dict_covid[node[0]] = node[2]

print(nx.info(G_covid))
# Number of nodes: 16007
# Number of edges: 53917
# Average degree:   6.7367

# print(nx.degree(G_covid)) # this is the list containing as many tuples as nodes, indicating each node's degree.

average_degree = statistics.mean([ tpl[1] for tpl in nx.degree(G_covid) ])
# print(average_degree) # 6.736677703504717

proteins_degrees = nx.degree(G_covid, nbunch=[n for n in G_covid.nodes if descr_dict_covid[n] == 'Human PPI (target)'])
# print("degree of protein nodes only", proteins_degrees)

average_degree_of_proteins = statistics.mean([ tpl[1] for tpl in proteins_degrees ])
print("average degree of proteins nodes", average_degree_of_proteins) # 

# Number of nodes: 16007
# Number of edges: 53917
# Average degree:   6.7367
# average degree of proteins nodes 115.66739606126914


# # # # #####################################################################################################
# # # # ################################# CALCULATION ON EDGES WITH FILTERING PROT-target #####################
# # # # #####################################################################################################

target_list = ['GO', 'Symptom', 'Disease','Drug', 'Human PPI (target)']

for target in target_list :
	print("Searching Human PPI (target) - %s links\n" % target)

	subfoldername = 'full_analysis_for_proteins/prot-%s' % (target)
	os.mkdir(subfoldername)
	# now we keep only the edges of proteins that interact with targets
	prots_linked_to_targets = []
	for u,v,c in G_covid.edges(data=True) :	
		if (descr_dict_covid[u] == 'Human PPI (target)' and descr_dict_covid[v] == target) :
			prots_linked_to_targets.append(u)
	print("len(prots_linked_to_targets), len(set(prots_linked_to_targets))", len(prots_linked_to_targets), len(set(prots_linked_to_targets))) # # 9208 435

	protein_counts = sorted(prots_linked_to_targets)

	protein_counts_dict = defaultdict( int )
	for protein in protein_counts:
	    protein_counts_dict[protein] += 1

	my_proteins_nodes = [n for n in G_covid.nodes if descr_dict_covid[n] == 'Human PPI (target)']
	for protein in my_proteins_nodes :
		if protein_counts_dict.get(protein) == None :
			protein_counts_dict[protein]= 0

	average_degree_of_my_proteins_nodes = statistics.mean([ protein_counts_dict[protein] for protein in protein_counts_dict ])
	print("average_degree_of_my_proteins_nodes", average_degree_of_my_proteins_nodes)

	max_degree_of_my_proteins_nodes = max([ protein_counts_dict[protein] for protein in protein_counts_dict ])
	print("max_degree_of_my_proteins_nodes", max_degree_of_my_proteins_nodes)

	ordered_protein_counts_dict = OrderedDict(sorted(protein_counts_dict.items()))


	with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\Further_Analysis\\Proteins\\full_analysis_for_proteins\\COVID19_GDDS_proteins_degrees_to_%s.tsv' % target, 'a+') as file_write:
		for protein in ordered_protein_counts_dict :
			file_write.write("%s\t%s\t%s\n" % (protein, ordered_protein_counts_dict[protein], ordered_protein_counts_dict[protein] / sum(ordered_protein_counts_dict.values())))


	# # # # ##################################################################
	# # # # ############ SAME ANALYSIS FOR RANDOM NETWORKS ###################
	# # # # # ##################################################################

	with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\Further_Analysis\\Proteins\\full_analysis_for_proteins\\Log_degree_analysis_on_random_networks_for_proteins_%s.txt' % target, 'a+') as log_write:
		log_write.write("In covid net : len(prots_linked_to_targets) and len(set(prots_linked_to_targets)) are : %s and %s \n" % (len(prots_linked_to_targets), len(set(prots_linked_to_targets))) ) # # 9208 435
		log_write.write("In covid net : average_degree_of_my_proteins_nodes = %s\n" % average_degree_of_my_proteins_nodes)
		log_write.write("In covid net : max_degree_of_my_proteins_nodes = %s\n" % max_degree_of_my_proteins_nodes)

		nodefiles_numbers = []
		edgefiles_numbers = []
		for filename in os.listdir("A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\random_networks\\Mock_networks_21Apr"):
			if 'nodes_' in filename :
				filename = filename.split(".csv")[0]
				nodefiles_numbers.append(int(filename[19:]))
			if 'edges_' in filename :
				filename = filename.split(".csv")[0]
				edgefiles_numbers.append(int(filename[19:]))

		# log_write.write("len(nodefiles_numbers) = %s\n" % len(nodefiles_numbers)) #2051
		# log_write.write("len(edgefiles_numbers) = %s\n" % len(edgefiles_numbers)) #2051
		# log_write.write("difference between edgefiles_numbers and nodefiles_numbers = %s\n" % len(set(nodefiles_numbers).difference(set(edgefiles_numbers)))) #0

		nodefiles_numbers = sorted(nodefiles_numbers)
		edgefiles_numbers = sorted(edgefiles_numbers)

		networks_to_avoid = [3, 11, 22, 23, 33, 34, 38, 57, 72, 77, 123, 127, 148, 156, 157, 183, 205, 209, 210, 217, 232, 252, 268, 281, 291, 312, 341, 353, 362, 364, 413, 421, 433, 443, 451, 494, 525, 532, 569, 592, 629, 658, 662, 674, 676, 701, 705, 712, 736, 757, 788, 796, 825, 827, 838, 873, 876, 908, 910, 916, 917, 918, 938, 958, 965, 1088, 1125, 1128, 1145, 1146, 1165, 1171, 1183, 1185, 1192, 1195, 1211, 1244, 1254, 1257, 1300, 1353, 1354, 1364, 1372, 1394, 1395, 1401, 1410, 1427, 1454, 1477, 1497, 1517, 1532, 1549, 1551, 1575, 1578, 1611, 1624, 1629, 1637, 1650, 1655, 1697, 1700, 1713, 1736, 1745, 1762, 1790, 1800, 1815, 1824, 1827, 1866, 1872, 1886, 1887, 2053, 2072, 2106, 2115, 2128, 2160, 2217, 2343, 2383, 2388, 2413, 2458, 2465]
		nodefiles_numbers = [item for item in nodefiles_numbers if item not in networks_to_avoid]

		for iteration in nodefiles_numbers[1:] :
			print("\n########################## ITERATION %s -- %s ###########################\n" % (iteration, target))
			log_write.write("\n########################## ITERATION %s ###########################\n" % iteration)
		
			with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\random_networks\\Mock_networks_21Apr\\COVID19_GDDS_nodes_%s.csv' % str(iteration), 'r') as nodecsv: # Open the file                   
				nodereader = csv.reader(nodecsv) # Read the csv  
				nodes = [n for n in nodereader][1:]                     
				node_names = [n[0] for n in nodes] # Get a list of only the node names 

			with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\random_networks\\Mock_networks_21Apr\\COVID19_GDDS_edges_%s.csv' % str(iteration), 'r') as edgecsv: # Open the file
				edgereader = csv.reader(edgecsv) # Read the csv     
				edges = [tuple(e) for e in edgereader][1:] # Retrieve the data

				log_write.write("len(node_names) : %s\n" % len(node_names)) #16007 --> 13448
				log_write.write("len(edges) : %s\n" % len(edges)) #53917 --> 37835

			G_random = nx.Graph()
			G_random.add_nodes_from(node_names)
			G_random.add_edges_from(edges)

			descr_dict_random = {}
			description_random_set=set()
			for node in nodes: 
				descr_dict_random[node[0]] = node[2]
				description_random_set.add(node[2])

			# print(nx.info(G_random))
			log_write.write("\nRANDOM GRAPH %s\n" % iteration)
			log_write.write("%s\n" % nx.info(G_random))

			#####################################################################################################
			################################# GET EDGES CONNECTING proteinS TO target ################
			#####################################################################################################

			prots_linked_to_targets = []
			for u,v,c in G_random.edges(data=True) :	
				if (descr_dict_random[u] == 'Human PPI (target)' and descr_dict_random[v] == target) :
					prots_linked_to_targets.append(u)
			
			protein_counts = sorted(prots_linked_to_targets)

			protein_counts_dict = defaultdict( int )
			for protein in protein_counts:
			    protein_counts_dict[protein] += 1

			my_proteins_nodes = [n for n in G_random.nodes if descr_dict_random[n] == 'Human PPI (target)']
			for protein in my_proteins_nodes :
				if protein_counts_dict.get(protein) == None :
					# print("%s is missing in the linked proteins : Adding it with a frequency of 0.\n" % protein)
					protein_counts_dict[protein]= 0

			average_degree_of_my_proteins_nodes = statistics.mean([ protein_counts_dict[protein] for protein in protein_counts_dict ])
			log_write.write("Average degree of Human PPI (target) linked to %s : %s\n" % (target, average_degree_of_my_proteins_nodes))

			ordered_protein_counts_dict = OrderedDict(sorted(protein_counts_dict.items()))

			with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\Further_Analysis\\Proteins\\full_analysis_for_proteins\\prot-%s\\COVID19_GDDS_proteins_degrees_to_%s_%s.tsv' % (target, target, str(iteration)), 'a+') as file_write:
				for protein in ordered_protein_counts_dict :
					if (sum(ordered_protein_counts_dict.values()) != 0) :
						file_write.write("%s\t%s\t%s\n" % (protein, ordered_protein_counts_dict[protein], ordered_protein_counts_dict[protein] / sum(ordered_protein_counts_dict.values())))
					else :
						file_write.write("%s\t%s\t%s\n" % (protein, ordered_protein_counts_dict[protein], ordered_protein_counts_dict[protein]))



	# # # # ##################################################################
	# # # # ############ COMPUTE DEGREES MEANS AND Z-SCORES ##################
	# # # # ######## INCLUDING ABSENCE OF target AS ZERO DEGREES ###########
	# # # # ##################################################################


	# # For covid network, compute the mean and the standard deviation of proteins degrees and normalized degrees
	covid_proteins = []
	covid_protein_degree_dict = {}
	covid_protein_norm_degree_dict = {}

	total_protein_degree_dict = defaultdict(list)
	total_protein_normalized_degree_dict = defaultdict(list)

	with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\Further_Analysis\\Proteins\\full_analysis_for_proteins\\COVID19_GDDS_proteins_degrees_to_%s.tsv' % target, 'r') as data :
		for line in data :
			line = line.split('\t')
			protein = line[0]
			covid_proteins.append(protein)
			protein_degree = int(line[1])
			protein_norm_degree = float(line[2].replace('\n', ''))
			if covid_protein_degree_dict.get(protein) == None :
				covid_protein_degree_dict[protein] = protein_degree
			else :
				print("warning : %s duplicate ?!" % protein)
			covid_protein_norm_degree_dict[protein] = protein_norm_degree

			if total_protein_degree_dict.get(protein) == None :
				total_protein_degree_dict[protein] = [protein_degree]
			else :
				total_protein_degree_dict[protein].append(protein_degree)

			if total_protein_normalized_degree_dict.get(protein) == None :
				total_protein_normalized_degree_dict[protein] = [protein_norm_degree]
			else :
				total_protein_normalized_degree_dict[protein].append(protein_norm_degree)

	mean_protein_deg = statistics.mean(covid_protein_degree_dict.values())
	mean_norm_protein_deg = statistics.mean(covid_protein_norm_degree_dict.values())
	sd_protein_deg = statistics.pstdev(covid_protein_degree_dict.values())
	sd_norm_protein_deg = statistics.pstdev(covid_protein_norm_degree_dict.values())
	print(mean_protein_deg, sd_protein_deg, mean_norm_protein_deg, sd_norm_protein_deg) #
	print(len(covid_proteins),len(covid_protein_degree_dict),len(covid_protein_norm_degree_dict)) #
	# 20.14879649890591 18.9696222570052 0.002188183807439825 0.00206012405050013
	# 457 457 457






	nodefiles_numbers = []
	edgefiles_numbers = []
	for filename in os.listdir("A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\random_networks\\Mock_networks_21Apr"):
		if 'nodes_' in filename :
			filename = filename.split(".csv")[0]
			nodefiles_numbers.append(int(filename[19:]))
		if 'edges_' in filename :
			filename = filename.split(".csv")[0]
			edgefiles_numbers.append(int(filename[19:]))

	# log_write.write("len(nodefiles_numbers) = %s\n" % len(nodefiles_numbers)) #2051
	# log_write.write("len(edgefiles_numbers) = %s\n" % len(edgefiles_numbers)) #2051
	# log_write.write("difference between edgefiles_numbers and nodefiles_numbers = %s\n" % len(set(nodefiles_numbers).difference(set(edgefiles_numbers)))) #0
	print("len(nodefiles_numbers) = %s\n" % len(nodefiles_numbers)) #2051
	print("len(edgefiles_numbers) = %s\n" % len(edgefiles_numbers)) #2051
	print("difference between edgefiles_numbers and nodefiles_numbers = %s\n" % len(set(nodefiles_numbers).difference(set(edgefiles_numbers)))) #0

	# nodefiles_numbers = sorted(nodefiles_numbers)
	# edgefiles_numbers = sorted(edgefiles_numbers)

	networks_to_avoid = [3, 11, 22, 23, 33, 34, 38, 57, 72, 77, 123, 127, 148, 156, 157, 183, 205, 209, 210, 217, 232, 252, 268, 281, 291, 312, 341, 353, 362, 364, 413, 421, 433, 443, 451, 494, 525, 532, 569, 592, 629, 658, 662, 674, 676, 701, 705, 712, 736, 757, 788, 796, 825, 827, 838, 873, 876, 908, 910, 916, 917, 918, 938, 958, 965, 1088, 1125, 1128, 1145, 1146, 1165, 1171, 1183, 1185, 1192, 1195, 1211, 1244, 1254, 1257, 1300, 1353, 1354, 1364, 1372, 1394, 1395, 1401, 1410, 1427, 1454, 1477, 1497, 1517, 1532, 1549, 1551, 1575, 1578, 1611, 1624, 1629, 1637, 1650, 1655, 1697, 1700, 1713, 1736, 1745, 1762, 1790, 1800, 1815, 1824, 1827, 1866, 1872, 1886, 1887, 2053, 2072, 2106, 2115, 2128, 2160, 2217, 2343, 2383, 2388, 2413, 2458, 2465]
	nodefiles_numbers = [item for item in nodefiles_numbers if item not in networks_to_avoid]

	for iteration in nodefiles_numbers[1:] :
		# log_write.write("\n########################## ITERATION %s ###########################\n" % iteration)
		print("\n########################## ITERATION %s ###########################\n" % iteration)

		with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\Further_Analysis\\Proteins\\full_analysis_for_proteins\\prot-%s\\COVID19_GDDS_proteins_degrees_to_%s_%s.tsv' % (target, target, str(iteration)), 'r') as data: # Open the file  
			for line in data :
				line = line.split('\t')
				protein = line[0]
				degree = float(line[1])
				normalized_degree = float(line[2].replace('\n', ''))

				if total_protein_degree_dict.get(protein) == None :
					total_protein_degree_dict[protein] = [degree]
				else :
					total_protein_degree_dict[protein].append(degree)

				if total_protein_normalized_degree_dict.get(protein) == None :
					total_protein_normalized_degree_dict[protein] = [normalized_degree]
				else :
					total_protein_normalized_degree_dict[protein].append(normalized_degree)

	print(len(total_protein_degree_dict)) # 19928

	print(len(total_protein_normalized_degree_dict)) # 19928

	### make the files with all calculations
	with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\Further_Analysis\\Proteins\\full_analysis_for_proteins\\COVID19_GDDS_Results_Proteins_to_%s-withZeros.tsv' % target, 'a+') as file_write:
		file_write.write("protein\t# of nets it is found in\tmin-degree\tmax-degree\tmean-degree\tmedian-degree\tsd-degree\n")
		for protein in total_protein_degree_dict :
			nets_number = len(total_protein_degree_dict[protein])
			min_deg = min(total_protein_degree_dict[protein])
			max_deg = max(total_protein_degree_dict[protein])
			mean_deg = statistics.mean(total_protein_degree_dict[protein])
			median_deg = statistics.median(total_protein_degree_dict[protein])
			sd_deg = statistics.pstdev(total_protein_degree_dict[protein])
			file_write.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (protein, nets_number, min_deg, max_deg, mean_deg, median_deg, sd_deg))


	with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\Further_Analysis\\Proteins\\full_analysis_for_proteins\\COVID19_GDDS_Results_Proteins_to_%s_Relative-withZeros.tsv' % target, 'a+') as file_write:
		file_write.write("protein\t# of nets it is found in\tmin-degree-rel\tmax-degree-rel\tmean-degree-rel\tmedian-degree-rel\tsd-degree-rel\n")
		for protein in total_protein_normalized_degree_dict :
			nets_number = len(total_protein_normalized_degree_dict[protein])
			min_deg = min(total_protein_normalized_degree_dict[protein])
			max_deg = max(total_protein_normalized_degree_dict[protein])
			mean_deg = statistics.mean(total_protein_normalized_degree_dict[protein])
			median_deg = statistics.median(total_protein_normalized_degree_dict[protein])
			sd_deg = statistics.pstdev(total_protein_normalized_degree_dict[protein])
			file_write.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (protein, nets_number, min_deg, max_deg, mean_deg, median_deg, sd_deg))

	### compute Z-score for the covid proteins
	# For each protein found in covid network, compute 2 Z-scores.

	with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\Further_Analysis\\Proteins\\full_analysis_for_proteins\\COVID19_GDDS_Results_ForCovidProteins_to_%s-withZeros.tsv' % target, 'a+') as file_write:
		file_write.write("covid protein\t# of nets it is found in\tmin-degree\tmax-degree\tmean-degree\tmedian-degree\tsd-degree\tcovid-degree\tZ-score\n")
		for protein in covid_proteins :
			nets_number = len(total_protein_degree_dict[protein])
			min_deg = min(total_protein_degree_dict[protein])
			max_deg = max(total_protein_degree_dict[protein])
			mean_deg = statistics.mean(total_protein_degree_dict[protein])
			median_deg = statistics.median(total_protein_degree_dict[protein])
			sd_deg = statistics.pstdev(total_protein_degree_dict[protein])
			covid_deg = covid_protein_degree_dict[protein]
			if sd_deg != 0 :
				z_score = ( covid_deg - mean_deg) / sd_deg
			else : 
				z_score = 'NaN'
			file_write.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (protein, nets_number, min_deg, max_deg, mean_deg, median_deg, sd_deg, covid_deg, z_score))


	with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\Further_Analysis\\Proteins\\full_analysis_for_proteins\\COVID19_GDDS_Results_ForCovidProteins_to_%s_Relative-withZeros.tsv' % target, 'a+') as file_write:
		file_write.write("covid protein\t# of nets it is found in\tmin-degree-rel\tmax-degree-rel\tmean-degree-rel\tmedian-degree-rel\tsd-degree-rel\tcovid-degree-rel\tZ-score-rel\n")
		for protein in covid_proteins :
			nets_number = len(total_protein_normalized_degree_dict[protein])
			min_deg = min(total_protein_normalized_degree_dict[protein])
			max_deg = max(total_protein_normalized_degree_dict[protein])
			mean_deg = statistics.mean(total_protein_normalized_degree_dict[protein])
			median_deg = statistics.median(total_protein_normalized_degree_dict[protein])
			sd_deg = statistics.pstdev(total_protein_normalized_degree_dict[protein])
			covid_deg = covid_protein_norm_degree_dict[protein]
			if sd_deg != 0 :
				z_score = ( covid_deg - mean_deg) / sd_deg
			else : 
				z_score = 'NaN'
			file_write.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (protein, nets_number, min_deg, max_deg, mean_deg, median_deg, sd_deg, covid_deg, z_score))










stop = timeit.default_timer()
print(stop - start)
  
pulling_results.py,
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
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
from collections import Counter
from collections import defaultdict
from collections import OrderedDict
import statistics 
import re
import gzip
import math

# COVID19_GDDS_proteins_degrees_to_Human PPI (target)

# target_list = ['GO', 'Symptom', 'Disease','Drug', 'Human PPI (target)']

# # degrees in covid-net
# deg_dict = {}
# for target in target_list :

# 	with open('full_analysis_for_proteins\\COVID19_GDDS_proteins_degrees_to_%s.tsv' % target, 'r') as file_read :
# 		data = file_read.readlines()
# 		for line in data :
# 			line = line.split('\t')
# 			protein = line[0]
# 			degree = float(line[1])
# 			degree_rel = float(line[2].replace('\n',''))
# 			if protein not in deg_dict :
# 				deg_dict[protein] = {}
# 			deg_dict[protein][target] = degree
 
# for protein in deg_dict :
# 	print(protein, deg_dict[protein])		

# with open('pulled_results_COVID19_GDDS_proteins_degrees_to_targets.tsv', 'a+') as file_write :
# 	#heading
# 	file_write.write('protein\t')
# 	for target in target_list[:-1] :
# 		file_write.write('%s\t' % target)
# 	file_write.write('%s\n' % target_list[-1])

# 	for protein in deg_dict :
# 		file_write.write('%s\t' % protein)
# 		for target in target_list[:-1] :
# 			file_write.write('%s\t' % deg_dict[protein][target])
# 		file_write.write('%s\n' % deg_dict[protein][target_list[-1]])








# target_list = ['GO', 'Symptom', 'Disease','Drug', 'Human PPI (target)']

# # degrees in covid-net
# score_dict = {}
# for target in target_list :

# 	with open('full_analysis_for_proteins\\COVID19_GDDS_Results_ForCovidProteins_to_%s_Relative-withZeros.tsv' % target, 'r') as file_read :
# 		data = file_read.readlines()
# 		for line in data[1:] :
# 			line = line.split('\t')
# 			protein = line[0]
# 			zscore = float(line[-1].replace('\n',''))
# 			if protein not in score_dict :
# 				score_dict[protein] = {}
# 			score_dict[protein][target] = zscore
 
# for protein in score_dict :
# 	print(protein, score_dict[protein])		

# with open('pulled_results_COVID19_GDDS_proteins_relative_zscores_to_targets.tsv', 'a+') as file_write :
# 	#heading
# 	file_write.write('protein\t')
# 	for target in target_list[:-1] :
# 		file_write.write('%s\t' % target)
# 	file_write.write('%s\n' % target_list[-1])

# 	for protein in score_dict :
# 		file_write.write('%s\t' % protein)
# 		for target in target_list[:-1] :
# 			file_write.write('%s\t' % score_dict[protein][target])
# 		file_write.write('%s\n' % score_dict[protein][target_list[-1]])





target_list = ['GO', 'Symptom', 'Disease','Drug', 'Human PPI (target)']

# degrees in covid-net
score_dict = {}
for target in target_list :

	with open('full_analysis_for_proteins\\COVID19_GDDS_Results_ForCovidProteins_to_%s-withZeros.tsv' % target, 'r') as file_read :
		data = file_read.readlines()
		for line in data[1:] :
			line = line.split('\t')
			protein = line[0]
			zscore = float(line[-1].replace('\n',''))
			if math.isnan(zscore) :
				zscore = 0
			if protein not in score_dict :
				score_dict[protein] = {}
			score_dict[protein][target] = zscore
 
# for protein in score_dict :
# 	print(protein, score_dict[protein])		

with open('pulled_results_COVID19_GDDS_proteins_absolute_zscores_to_targets.tsv', 'a+') as file_write :
	#heading
	file_write.write('protein\t')
	for target in target_list[:-1] :
		file_write.write('%s\t' % target)
	file_write.write('%s\n' % target_list[-1])

	for protein in score_dict :
		file_write.write('%s\t' % protein)
		for target in target_list[:-1] :
			file_write.write('%s\t' % score_dict[protein][target])
		file_write.write('%s\n' % score_dict[protein][target_list[-1]])



data = pd.read_csv('pulled_results_COVID19_GDDS_proteins_absolute_zscores_to_targets.tsv', sep='\t')

# Correlation Matrix Heatmap
# f, ax = plt.subplots(figsize=(10, 6))
# corr = data.corr(method='kendall')
# hm = sns.heatmap(round(corr,2), annot=True, ax=ax, cmap="coolwarm",fmt='.2f',
#                  linewidths=.05)
# f.subplots_adjust(top=0.93)
# t= f.suptitle('Protein Absolute Zscores of Degrees Correlation Heatmap', fontsize=14)
# plt.show()


g = sns.pairplot(data, corner=True, diag_kind="kde", kind="reg", dropna=True)
# g.fig.suptitle('Protein Abs Degree Zscores Pairplot', x=0.2, y=1, ha='left')
plt.show()




# # Pair-wise Scatter Plots
# cols = ['GO', 'Symptom', 'Disease','Drug', 'Human PPI (target)']
# pp = sns.pairplot(data[cols], height=1.8, aspect=1.8,
#                   plot_kws=dict(edgecolor="k", linewidth=0.5),
#                   diag_kind="kde", diag_kws=dict(shade=True))

# fig = pp.fig 
# fig.subplots_adjust(top=0.93, wspace=0.3)
# t = fig.suptitle('Protein Relative Degree Zscores Pairwise Plots', fontsize=14)
# fig.show()


# sns.set(style="ticks", color_codes=True)
# iris = sns.load_dataset("iris")
# print(iris)
# g = sns.pairplot(iris)
# g = sns.pairplot(iris, hue="species")
# g = sns.pairplot(data, corner=True, diag_kind="kde", kind="reg")
# g.fig.suptitle('Protein Abs Degrees Pairplot', x=0.2, y=1, ha='left')
# plt.show()

# sns.set(style="ticks", color_codes=True)
# # g = sns.pairplot(iris)
# # g = sns.pairplot(iris, hue="species")
# # g = sns.pairplot(iris, corner=True, diag_kind="kde", kind="reg")
# g = sns.PairGrid(data)
# g = g.map(plt.scatter)
# plt.show()














# targets = [1,2,3,4,5]
# b = {}
# data = ['prot1','prot2','prot3']
# for target in targets :
# 	for prot in data :
# 		deg = 'deg_%s_%s' % (prot,target)
# 		degrel = 'degrel_%s_%s' % (prot,target)
# 		b[prot] = {target : (deg, degrel)}

# for i in b :
# 	print(i, b[i])







# with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\COVID19_GDDS_nodes.csv', 'r') as nodecsv: # Open the file                       
# 	nodereader = csv.reader(nodecsv) # Read the csv  
# 	nodes = [n for n in nodereader][1:]                     
# 	node_names = [n[0] for n in nodes] # Get a list of only the node names 

# with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\COVID19_GDDS_edges.csv', 'r') as edgecsv: # Open the file
# 	edgereader = csv.reader(edgecsv) # Read the csv     
# 	edges = [tuple(e) for e in edgereader][1:] # Retrieve the data

# G = nx.Graph()
# G.add_edges_from(edges)


# descr_dict = {}
# for node in nodes: 
# 	descr_dict[node[0]] = node[2]

# prots_linked_to_targets = []
# for u,v,c in G.edges(data=True) :	
# 	if (descr_dict[u] == 'Human PPI (target)' and descr_dict[v] == 'Drug') :
# 		prots_linked_to_targets.append(u)
# if 'ACPP' in prots_linked_to_targets :
# 	print('found ACPP')






stop = timeit.default_timer()
print(stop - start)
  
pulling_results_undirected.py,
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
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
from collections import Counter
from collections import defaultdict
from collections import OrderedDict
import statistics 
import re
import gzip
import math

# COVID19_GDDS_proteins_degrees_to_Human PPI (target)

# target_list = ['GO', 'Symptom', 'Disease','Drug', 'Human PPI (target)']

# # degrees in covid-net
# deg_dict = {}
# for target in target_list :

# 	with open('full_analysis_for_proteins_undirected\\COVID19_GDDS_proteins_degrees_to_%s.tsv' % target, 'r') as file_read :
# 		data = file_read.readlines()
# 		for line in data :
# 			line = line.split('\t')
# 			protein = line[0]
# 			degree = float(line[1])
# 			degree_rel = float(line[2].replace('\n',''))
# 			if protein not in deg_dict :
# 				deg_dict[protein] = {}
# 			deg_dict[protein][target] = degree
 
# for protein in deg_dict :
# 	print(protein, deg_dict[protein])		

# with open('pulled_results_COVID19_GDDS_proteins_degrees_to_targets_undirected.tsv', 'a+') as file_write :
# 	#heading
# 	file_write.write('protein\t')
# 	for target in target_list[:-1] :
# 		file_write.write('%s\t' % target)
# 	file_write.write('%s\n' % target_list[-1])

# 	for protein in deg_dict :
# 		file_write.write('%s\t' % protein)
# 		for target in target_list[:-1] :
# 			file_write.write('%s\t' % deg_dict[protein][target])
# 		file_write.write('%s\n' % deg_dict[protein][target_list[-1]])









# target_list = ['GO', 'Symptom', 'Disease','Drug', 'Human PPI (target)']

# # degrees in covid-net
# deg_dict = {}
# for target in target_list :

# 	with open('full_analysis_for_proteins_undirected\\COVID19_GDDS_proteins_degrees_to_%s.tsv' % target, 'r') as file_read :
# 		data = file_read.readlines()
# 		for line in data :
# 			line = line.split('\t')
# 			protein = line[0]
# 			degree = float(line[1])
# 			degree_rel = float(line[2].replace('\n',''))
# 			if protein not in deg_dict :
# 				deg_dict[protein] = {}
# 			deg_dict[protein][target] = degree_rel
 
# for protein in deg_dict :
# 	print(protein, deg_dict[protein])		

# with open('pulled_results_COVID19_GDDS_proteins_relative_degrees_to_targets_undirected.tsv', 'a+') as file_write :
# 	#heading
# 	file_write.write('protein\t')
# 	for target in target_list[:-1] :
# 		file_write.write('%s\t' % target)
# 	file_write.write('%s\n' % target_list[-1])

# 	for protein in deg_dict :
# 		file_write.write('%s\t' % protein)
# 		for target in target_list[:-1] :
# 			file_write.write('%s\t' % deg_dict[protein][target])
# 		file_write.write('%s\n' % deg_dict[protein][target_list[-1]])









# target_list = ['GO', 'Symptom', 'Disease','Drug', 'Human PPI (target)']

# # degrees in covid-net
# score_dict = {}
# for target in target_list :

# 	with open('full_analysis_for_proteins_undirected\\COVID19_GDDS_Results_ForCovidProteins_to_%s_Relative-withZeros.tsv' % target, 'r') as file_read :
# 		data = file_read.readlines()
# 		for line in data[1:] :
# 			line = line.split('\t')
# 			protein = line[0]
# 			zscore = float(line[-1].replace('\n',''))
# 			if protein not in score_dict :
# 				score_dict[protein] = {}
# 			score_dict[protein][target] = zscore
 
# for protein in score_dict :
# 	print(protein, score_dict[protein])		

# with open('pulled_results_COVID19_GDDS_proteins_relative_zscores_to_targets_undirected.tsv', 'a+') as file_write :
# 	#heading
# 	file_write.write('protein\t')
# 	for target in target_list[:-1] :
# 		file_write.write('%s\t' % target)
# 	file_write.write('%s\n' % target_list[-1])

# 	for protein in score_dict :
# 		file_write.write('%s\t' % protein)
# 		for target in target_list[:-1] :
# 			file_write.write('%s\t' % score_dict[protein][target])
# 		file_write.write('%s\n' % score_dict[protein][target_list[-1]])









from scipy.stats import spearmanr
def corrfunc(x,y, ax=None, **kws):
    """Plot the correlation coefficient in the top left hand corner of a plot."""
    r, _ = spearmanr(x, y)
    ax = ax or plt.gca()
    # Unicode for lowercase rho (ρ)
    rho = '\u03C1'
    ax.annotate(f'{rho}_spearman = {r:.2f}', xy=(.1, .9), xycoords=ax.transAxes)

# data = pd.read_csv('pulled_results_COVID19_GDDS_proteins_relative_degrees_to_targets_undirected.tsv', usecols=['Drug', 'Human PPI (target)'], sep='\t')
data = pd.read_csv('pulled_results_COVID19_GDDS_proteins_relative_degrees_to_targets_undirected.tsv', sep='\t')

# Correlation Matrix Heatmap
# f, ax = plt.subplots(figsize=(10, 6))
# corr = data.corr(method='spearman')
# hm = sns.heatmap(round(corr,2), annot=True, ax=ax, cmap="coolwarm",fmt='.2f',
#                  linewidths=.05)
# f.subplots_adjust(top=0.93)
# t= f.suptitle('Protein Abs Degrees Undirected Correlation Heatmap', fontsize=14)
# plt.show()


g = sns.pairplot(data, corner=True, diag_kind="kde", kind="reg", dropna=True)
# g.map_lower(corrfunc)
g.fig.suptitle('Protein Relative Degrees Undirected Pairplot', x=0.2, y=1, ha='left')
plt.show()




# # Pair-wise Scatter Plots
# cols = ['GO', 'Symptom', 'Disease','Drug', 'Human PPI (target)']
# pp = sns.pairplot(data[cols], height=1.8, aspect=1.8,
#                   plot_kws=dict(edgecolor="k", linewidth=0.5),
#                   diag_kind="kde", diag_kws=dict(shade=True))

# fig = pp.fig 
# fig.subplots_adjust(top=0.93, wspace=0.3)
# t = fig.suptitle('Protein Relative Degree Zscores Pairwise Plots', fontsize=14)
# fig.show()


# sns.set(style="ticks", color_codes=True)
# iris = sns.load_dataset("iris")
# print(iris)
# g = sns.pairplot(iris)
# g = sns.pairplot(iris, hue="species")
# g = sns.pairplot(data, corner=True, diag_kind="kde", kind="reg")
# g.fig.suptitle('Protein Abs Degrees Pairplot', x=0.2, y=1, ha='left')
# plt.show()

# sns.set(style="ticks", color_codes=True)
# # g = sns.pairplot(iris)
# # g = sns.pairplot(iris, hue="species")
# # g = sns.pairplot(iris, corner=True, diag_kind="kde", kind="reg")
# g = sns.PairGrid(data)
# g = g.map(plt.scatter)
# plt.show()














# targets = [1,2,3,4,5]
# b = {}
# data = ['prot1','prot2','prot3']
# for target in targets :
# 	for prot in data :
# 		deg = 'deg_%s_%s' % (prot,target)
# 		degrel = 'degrel_%s_%s' % (prot,target)
# 		b[prot] = {target : (deg, degrel)}

# for i in b :
# 	print(i, b[i])







# with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\COVID19_GDDS_nodes.csv', 'r') as nodecsv: # Open the file                       
# 	nodereader = csv.reader(nodecsv) # Read the csv  
# 	nodes = [n for n in nodereader][1:]                     
# 	node_names = [n[0] for n in nodes] # Get a list of only the node names 

# with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\COVID19_GDDS_edges.csv', 'r') as edgecsv: # Open the file
# 	edgereader = csv.reader(edgecsv) # Read the csv     
# 	edges = [tuple(e) for e in edgereader][1:] # Retrieve the data

# G = nx.Graph()
# G.add_edges_from(edges)


# descr_dict = {}
# for node in nodes: 
# 	descr_dict[node[0]] = node[2]

# prots_linked_to_targets = []
# for u,v,c in G.edges(data=True) :	
# 	if (descr_dict[u] == 'Human PPI (target)' and descr_dict[v] == 'Drug') :
# 		prots_linked_to_targets.append(u)
# if 'ACPP' in prots_linked_to_targets :
# 	print('found ACPP')






stop = timeit.default_timer()
print(stop - start)
  
random-network-analysis-disease-ranking-filtering.py,
#!/usr/bin/python
# -*- coding: utf-8 -*-
# jupyter notebook --browser='C:/Program Files (x86)/Google/Chrome/Application/chrome.exe %s'


import timeit
start = timeit.default_timer()

import os
import csv
from operator import itemgetter
import networkx as nx
from networkx.algorithms import community
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
from collections import Counter
from collections import defaultdict
from collections import OrderedDict
import statistics 
import re
import gzip





#####################################################################################################
#########################  Compute network Diseases degree distribution #################################
##########################################  based on all edges  #####################################
#####################################################################################################

with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\COVID19_GDDS_nodes.csv', 'r') as nodecsv: # Open the file                       
	nodereader = csv.reader(nodecsv) # Read the csv  
	nodes = [n for n in nodereader][1:]                     
	node_names = [n[0] for n in nodes] # Get a list of only the node names 

with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\COVID19_GDDS_edges.csv', 'r') as edgecsv: # Open the file
	edgereader = csv.reader(edgecsv) # Read the csv     
	edges = [tuple(e) for e in edgereader][1:] # Retrieve the data

G = nx.Graph()
G.add_nodes_from(node_names)
G.add_edges_from(edges)

descr_dict = {}
for node in nodes: 
	descr_dict[node[0]] = node[2]

print(nx.info(G))
# Number of nodes: 16007
# Number of edges: 53917
# Average degree:   6.7367

# print(nx.degree(G)) # this is the list containing as many tuples as nodes, indicating each node's degree.

average_degree = statistics.mean([ tpl[1] for tpl in nx.degree(G) ])
# print(average_degree) # 6.736677703504717

diseases_degrees = nx.degree(G, nbunch=[n for n in G.nodes if descr_dict[n] == 'Disease'])
# print("degree of disease nodes only", diseases_degrees)

average_degree_of_diseases = statistics.mean([ tpl[1] for tpl in diseases_degrees ])
print("average degree of diseases nodes", average_degree_of_diseases) # 4.329980842911877

# #####################################################################################################
# ################################# CALCULATION ON EDGES WITH FILTERING  ###########################
# #####################################################################################################
# # let's count edges only in order to calculate the nodes' degrees based on the number of edges they make in the graph that connect them to proteins that are themselves connected to viral genes

#first we have to build the set of proteins of interest
proteins_of_interest = set()
for u,v,c in G.edges(data=True) :	
	if (descr_dict[u] == 'Viral Gene' and descr_dict[v] == 'Human PPI (target)') :
		proteins_of_interest.add(v)
print(len(proteins_of_interest)) # 332 OK

# now we keep only the disease that interact with those proteins of interest
diseases_linked_to_some_proteins_of_interest = []
for u,v,c in G.edges(data=True) :	
	if (descr_dict[u] == 'Disease' and v in proteins_of_interest) :
		diseases_linked_to_some_proteins_of_interest.append(u) # this never happens
	if (descr_dict[v] == 'Disease' and u in proteins_of_interest) :
		diseases_linked_to_some_proteins_of_interest.append(v)
print(len(diseases_linked_to_some_proteins_of_interest), len(set(diseases_linked_to_some_proteins_of_interest))) # 9430 3164
# ({'Drug': 5703, 'Disease': 4176, 'GO': 3487, 'Symptom': 2157, 'Human PPI (target)': 457, 'Viral Gene': 27})

# # print(sorted(diseases_linked_to_some_proteins_of_interest))

disease_counts = sorted(diseases_linked_to_some_proteins_of_interest)

# # # for the diseases
disease_counts_dict = defaultdict( int )
for disease in disease_counts:
    disease_counts_dict[disease] += 1

my_diseases_nodes = [n for n in G.nodes if descr_dict[n] == 'Disease']
for disease in my_diseases_nodes :
	if disease_counts_dict.get(disease) == None :
		# print("%s is missing in the linked diseases : Adding it with a frequency of 0.\n" % disease)
		disease_counts_dict[disease]= 0

average_degree_of_my_diseases_nodes = statistics.mean([ disease_counts_dict[disease] for disease in disease_counts_dict ])
print("average_degree_of_my_diseases_nodes", average_degree_of_my_diseases_nodes)
# 2.2581417624521074

max_degree_of_my_diseases_nodes = max([ disease_counts_dict[disease] for disease in disease_counts_dict ])
print("max_degree_of_my_diseases_nodes", max_degree_of_my_diseases_nodes) #93


ordered_disease_counts_dict = OrderedDict(sorted(disease_counts_dict.items()))

# with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\Further_Analysis\\Diseases\\First_order_proteins\\COVID19_GDDS_diseases_degrees_to_first_order_proteins.tsv', 'a+') as file_write:
# 	for disease in ordered_disease_counts_dict :
# 		file_write.write("%s\t%s\t%s\n" % (disease, ordered_disease_counts_dict[disease], ordered_disease_counts_dict[disease] / sum(ordered_disease_counts_dict.values())))


# # ##################################################################
# # ############ SAME ANALYSIS FOR RANDOM NETWORKS ###################
# # # ##################################################################

# with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\Further_Analysis\\Diseases\\First_order_proteins\\Log_degree_analysis_on_random_networks_for_diseases_proteins_of_interest.txt', 'a+') as log_write:

# 	nodefiles_numbers = []
# 	edgefiles_numbers = []
# 	for filename in os.listdir("A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\random_networks\\Mock_networks_21Apr"):
# 		if 'nodes_' in filename :
# 			filename = filename.split(".csv")[0]
# 			nodefiles_numbers.append(int(filename[19:]))
# 		if 'edges_' in filename :
# 			filename = filename.split(".csv")[0]
# 			edgefiles_numbers.append(int(filename[19:]))

# 	log_write.write("len(nodefiles_numbers) = %s\n" % len(nodefiles_numbers)) #2051
# 	log_write.write("len(edgefiles_numbers) = %s\n" % len(edgefiles_numbers)) #2051
# 	log_write.write("difference between edgefiles_numbers and nodefiles_numbers = %s\n" % len(set(nodefiles_numbers).difference(set(edgefiles_numbers)))) #0

# 	nodefiles_numbers = sorted(nodefiles_numbers)
# 	edgefiles_numbers = sorted(edgefiles_numbers)

# 	for iteration in nodefiles_numbers[1:] :
# 		print("\n########################## ITERATION %s ###########################\n" % iteration)
# 		log_write.write("\n########################## ITERATION %s ###########################\n" % iteration)
	
# 		with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\random_networks\\Mock_networks_21Apr\\COVID19_GDDS_nodes_%s.csv' % str(iteration), 'r') as nodecsv: # Open the file                   
# 			nodereader = csv.reader(nodecsv) # Read the csv  
# 			nodes = [n for n in nodereader][1:]                     
# 			node_names = [n[0] for n in nodes] # Get a list of only the node names 

# 		with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\random_networks\\Mock_networks_21Apr\\COVID19_GDDS_edges_%s.csv' % str(iteration), 'r') as edgecsv: # Open the file
# 			edgereader = csv.reader(edgecsv) # Read the csv     
# 			edges = [tuple(e) for e in edgereader][1:] # Retrieve the data

# 			log_write.write("len(node_names) : %s\n" % len(node_names)) #16007 --> 13448
# 			log_write.write("len(edges) : %s\n" % len(edges)) #53917 --> 37835

# 		G = nx.Graph()
# 		G.add_nodes_from(node_names)
# 		G.add_edges_from(edges)

# 		descr_dict = {}
# 		description_set=set()
# 		for node in nodes: 
# 			descr_dict[node[0]] = node[2]
# 			description_set.add(node[2])

# 		# print(description_set) #{'Symptom', 'Human PPI (target)', 'GO', 'Disease', 'Drug'}
# 		# recounted = Counter(list(descr_dict.values()))
# 		# print("%s\n" % recounted)

# 		# print(nx.info(G))
# 		log_write.write("\nRANDOM GRAPH %s\n" % iteration)
# 		log_write.write("%s\n" % nx.info(G))

# 		#####################################################################################################
# 		################################# GET EDGES CONNECTING DISEASES TO FIRST ORDER PROTEINS ################
# 		#####################################################################################################
# 		#first we have to build the set of proteins of interest
# 		proteins_of_interest = set([n for n in G.nodes if descr_dict[n] == 'Human PPI (target)'])
# 		# print(proteins_of_interest)

# 		# now we keep only the diseases that interact with those proteins of interest
# 		diseases_linked_to_some_proteins_of_interest = []
# 		for u,v,c in G.edges(data=True) :	
# 			if (descr_dict[u] == 'Disease' and v in proteins_of_interest) :
# 				diseases_linked_to_some_proteins_of_interest.append(u) # this never happens
# 			if (descr_dict[v] == 'Disease' and u in proteins_of_interest) :
# 				diseases_linked_to_some_proteins_of_interest.append(v)
# 		# print(len(diseases_linked_to_some_proteins_of_interest), len(set(diseases_linked_to_some_proteins_of_interest))) # 8200 3740
		
# 		disease_counts = sorted(diseases_linked_to_some_proteins_of_interest)

# 		# for the diseases
# 		disease_counts_dict = defaultdict( int )
# 		for disease in disease_counts:
# 		    disease_counts_dict[disease] += 1

# 		my_diseases_nodes = [n for n in G.nodes if descr_dict[n] == 'Disease']
# 		for disease in my_diseases_nodes :
# 			if disease_counts_dict.get(disease) == None :
# 				# print("%s is missing in the linked disease : Adding it with a frequency of 0.\n" % disease)
# 				disease_counts_dict[disease]= 0

# 		average_degree_of_my_diseases_nodes = statistics.mean([ disease_counts_dict[disease] for disease in disease_counts_dict ])
# 		# print("average_degree_of_my_diseases_nodes", average_degree_of_my_diseases_nodes)
# 		log_write.write("Average degree of diseases to a protein : %s\n" % average_degree_of_my_diseases_nodes)

# 		ordered_disease_counts_dict = OrderedDict(sorted(disease_counts_dict.items()))

# 		with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\Further_Analysis\\Random_Networks\\Diseases\\First_order_proteins\\COVID19_GDDS_diseases_degrees_to_first_order_proteins_%s.tsv' % str(iteration), 'a+') as file_write:
# 			for disease in ordered_disease_counts_dict :
# 				if (sum(ordered_disease_counts_dict.values()) != 0) :
# 					file_write.write("%s\t%s\t%s\n" % (disease, ordered_disease_counts_dict[disease], ordered_disease_counts_dict[disease] / sum(ordered_disease_counts_dict.values())))
# 				else :
# 					file_write.write("%s\t%s\t%s\n" % (disease, ordered_disease_counts_dict[disease], ordered_disease_counts_dict[disease]))



# # ##################################################################
# # ############ COMPUTE DEGREES MEANS AND Z-SCORES ##################
# # ######## INCLUDING ABSENCE OF DISEASES AS ZERO DEGREES ###########
# # ##################################################################


# # For covid network, compute the mean and the standard deviation of diseases degrees and normalized degrees
covid_diseases = []
covid_disease_degree_dict = {}
covid_disease_norm_degree_dict = {}

total_disease_degree_dict = defaultdict(list)
total_disease_normalized_degree_dict = defaultdict(list)

with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\Further_Analysis\\Diseases\\First_order_proteins\\COVID19_GDDS_diseases_degrees_to_first_order_proteins.tsv', 'r') as data :
	for line in data :
		line = line.split('\t')
		disease = line[0]
		covid_diseases.append(disease)
		disease_degree = int(line[1])
		disease_norm_degree = float(line[2].replace('\n', ''))
		if covid_disease_degree_dict.get(disease) == None :
			covid_disease_degree_dict[disease] = disease_degree
		else :
			print("warning : %s duplicate ?!" % disease)
		covid_disease_norm_degree_dict[disease] = disease_norm_degree

		if total_disease_degree_dict.get(disease) == None :
			total_disease_degree_dict[disease] = [disease_degree]
		else :
			total_disease_degree_dict[disease].append(disease_degree)

		if total_disease_normalized_degree_dict.get(disease) == None :
			total_disease_normalized_degree_dict[disease] = [disease_norm_degree]
		else :
			total_disease_normalized_degree_dict[disease].append(disease_norm_degree)

mean_disease_deg = statistics.mean(covid_disease_degree_dict.values())
mean_norm_disease_deg = statistics.mean(covid_disease_norm_degree_dict.values())
sd_disease_deg = statistics.pstdev(covid_disease_degree_dict.values())
sd_norm_disease_deg = statistics.pstdev(covid_disease_norm_degree_dict.values())
print(mean_disease_deg, sd_disease_deg, mean_norm_disease_deg, sd_norm_disease_deg) #1.4378397334736104 2.251958221466831 0.0001753463089601964 0.00027462905139839403

print(len(covid_diseases),len(covid_disease_degree_dict),len(covid_disease_norm_degree_dict)) #5703 5703 5703


nodefiles_numbers = []
edgefiles_numbers = []
for filename in os.listdir("A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\random_networks\\Mock_networks_21Apr"):
	if 'nodes_' in filename :
		filename = filename.split(".csv")[0]
		nodefiles_numbers.append(int(filename[19:]))
	if 'edges_' in filename :
		filename = filename.split(".csv")[0]
		edgefiles_numbers.append(int(filename[19:]))

# log_write.write("len(nodefiles_numbers) = %s\n" % len(nodefiles_numbers)) #2051
# log_write.write("len(edgefiles_numbers) = %s\n" % len(edgefiles_numbers)) #2051
# log_write.write("difference between edgefiles_numbers and nodefiles_numbers = %s\n" % len(set(nodefiles_numbers).difference(set(edgefiles_numbers)))) #0
print("len(nodefiles_numbers) = %s\n" % len(nodefiles_numbers)) #2051
print("len(edgefiles_numbers) = %s\n" % len(edgefiles_numbers)) #2051
print("difference between edgefiles_numbers and nodefiles_numbers = %s\n" % len(set(nodefiles_numbers).difference(set(edgefiles_numbers)))) #0

nodefiles_numbers = sorted(nodefiles_numbers)
edgefiles_numbers = sorted(edgefiles_numbers)


for iteration in nodefiles_numbers[1:] :
	# log_write.write("\n########################## ITERATION %s ###########################\n" % iteration)
	print("\n########################## ITERATION %s ###########################\n" % iteration)

	with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\Further_Analysis\\Random_Networks\\Diseases\\First_order_proteins\\COVID19_GDDS_diseases_degrees_to_first_order_proteins_%s.tsv' % str(iteration), 'r') as data: # Open the file  
		for line in data :
			line = line.split('\t')
			disease = line[0]
			degree = float(line[1])
			normalized_degree = float(line[2].replace('\n', ''))

			if total_disease_degree_dict.get(disease) == None :
				total_disease_degree_dict[disease] = [degree]
			else :
				total_disease_degree_dict[disease].append(degree)

			if total_disease_normalized_degree_dict.get(disease) == None :
				total_disease_normalized_degree_dict[disease] = [normalized_degree]
			else :
				total_disease_normalized_degree_dict[disease].append(normalized_degree)

print(len(total_disease_degree_dict)) # 72069
print(len(total_disease_normalized_degree_dict)) # 72069

### make the files with all calculations
with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\Further_Analysis\\Diseases\\First_order_proteins\\COVID19_GDDS_Results_Diseases_first_order_proteins-withZeros.tsv', 'a+') as file_write:
	file_write.write("disease\t# of nets it is found in\tmin-degree\tmax-degree\tmean-degree\tmedian-degree\tsd-degree\n")
	for disease in total_disease_degree_dict :
		nets_number = len(total_disease_degree_dict[disease])
		min_deg = min(total_disease_degree_dict[disease])
		max_deg = max(total_disease_degree_dict[disease])
		mean_deg = statistics.mean(total_disease_degree_dict[disease])
		median_deg = statistics.median(total_disease_degree_dict[disease])
		sd_deg = statistics.pstdev(total_disease_degree_dict[disease])
		file_write.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (disease, nets_number, min_deg, max_deg, mean_deg, median_deg, sd_deg))


with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\Further_Analysis\\Diseases\\\First_order_proteins\\COVID19_GDDS_Results_Diseases_first_order_proteins_Relative-withZeros.tsv', 'a+') as file_write:
	file_write.write("disease\t# of nets it is found in\tmin-degree-rel\tmax-degree-rel\tmean-degree-rel\tmedian-degree-rel\tsd-degree-rel\n")
	for disease in total_disease_normalized_degree_dict :
		nets_number = len(total_disease_normalized_degree_dict[disease])
		min_deg = min(total_disease_normalized_degree_dict[disease])
		max_deg = max(total_disease_normalized_degree_dict[disease])
		mean_deg = statistics.mean(total_disease_normalized_degree_dict[disease])
		median_deg = statistics.median(total_disease_normalized_degree_dict[disease])
		sd_deg = statistics.pstdev(total_disease_normalized_degree_dict[disease])
		file_write.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (disease, nets_number, min_deg, max_deg, mean_deg, median_deg, sd_deg))

### compute Z-score for the covid diseases
# For each disease found in covid network, compute 2 Z-scores.

with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\Further_Analysis\\Diseases\\First_order_proteins\\COVID19_GDDS_Results_ForCovidDiseases-first_order_proteins-withZeros.tsv', 'a+') as file_write:
	file_write.write("covid disease\t# of nets it is found in\tmin-degree\tmax-degree\tmean-degree\tmedian-degree\tsd-degree\tcovid-degree\tZ-score\n")
	for disease in covid_diseases :
		nets_number = len(total_disease_degree_dict[disease])
		min_deg = min(total_disease_degree_dict[disease])
		max_deg = max(total_disease_degree_dict[disease])
		mean_deg = statistics.mean(total_disease_degree_dict[disease])
		median_deg = statistics.median(total_disease_degree_dict[disease])
		sd_deg = statistics.pstdev(total_disease_degree_dict[disease])
		covid_deg = covid_disease_degree_dict[disease]
		if sd_deg != 0 :
			z_score = ( covid_deg - mean_deg) / sd_deg
		else : 
			z_score = 'NaN'
		file_write.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (disease, nets_number, min_deg, max_deg, mean_deg, median_deg, sd_deg, covid_deg, z_score))


with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\Further_Analysis\\Diseases\\First_order_proteins\\COVID19_GDDS_Results_ForCovidDiseases_first_order_proteins_Relative-withZeros.tsv', 'a+') as file_write:
	file_write.write("covid disease\t# of nets it is found in\tmin-degree-rel\tmax-degree-rel\tmean-degree-rel\tmedian-degree-rel\tsd-degree-rel\tcovid-degree-rel\tZ-score-rel\n")
	for disease in covid_diseases :
		nets_number = len(total_disease_normalized_degree_dict[disease])
		min_deg = min(total_disease_normalized_degree_dict[disease])
		max_deg = max(total_disease_normalized_degree_dict[disease])
		mean_deg = statistics.mean(total_disease_normalized_degree_dict[disease])
		median_deg = statistics.median(total_disease_normalized_degree_dict[disease])
		sd_deg = statistics.pstdev(total_disease_normalized_degree_dict[disease])
		covid_deg = covid_disease_norm_degree_dict[disease]
		if sd_deg != 0 :
			z_score = ( covid_deg - mean_deg) / sd_deg
		else : 
			z_score = 'NaN'
		file_write.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (disease, nets_number, min_deg, max_deg, mean_deg, median_deg, sd_deg, covid_deg, z_score))










stop = timeit.default_timer()
print(stop - start)
 
random-network-analysis-drugs-ranking-filtering.py,
#!/usr/bin/python
# -*- coding: utf-8 -*-
# jupyter notebook --browser='C:/Program Files (x86)/Google/Chrome/Application/chrome.exe %s'


import timeit
start = timeit.default_timer()

import os
import csv
from operator import itemgetter
import networkx as nx
from networkx.algorithms import community
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
from collections import Counter
from collections import defaultdict
from collections import OrderedDict
import statistics 
import re
import gzip





#####################################################################################################
#########################  Compute network Drug degree distribution #################################
##########################################  based on all edges  #####################################
#####################################################################################################

with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\COVID19_GDDS_nodes.csv', 'r') as nodecsv: # Open the file                       
	nodereader = csv.reader(nodecsv) # Read the csv  
	nodes = [n for n in nodereader][1:]                     
	node_names = [n[0] for n in nodes] # Get a list of only the node names 

with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\COVID19_GDDS_edges.csv', 'r') as edgecsv: # Open the file
	edgereader = csv.reader(edgecsv) # Read the csv     
	edges = [tuple(e) for e in edgereader][1:] # Retrieve the data

G = nx.Graph()
G.add_nodes_from(node_names)
G.add_edges_from(edges)

descr_dict = {}
for node in nodes: 
	descr_dict[node[0]] = node[2]

print(nx.info(G))
# Number of nodes: 16007
# Number of edges: 53917
# Average degree:   6.7367

# print(nx.degree(G)) # this is the list containing as many tuples as nodes, indicating each node's degree.

average_degree = statistics.mean([ tpl[1] for tpl in nx.degree(G) ])
# print(average_degree) # 6.736677703504717

drugs_degrees = nx.degree(G, nbunch=[n for n in G.nodes if descr_dict[n] == 'Drug'])
# print("degree of drugs nodes only", drugs_degrees)

average_degree_of_drugs = statistics.mean([ tpl[1] for tpl in drugs_degrees ])
# print("average degree of drugs nodes", average_degree_of_drugs) # 2.2826582500438364

# #####################################################################################################
# ################################# CALCULATION ON EDGES WITH FILTERING  ###########################
# #####################################################################################################
# # let's count edges only in order to calculate the nodes' degrees based on the number of edges they make in the graph that connect them to proteins that are themselves connected to viral genes

#first we have to build the set of proteins of interest
proteins_of_interest = set()
for u,v,c in G.edges(data=True) :	
	if (descr_dict[u] == 'Viral Gene' and descr_dict[v] == 'Human PPI (target)') :
		proteins_of_interest.add(v)
print(len(proteins_of_interest)) # 332 OK

# now we keep only the drug that interact with those proteins of interest
drugs_linked_to_some_proteins_of_interest = []
for u,v,c in G.edges(data=True) :	
	if (descr_dict[u] == 'Drug' and v in proteins_of_interest) :
		drugs_linked_to_some_proteins_of_interest.append(u) # this never happens
	if (descr_dict[v] == 'Drug' and u in proteins_of_interest) :
		drugs_linked_to_some_proteins_of_interest.append(v)
print(len(drugs_linked_to_some_proteins_of_interest), len(set(drugs_linked_to_some_proteins_of_interest))) # 8200 3740
# ({'Drug': 5703, 'Disease': 4176, 'GO': 3487, 'Symptom': 2157, 'Human PPI (target)': 457, 'Viral Gene': 27})

# # print(sorted(drugs_linked_to_some_proteins_of_interest))

drug_counts = sorted(drugs_linked_to_some_proteins_of_interest)

# # # for the drugs
drug_counts_dict = defaultdict( int )
for drug in drug_counts:
    drug_counts_dict[drug] += 1

my_drugs_nodes = [n for n in G.nodes if descr_dict[n] == 'Drug']
for drug in my_drugs_nodes :
	if drug_counts_dict.get(drug) == None :
		# print("%s is missing in the linked drugs : Adding it with a frequency of 0.\n" % drug)
		drug_counts_dict[drug]= 0

average_degree_of_my_drugs_nodes = statistics.mean([ drug_counts_dict[drug] for drug in drug_counts_dict ])
print("average_degree_of_my_drugs_nodes", average_degree_of_my_drugs_nodes)
# 1.4378397334736104

max_degree_of_my_drugs_nodes = max([ drug_counts_dict[drug] for drug in drug_counts_dict ])
print("max_degree_of_my_drugs_nodes", max_degree_of_my_drugs_nodes) #14


ordered_drug_counts_dict = OrderedDict(sorted(drug_counts_dict.items()))

# with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\Further_Analysis\\Drugs\\COVID19_GDDS_drugs_degrees_to_proteins_of_interest.tsv', 'a+') as file_write:
# 	for drug in ordered_drug_counts_dict :
# 		file_write.write("%s\t%s\t%s\n" % (drug, ordered_drug_counts_dict[drug], ordered_drug_counts_dict[drug] / sum(ordered_drug_counts_dict.values())))


# # ##################################################################
# # ############ SAME ANALYSIS FOR RANDOM NETWORKS ###################
# # ##################################################################

# with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\Further_Analysis\\Drugs\\Log_degree_analysis_on_random_networks_for_drugs_proteins_of_interest.txt', 'a+') as log_write:

# 	nodefiles_numbers = []
# 	edgefiles_numbers = []
# 	for filename in os.listdir("A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\random_networks\\Mock_networks_21Apr"):
# 		if 'nodes_' in filename :
# 			filename = filename.split(".csv")[0]
# 			nodefiles_numbers.append(int(filename[19:]))
# 		if 'edges_' in filename :
# 			filename = filename.split(".csv")[0]
# 			edgefiles_numbers.append(int(filename[19:]))

# 	log_write.write("len(nodefiles_numbers) = %s\n" % len(nodefiles_numbers)) #2051
# 	log_write.write("len(edgefiles_numbers) = %s\n" % len(edgefiles_numbers)) #2051
# 	log_write.write("difference between edgefiles_numbers and nodefiles_numbers = %s\n" % len(set(nodefiles_numbers).difference(set(edgefiles_numbers)))) #0

# 	nodefiles_numbers = sorted(nodefiles_numbers)
# 	edgefiles_numbers = sorted(edgefiles_numbers)

# 	for iteration in nodefiles_numbers[1:] :
# 		print("\n########################## ITERATION %s ###########################\n" % iteration)
# 		log_write.write("\n########################## ITERATION %s ###########################\n" % iteration)
	
# 		with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\random_networks\\Mock_networks_21Apr\\COVID19_GDDS_nodes_%s.csv' % str(iteration), 'r') as nodecsv: # Open the file                   
# 			nodereader = csv.reader(nodecsv) # Read the csv  
# 			nodes = [n for n in nodereader][1:]                     
# 			node_names = [n[0] for n in nodes] # Get a list of only the node names 

# 		with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\random_networks\\Mock_networks_21Apr\\COVID19_GDDS_edges_%s.csv' % str(iteration), 'r') as edgecsv: # Open the file
# 			edgereader = csv.reader(edgecsv) # Read the csv     
# 			edges = [tuple(e) for e in edgereader][1:] # Retrieve the data

# 			log_write.write("len(node_names) : %s\n" % len(node_names)) #16007 --> 13448
# 			log_write.write("len(edges) : %s\n" % len(edges)) #53917 --> 37835

# 		G = nx.Graph()
# 		G.add_nodes_from(node_names)
# 		G.add_edges_from(edges)

# 		descr_dict = {}
# 		description_set=set()
# 		for node in nodes: 
# 			descr_dict[node[0]] = node[2]
# 			description_set.add(node[2])

# 		# print(description_set) #{'Symptom', 'Human PPI (target)', 'GO', 'Disease', 'Drug'}
# 		# recounted = Counter(list(descr_dict.values()))
# 		# print("%s\n" % recounted)

# 		# print(nx.info(G))
# 		log_write.write("\nRANDOM GRAPH %s\n" % iteration)
# 		log_write.write("%s\n" % nx.info(G))

# 		#####################################################################################################
# 		################################# GET EDGES CONNECTING DRUGS TO NODES  OF INTEREST ################
# 		#####################################################################################################

########### I am forced to take all proteins because there are no first order proteins in the random networks #########

# 		#first we have to build the set of proteins of interest
# 		proteins_of_interest = set([n for n in G.nodes if descr_dict[n] == 'Human PPI (target)'])
# 		# print(proteins_of_interest)

# 		# now we keep only the drug that interact with those proteins of interest
# 		drugs_linked_to_some_proteins_of_interest = []
# 		for u,v,c in G.edges(data=True) :	
# 			if (descr_dict[u] == 'Drug' and v in proteins_of_interest) :
# 				drugs_linked_to_some_proteins_of_interest.append(u) # this never happens
# 			if (descr_dict[v] == 'Drug' and u in proteins_of_interest) :
# 				drugs_linked_to_some_proteins_of_interest.append(v)
# 		# print(len(drugs_linked_to_some_proteins_of_interest), len(set(drugs_linked_to_some_proteins_of_interest))) # 8200 3740
		
# 		drug_counts = sorted(drugs_linked_to_some_proteins_of_interest)

# 		# for the drugs
# 		drug_counts_dict = defaultdict( int )
# 		for drug in drug_counts:
# 		    drug_counts_dict[drug] += 1

# 		my_drugs_nodes = [n for n in G.nodes if descr_dict[n] == 'Drug']
# 		for drug in my_drugs_nodes :
# 			if drug_counts_dict.get(drug) == None :
# 				# print("%s is missing in the linked drugs : Adding it with a frequency of 0.\n" % drug)
# 				drug_counts_dict[drug]= 0

# 		average_degree_of_my_drugs_nodes = statistics.mean([ drug_counts_dict[drug] for drug in drug_counts_dict ])
# 		# print("average_degree_of_my_drugs_nodes", average_degree_of_my_drugs_nodes)
# 		log_write.write("Average degree of drugs to a viral interactor : %s\n" % average_degree_of_my_drugs_nodes)

# 		ordered_drug_counts_dict = OrderedDict(sorted(drug_counts_dict.items()))

# 		with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\Further_Analysis\\Random_Networks\\Drugs\\viral_interactors\\COVID19_GDDS_drugs_degrees_to_proteins_of_interest_%s.tsv' % str(iteration), 'a+') as file_write:
# 			for drug in ordered_drug_counts_dict :
# 				if (sum(ordered_drug_counts_dict.values()) != 0) :
# 					file_write.write("%s\t%s\t%s\n" % (drug, ordered_drug_counts_dict[drug], ordered_drug_counts_dict[drug] / sum(ordered_drug_counts_dict.values())))
# 				else :
# 					file_write.write("%s\t%s\t%s\n" % (drug, ordered_drug_counts_dict[drug], ordered_drug_counts_dict[drug]))



# # ##################################################################
# # ############ COMPUTE DEGREES MEANS AND Z-SCORES ##################
# # ######## INCLUDING ABSENCE OF DISEASES AS ZERO DEGREES ###########
# # ##################################################################


# # For covid network, compute the mean and the standard deviation of drugs degrees and normalized degrees
covid_drugs = []
covid_drug_degree_dict = {}
covid_drug_norm_degree_dict = {}

total_drug_degree_dict = defaultdict(list)
total_drug_normalized_degree_dict = defaultdict(list)

with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\Further_Analysis\\Drugs\\COVID19_GDDS_drugs_degrees_to_proteins_of_interest.tsv', 'r') as data :
	for line in data :
		line = line.split('\t')
		drug = line[0]
		covid_drugs.append(drug)
		drug_degree = int(line[1])
		drug_norm_degree = float(line[2].replace('\n', ''))
		if covid_drug_degree_dict.get(drug) == None :
			covid_drug_degree_dict[drug] = drug_degree
		else :
			print("warning : %s duplicate ?!" % drug)
		covid_drug_norm_degree_dict[drug] = drug_norm_degree

		if total_drug_degree_dict.get(drug) == None :
			total_drug_degree_dict[drug] = [drug_degree]
		else :
			total_drug_degree_dict[drug].append(drug_degree)

		if total_drug_normalized_degree_dict.get(drug) == None :
			total_drug_normalized_degree_dict[drug] = [drug_norm_degree]
		else :
			total_drug_normalized_degree_dict[drug].append(drug_norm_degree)

mean_drug_deg = statistics.mean(covid_drug_degree_dict.values())
mean_norm_drug_deg = statistics.mean(covid_drug_norm_degree_dict.values())
sd_drug_deg = statistics.pstdev(covid_drug_degree_dict.values())
sd_norm_drug_deg = statistics.pstdev(covid_drug_norm_degree_dict.values())
print(mean_drug_deg, sd_drug_deg, mean_norm_drug_deg, sd_norm_drug_deg) #1.4378397334736104 2.251958221466831 0.0001753463089601964 0.00027462905139839403

print(len(covid_drugs),len(covid_drug_degree_dict),len(covid_drug_norm_degree_dict)) #5703 5703 5703


nodefiles_numbers = []
edgefiles_numbers = []
for filename in os.listdir("A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\random_networks\\Mock_networks_21Apr"):
	if 'nodes_' in filename :
		filename = filename.split(".csv")[0]
		nodefiles_numbers.append(int(filename[19:]))
	if 'edges_' in filename :
		filename = filename.split(".csv")[0]
		edgefiles_numbers.append(int(filename[19:]))

# log_write.write("len(nodefiles_numbers) = %s\n" % len(nodefiles_numbers)) #2051
# log_write.write("len(edgefiles_numbers) = %s\n" % len(edgefiles_numbers)) #2051
# log_write.write("difference between edgefiles_numbers and nodefiles_numbers = %s\n" % len(set(nodefiles_numbers).difference(set(edgefiles_numbers)))) #0
print("len(nodefiles_numbers) = %s\n" % len(nodefiles_numbers)) #2051
print("len(edgefiles_numbers) = %s\n" % len(edgefiles_numbers)) #2051
print("difference between edgefiles_numbers and nodefiles_numbers = %s\n" % len(set(nodefiles_numbers).difference(set(edgefiles_numbers)))) #0

nodefiles_numbers = sorted(nodefiles_numbers)
edgefiles_numbers = sorted(edgefiles_numbers)


for iteration in nodefiles_numbers[1:] :
	# log_write.write("\n########################## ITERATION %s ###########################\n" % iteration)
	print("\n########################## ITERATION %s ###########################\n" % iteration)

	with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\Further_Analysis\\Random_Networks\\Drugs\\viral_interactors\\COVID19_GDDS_drugs_degrees_to_proteins_of_interest_%s.tsv' % str(iteration), 'r') as data: # Open the file  
		for line in data :
			line = line.split('\t')
			drug = line[0]
			degree = float(line[1])
			normalized_degree = float(line[2].replace('\n', ''))

			if total_drug_degree_dict.get(drug) == None :
				total_drug_degree_dict[drug] = [degree]
			else :
				total_drug_degree_dict[drug].append(degree)

			if total_drug_normalized_degree_dict.get(drug) == None :
				total_drug_normalized_degree_dict[drug] = [normalized_degree]
			else :
				total_drug_normalized_degree_dict[drug].append(normalized_degree)

print(len(total_drug_degree_dict)) # 72069
print(len(total_drug_normalized_degree_dict)) # 72069

### make the files with all calculations
with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\Further_Analysis\\Drugs\\COVID19_GDDS_Results_Drugs_viral_interactors-withZeros.tsv', 'a+') as file_write:
	file_write.write("drug\t# of nets it is found in\tmin-degree\tmax-degree\tmean-degree\tmedian-degree\tsd-degree\n")
	for drug in total_drug_degree_dict :
		nets_number = len(total_drug_degree_dict[drug])
		min_deg = min(total_drug_degree_dict[drug])
		max_deg = max(total_drug_degree_dict[drug])
		mean_deg = statistics.mean(total_drug_degree_dict[drug])
		median_deg = statistics.median(total_drug_degree_dict[drug])
		sd_deg = statistics.pstdev(total_drug_degree_dict[drug])
		file_write.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (drug, nets_number, min_deg, max_deg, mean_deg, median_deg, sd_deg))


with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\Further_Analysis\\Drugs\\COVID19_GDDS_Results_Drugs_viral_interactors_Relative-withZeros.tsv', 'a+') as file_write:
	file_write.write("drug\t# of nets it is found in\tmin-degree-rel\tmax-degree-rel\tmean-degree-rel\tmedian-degree-rel\tsd-degree-rel\n")
	for drug in total_drug_normalized_degree_dict :
		nets_number = len(total_drug_normalized_degree_dict[drug])
		min_deg = min(total_drug_normalized_degree_dict[drug])
		max_deg = max(total_drug_normalized_degree_dict[drug])
		mean_deg = statistics.mean(total_drug_normalized_degree_dict[drug])
		median_deg = statistics.median(total_drug_normalized_degree_dict[drug])
		sd_deg = statistics.pstdev(total_drug_normalized_degree_dict[drug])
		file_write.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (drug, nets_number, min_deg, max_deg, mean_deg, median_deg, sd_deg))

### compute Z-score for the covid diseases
# For each drug found in covid network, compute 2 Z-scores.

with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\Further_Analysis\\Drugs\\COVID19_GDDS_Results_ForCovidDrugs-viral_interactors-withZeros.tsv', 'a+') as file_write:
	file_write.write("covid drug\t# of nets it is found in\tmin-degree\tmax-degree\tmean-degree\tmedian-degree\tsd-degree\tcovid-degree\tZ-score\n")
	for drug in covid_drugs :
		nets_number = len(total_drug_degree_dict[drug])
		min_deg = min(total_drug_degree_dict[drug])
		max_deg = max(total_drug_degree_dict[drug])
		mean_deg = statistics.mean(total_drug_degree_dict[drug])
		median_deg = statistics.median(total_drug_degree_dict[drug])
		sd_deg = statistics.pstdev(total_drug_degree_dict[drug])
		covid_deg = covid_drug_degree_dict[drug]
		if sd_deg != 0 :
			z_score = ( covid_deg - mean_deg) / sd_deg
		else : 
			z_score = 'NaN'
		file_write.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (drug, nets_number, min_deg, max_deg, mean_deg, median_deg, sd_deg, covid_deg, z_score))


with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\Further_Analysis\\Drugs\\COVID19_GDDS_Results_ForCovidDrugs_viral_interactors_Relative-withZeros.tsv', 'a+') as file_write:
	file_write.write("covid drug\t# of nets it is found in\tmin-degree-rel\tmax-degree-rel\tmean-degree-rel\tmedian-degree-rel\tsd-degree-rel\tcovid-degree-rel\tZ-score-rel\n")
	for drug in covid_drugs :
		nets_number = len(total_drug_normalized_degree_dict[drug])
		min_deg = min(total_drug_normalized_degree_dict[drug])
		max_deg = max(total_drug_normalized_degree_dict[drug])
		mean_deg = statistics.mean(total_drug_normalized_degree_dict[drug])
		median_deg = statistics.median(total_drug_normalized_degree_dict[drug])
		sd_deg = statistics.pstdev(total_drug_normalized_degree_dict[drug])
		covid_deg = covid_drug_norm_degree_dict[drug]
		if sd_deg != 0 :
			z_score = ( covid_deg - mean_deg) / sd_deg
		else : 
			z_score = 'NaN'
		file_write.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (drug, nets_number, min_deg, max_deg, mean_deg, median_deg, sd_deg, covid_deg, z_score))










stop = timeit.default_timer()
print(stop - start)
 
random-network-analysis-drugs-ranking.py,
#!/usr/bin/python
# -*- coding: utf-8 -*-
# jupyter notebook --browser='C:/Program Files (x86)/Google/Chrome/Application/chrome.exe %s'


import timeit
start = timeit.default_timer()

import os
import csv
from operator import itemgetter
import networkx as nx
from networkx.algorithms import community
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
from collections import Counter
from collections import defaultdict
from collections import OrderedDict
import statistics 
import re
import gzip





#####################################################################################################
#########################  Compute network Drug degree distribution #################################
##########################################  based on all edges  #####################################
#####################################################################################################

with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\COVID19_GDDS_nodes.csv', 'r') as nodecsv: # Open the file                       
	nodereader = csv.reader(nodecsv) # Read the csv  
	nodes = [n for n in nodereader][1:]                     
	node_names = [n[0] for n in nodes] # Get a list of only the node names 

with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\COVID19_GDDS_edges.csv', 'r') as edgecsv: # Open the file
	edgereader = csv.reader(edgecsv) # Read the csv     
	edges = [tuple(e) for e in edgereader][1:] # Retrieve the data

G = nx.Graph()
G.add_nodes_from(node_names)
G.add_edges_from(edges)

descr_dict = {}
for node in nodes: 
	descr_dict[node[0]] = node[2]

print(nx.info(G))
# Number of nodes: 16007
# Number of edges: 53917
# Average degree:   6.7367

# print(nx.degree(G)) # this is the list containing as many tuples as nodes, indicating each node's degree.

average_degree = statistics.mean([ tpl[1] for tpl in nx.degree(G) ])
print(average_degree) # 6.736677703504717

drugs_degrees = nx.degree(G, nbunch=[n for n in G.nodes if descr_dict[n] == 'Drug'])
# print("degree of drugs nodes only", drugs_degrees)

average_degree_of_drugs = statistics.mean([ tpl[1] for tpl in drugs_degrees ])
print("average degree of drugs nodes", average_degree_of_drugs) # 2.2826582500438364

# #####################################################################################################
# ################################# CALCULATION ON EDGES WITHOUT FILTERING  ###########################
# #####################################################################################################
# # let's count edges only in order to calculate the nodes' degrees based on the number of edges they make in the graph

drugs_linked_to_some_nodes = []
for u,v,c in G.edges(data=True) :	
	if (descr_dict[u] == 'Drug') :
		drugs_linked_to_some_nodes.append(u) # this never happens
	if (descr_dict[v] == 'Drug') :
		drugs_linked_to_some_nodes.append(v)
print(len(drugs_linked_to_some_nodes), len(set(drugs_linked_to_some_nodes))) # 13018 5703
# ({'Drug': 5703, 'Disease': 4176, 'GO': 3487, 'Symptom': 2157, 'Human PPI (target)': 457, 'Viral Gene': 27})

# print(sorted(drugs_linked_to_some_nodes))

drug_counts = sorted(drugs_linked_to_some_nodes)

# # for the drugs
drug_counts_dict = defaultdict( int )
for drug in drug_counts:
    drug_counts_dict[drug] += 1

my_drugs_nodes = [n for n in G.nodes if descr_dict[n] == 'Drug']
for drug in my_drugs_nodes :
	if drug_counts_dict.get(drug) == None :
		print("%s is missing in the linked drugs : Adding it with a frequency of 0.\n" % drug)
		drug_counts_dict[drug]= 0

average_degree_of_my_drugs_nodes = statistics.mean([ drug_counts_dict[drug] for drug in drug_counts_dict ])
# print("average_degree_of_my_drugs_nodes", average_degree_of_my_drugs_nodes)
# 2.2826582500438364

ordered_drug_counts_dict = OrderedDict(sorted(drug_counts_dict.items()))

# with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\Further_Analysis\\Drugs\\COVID19_GDDS_drugs_degrees_to_all_nodes.tsv', 'a+') as file_write:
# 	for drug in ordered_drug_counts_dict :
# 		file_write.write("%s\t%s\t%s\n" % (drug, ordered_drug_counts_dict[drug], ordered_drug_counts_dict[drug] / sum(ordered_drug_counts_dict.values())))


# ##################################################################
# ############ SAME ANALYSIS FOR RANDOM NETWORKS ###################
# ##################################################################

# with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\Further_Analysis\\Drugs\\Log_degree_analysis_on_random_networks_for_drugs_all_nodes.txt', 'a+') as log_write:

# 	nodefiles_numbers = []
# 	edgefiles_numbers = []
# 	for filename in os.listdir("A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\random_networks\\Mock_networks_21Apr"):
# 		if 'nodes_' in filename :
# 			filename = filename.split(".csv")[0]
# 			nodefiles_numbers.append(int(filename[19:]))
# 		if 'edges_' in filename :
# 			filename = filename.split(".csv")[0]
# 			edgefiles_numbers.append(int(filename[19:]))

# 	log_write.write("len(nodefiles_numbers) = %s\n" % len(nodefiles_numbers)) #2051
# 	log_write.write("len(edgefiles_numbers) = %s\n" % len(edgefiles_numbers)) #2051
# 	log_write.write("difference between edgefiles_numbers and nodefiles_numbers = %s\n" % len(set(nodefiles_numbers).difference(set(edgefiles_numbers)))) #0

# 	nodefiles_numbers = sorted(nodefiles_numbers)
# 	edgefiles_numbers = sorted(edgefiles_numbers)

# 	for iteration in nodefiles_numbers[1:] :
# 		print("\n########################## ITERATION %s ###########################\n" % iteration)
# 		log_write.write("\n########################## ITERATION %s ###########################\n" % iteration)
	
# 		with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\random_networks\\Mock_networks_21Apr\\COVID19_GDDS_nodes_%s.csv' % str(iteration), 'r') as nodecsv: # Open the file                   
# 			nodereader = csv.reader(nodecsv) # Read the csv  
# 			nodes = [n for n in nodereader][1:]                     
# 			node_names = [n[0] for n in nodes] # Get a list of only the node names 

# 		with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\random_networks\\Mock_networks_21Apr\\COVID19_GDDS_edges_%s.csv' % str(iteration), 'r') as edgecsv: # Open the file
# 			edgereader = csv.reader(edgecsv) # Read the csv     
# 			edges = [tuple(e) for e in edgereader][1:] # Retrieve the data

# 			log_write.write("len(node_names) : %s\n" % len(node_names)) #16007 --> 13448
# 			log_write.write("len(edges) : %s\n" % len(edges)) #53917 --> 37835

# 		G = nx.Graph()
# 		G.add_nodes_from(node_names)
# 		G.add_edges_from(edges)

# 		descr_dict = {}
# 		for node in nodes: 
# 			descr_dict[node[0]] = node[2]

# 		print(nx.info(G))
# 		log_write.write("\nRANDOM GRAPH %s\n" % iteration)
# 		log_write.write("%s\n" % nx.info(G))

# 		#####################################################################################################
# 		################################# GET EDGES CONNECTING DRUGS TO ANY NODES  ################
# 		#####################################################################################################
# 		# We fetch all edges based on drug-nodes or nodes-drug links, count these edges, 
# 		# in order to calculate the nodes' degrees based on their number.

# 		drugs_linked_to_some_nodes = []
# 		for u,v,c in G.edges(data=True) :	
# 			if (descr_dict[u] == 'Drug') :
# 				drugs_linked_to_some_nodes.append(u)
# 			if (descr_dict[v] == 'Drug') :
# 				drugs_linked_to_some_nodes.append(v)
# 		print(len(drugs_linked_to_some_nodes), len(set(drugs_linked_to_some_nodes))) #
# 		log_write.write("Number of edges linking drugs to any nodes : %s \n" % len(drugs_linked_to_some_nodes))
# 		log_write.write("Number of drugs with at least one edge to a node: %s \n" % len(set(drugs_linked_to_some_nodes)))
		
# 		drug_counts = sorted(drugs_linked_to_some_nodes)

# 		# for the drugs
# 		drug_counts_dict = defaultdict( int )
# 		for drug in drug_counts:
# 		    drug_counts_dict[drug] += 1

# 		my_drugs_nodes = [n for n in G.nodes if descr_dict[n] == 'Drug']
# 		for drug in my_drugs_nodes :
# 			if drug_counts_dict.get(drug) == None :
# 				# print("%s is missing in the linked drugs : Adding it with a frequency of 0.\n" % drug)
# 				drug_counts_dict[drug]= 0

# 		average_degree_of_my_drugs_nodes = statistics.mean([ drug_counts_dict[drug] for drug in drug_counts_dict ])
# 		print("average_degree_of_my_drugs_nodes", average_degree_of_my_drugs_nodes)
# 		log_write.write("Average degree of drugs to any nodes : %s\n" % average_degree_of_my_drugs_nodes)

# 		ordered_drug_counts_dict = OrderedDict(sorted(drug_counts_dict.items()))

# 		with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\Further_Analysis\\Random_Networks\\Drugs\\COVID19_GDDS_drugs_degrees_to_all_nodes_%s.tsv' % str(iteration), 'a+') as file_write:
# 			for drug in ordered_drug_counts_dict :
# 				file_write.write("%s\t%s\t%s\n" % (drug, ordered_drug_counts_dict[drug], ordered_drug_counts_dict[drug] / sum(ordered_drug_counts_dict.values())))




# ##################################################################
# ############ COMPUTE DEGREES MEANS AND Z-SCORES ##################
# ######## INCLUDING ABSENCE OF DISEASES AS ZERO DEGREES ###########
# ##################################################################


# For covid network, compute the mean and the standard deviation of drugs degrees and normalized degrees
covid_drugs = []
covid_drug_degree_dict = {}
covid_drug_norm_degree_dict = {}

total_drug_degree_dict = defaultdict(list)
total_drug_normalized_degree_dict = defaultdict(list)

with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\Further_Analysis\\Drugs\\COVID19_GDDS_drugs_degrees_to_all_nodes.tsv', 'r') as data :
	for line in data :
		line = line.split('\t')
		drug = line[0]
		covid_drugs.append(drug)
		drug_degree = int(line[1])
		drug_norm_degree = float(line[2].replace('\n', ''))
		if covid_drug_degree_dict.get(drug) == None :
			covid_drug_degree_dict[drug] = drug_degree
		else :
			print("warning : %s duplicate ?!" % drug)
		covid_drug_norm_degree_dict[drug] = drug_norm_degree

		if total_drug_degree_dict.get(drug) == None :
			total_drug_degree_dict[drug] = [drug_degree]
		else :
			total_drug_degree_dict[drug].append(drug_degree)

		if total_drug_normalized_degree_dict.get(drug) == None :
			total_drug_normalized_degree_dict[drug] = [drug_norm_degree]
		else :
			total_drug_normalized_degree_dict[drug].append(drug_norm_degree)

mean_drug_deg = statistics.mean(covid_drug_degree_dict.values())
mean_norm_drug_deg = statistics.mean(covid_drug_norm_degree_dict.values())
sd_drug_deg = statistics.pstdev(covid_drug_degree_dict.values())
sd_norm_drug_deg = statistics.pstdev(covid_drug_norm_degree_dict.values())
print(mean_drug_deg, sd_drug_deg, mean_norm_drug_deg, sd_norm_drug_deg) #2.2826582500438364 2.9197086133682784 0.0001753463089601964 0.0002242824253624426

print(len(covid_drugs),len(covid_drug_degree_dict),len(covid_drug_norm_degree_dict)) #5703 5703 5703


nodefiles_numbers = []
edgefiles_numbers = []
for filename in os.listdir("A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\random_networks\\Mock_networks_21Apr"):
	if 'nodes_' in filename :
		filename = filename.split(".csv")[0]
		nodefiles_numbers.append(int(filename[19:]))
	if 'edges_' in filename :
		filename = filename.split(".csv")[0]
		edgefiles_numbers.append(int(filename[19:]))

# log_write.write("len(nodefiles_numbers) = %s\n" % len(nodefiles_numbers)) #2051
# log_write.write("len(edgefiles_numbers) = %s\n" % len(edgefiles_numbers)) #2051
# log_write.write("difference between edgefiles_numbers and nodefiles_numbers = %s\n" % len(set(nodefiles_numbers).difference(set(edgefiles_numbers)))) #0
print("len(nodefiles_numbers) = %s\n" % len(nodefiles_numbers)) #2051
print("len(edgefiles_numbers) = %s\n" % len(edgefiles_numbers)) #2051
print("difference between edgefiles_numbers and nodefiles_numbers = %s\n" % len(set(nodefiles_numbers).difference(set(edgefiles_numbers)))) #0

nodefiles_numbers = sorted(nodefiles_numbers)
edgefiles_numbers = sorted(edgefiles_numbers)


for iteration in nodefiles_numbers[1:] :
	# log_write.write("\n########################## ITERATION %s ###########################\n" % iteration)
	print("\n########################## ITERATION %s ###########################\n" % iteration)

	with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\Further_Analysis\\Random_Networks\\Drugs\\COVID19_GDDS_drugs_degrees_to_all_nodes_%s.tsv' % str(iteration), 'r') as data: # Open the file  
		for line in data :
			line = line.split('\t')
			drug = line[0]
			degree = float(line[1])
			normalized_degree = float(line[2].replace('\n', ''))

			if total_drug_degree_dict.get(drug) == None :
				total_drug_degree_dict[drug] = [degree]
			else :
				total_drug_degree_dict[drug].append(degree)

			if total_drug_normalized_degree_dict.get(drug) == None :
				total_drug_normalized_degree_dict[drug] = [normalized_degree]
			else :
				total_drug_normalized_degree_dict[drug].append(normalized_degree)

print(len(total_drug_degree_dict)) # 72069
print(len(total_drug_normalized_degree_dict)) # 72069

### make the files with all calculations
with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\Further_Analysis\\Drugs\\COVID19_GDDS_Results_Drugs_All_Nodes-withZeros.tsv', 'a+') as file_write:
	file_write.write("drug\t# of nets it is found in\tmin-degree\tmax-degree\tmean-degree\tmedian-degree\tsd-degree\n")
	for drug in total_drug_degree_dict :
		nets_number = len(total_drug_degree_dict[drug])
		min_deg = min(total_drug_degree_dict[drug])
		max_deg = max(total_drug_degree_dict[drug])
		mean_deg = statistics.mean(total_drug_degree_dict[drug])
		median_deg = statistics.median(total_drug_degree_dict[drug])
		sd_deg = statistics.pstdev(total_drug_degree_dict[drug])
		file_write.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (drug, nets_number, min_deg, max_deg, mean_deg, median_deg, sd_deg))


with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\Further_Analysis\\Drugs\\COVID19_GDDS_Results_Drugs_All_Nodes_Relative-withZeros.tsv', 'a+') as file_write:
	file_write.write("drug\t# of nets it is found in\tmin-degree-rel\tmax-degree-rel\tmean-degree-rel\tmedian-degree-rel\tsd-degree-rel\n")
	for drug in total_drug_normalized_degree_dict :
		nets_number = len(total_drug_normalized_degree_dict[drug])
		min_deg = min(total_drug_normalized_degree_dict[drug])
		max_deg = max(total_drug_normalized_degree_dict[drug])
		mean_deg = statistics.mean(total_drug_normalized_degree_dict[drug])
		median_deg = statistics.median(total_drug_normalized_degree_dict[drug])
		sd_deg = statistics.pstdev(total_drug_normalized_degree_dict[drug])
		file_write.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (drug, nets_number, min_deg, max_deg, mean_deg, median_deg, sd_deg))

### compute Z-score for the covid diseases
# For each drug found in covid network, compute 2 Z-scores.

with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\Further_Analysis\\Drugs\\COVID19_GDDS_Results_ForCovidDrugs-withZeros.tsv', 'a+') as file_write:
	file_write.write("covid drug\t# of nets it is found in\tmin-degree\tmax-degree\tmean-degree\tmedian-degree\tsd-degree\tcovid-degree\tZ-score\n")
	for drug in covid_drugs :
		nets_number = len(total_drug_degree_dict[drug])
		min_deg = min(total_drug_degree_dict[drug])
		max_deg = max(total_drug_degree_dict[drug])
		mean_deg = statistics.mean(total_drug_degree_dict[drug])
		median_deg = statistics.median(total_drug_degree_dict[drug])
		sd_deg = statistics.pstdev(total_drug_degree_dict[drug])
		covid_deg = covid_drug_degree_dict[drug]
		if sd_deg != 0 :
			z_score = ( covid_deg - mean_deg) / sd_deg
		else : 
			z_score = 'NaN'
		file_write.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (drug, nets_number, min_deg, max_deg, mean_deg, median_deg, sd_deg, covid_deg, z_score))


with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\Further_Analysis\\Drugs\\COVID19_GDDS_Results_ForCovidDrugs_Relative-withZeros.tsv', 'a+') as file_write:
	file_write.write("covid drug\t# of nets it is found in\tmin-degree-rel\tmax-degree-rel\tmean-degree-rel\tmedian-degree-rel\tsd-degree-rel\tcovid-degree-rel\tZ-score-rel\n")
	for drug in covid_drugs :
		nets_number = len(total_drug_normalized_degree_dict[drug])
		min_deg = min(total_drug_normalized_degree_dict[drug])
		max_deg = max(total_drug_normalized_degree_dict[drug])
		mean_deg = statistics.mean(total_drug_normalized_degree_dict[drug])
		median_deg = statistics.median(total_drug_normalized_degree_dict[drug])
		sd_deg = statistics.pstdev(total_drug_normalized_degree_dict[drug])
		covid_deg = covid_drug_norm_degree_dict[drug]
		if sd_deg != 0 :
			z_score = ( covid_deg - mean_deg) / sd_deg
		else : 
			z_score = 'NaN'
		file_write.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (drug, nets_number, min_deg, max_deg, mean_deg, median_deg, sd_deg, covid_deg, z_score))










stop = timeit.default_timer()
print(stop - start)
 
source-protein-degrees-undirected-corrected-generic-script.py,
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
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
from collections import Counter
from collections import defaultdict
from collections import OrderedDict
import statistics 
import re
import gzip


# dirname = "full_analysis_for_source_to_proteins_undirected"
# os.mkdir(dirname)
source_list = ['Symptom', 'Disease','Drug', 'Human PPI (target)']
# source_list = ['GO']

with open('Log_full_analysis_for_source_to_proteins_undirected.txt', 'a+') as full_analysis_log_write :

	#####################################################################################################
	#########################  Compute network source to protein degree distribution ####################
	##########################################  based on edges to protein ###############################
	#####################################################################################################

	with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\COVID19_GDDS_nodes.csv', 'r') as nodecsv: # Open the file                       
		nodereader = csv.reader(nodecsv) # Read the csv  
		nodes = [n for n in nodereader][1:]                     
		node_names = [n[0] for n in nodes] # Get a list of only the node names 

	with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\COVID19_GDDS_edges.csv', 'r') as edgecsv: # Open the file
		edgereader = csv.reader(edgecsv) # Read the csv     
		edges = [tuple(e) for e in edgereader][1:] # Retrieve the data

	G_covid = nx.Graph()
	G_covid.add_nodes_from(node_names)
	G_covid.add_edges_from(edges)

	descr_dict_covid = {}
	for node in nodes: 
		descr_dict_covid[node[0]] = node[2]

	print(nx.info(G_covid))
	# full_analysis_log_write.write(nx.info(G_covid))
	# Number of nodes: 16007
	# Number of edges: 53917
	# Average degree:   6.7367

	# print(nx.degree(G_covid)) # this is the list containing as many tuples as nodes, indicating each node's degree.

	average_degree = statistics.mean([ tpl[1] for tpl in nx.degree(G_covid) ])
	# print(average_degree) # 6.736677703504717

	for source_node in source_list :

		source_degrees = nx.degree(G_covid, nbunch=[n for n in G_covid.nodes if descr_dict_covid[n] == source_node])
		# print("degree of source nodes only", source_degrees)

		average_degree_of_sources = statistics.mean([ tpl[1] for tpl in source_degrees ])
		print("average degree of %s nodes" % source_node, average_degree_of_sources) # 
		# full_analysis_log_write.write("average degree of %s nodes : %s\n" % (source_node, average_degree_of_sources))
	# average degree of GO nodes 2.64123888729567
	# average degree of Symptom nodes 6.644413537320353
	# average degree of Disease nodes 4.329980842911877
	# average degree of Drug nodes 2.2826582500438364
	# average degree of Human PPI (target) nodes 115.66739606126914



	# # # # # #####################################################################################################
	# # # # # ################################# CALCULATION ON EDGES WITH FILTERING source - PROT #####################
	# # # # # #####################################################################################################

	for source_node in source_list :
		subfoldername = 'full_analysis_for_source_to_proteins_undirected/%s-prot' % (source_node)
		os.mkdir(subfoldername)
		print("Searching %s - Human PPI (target) links\n" % source_node)
		full_analysis_log_write.write("Searching %s - Human PPI (target) links\n" % source_node)
		
		# now we keep only the edges of source_node that interact with proteins
		sources_linked_to_proteins = []
		for u,v,c in G_covid.edges(data=True) :	
			if (descr_dict_covid[u] == source_node and descr_dict_covid[v] == 'Human PPI (target)') :
				sources_linked_to_proteins.append(u)
			if (descr_dict_covid[v] == source_node and descr_dict_covid[u] == 'Human PPI (target)') :
				sources_linked_to_proteins.append(v)	
		print("len(%s_linked_to_proteins), len(set(%s_linked_to_proteins))" % (source_node, source_node), len(sources_linked_to_proteins), len(set(sources_linked_to_proteins))) #
		full_analysis_log_write.write("len(%s_linked_to_proteins), len(set(%s_linked_to_proteins)) : %s, %s\n" % (source_node, source_node, len(sources_linked_to_proteins), len(set(sources_linked_to_proteins)))) #


		source_counts = sorted(sources_linked_to_proteins)

		source_counts_dict = defaultdict( int )
		for source in  source_counts:
		    source_counts_dict[source] += 1

		my_sources_nodes = [n for n in G_covid.nodes if descr_dict_covid[n] == source_node]
		for source in  my_sources_nodes :
			if source_counts_dict.get(source) == None :
				source_counts_dict[source]= 0

		average_degree_of_my_sources_nodes = statistics.mean([ source_counts_dict[source] for source in  source_counts_dict ])
		print("average_degree_of_my_sources_nodes", average_degree_of_my_sources_nodes)
		full_analysis_log_write.write("average_degree_of_my_%s_nodes : %s\n" % (source_node, average_degree_of_my_sources_nodes))

		max_degree_of_my_sources_nodes = max([ source_counts_dict[source] for source in  source_counts_dict ])
		print("max_degree_of_my_sources_nodes", max_degree_of_my_sources_nodes)

		ordered_source_counts_dict = OrderedDict(sorted(source_counts_dict.items()))


		with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\Further_Analysis\\full_analysis_for_source_to_proteins_undirected\\COVID19_GDDS_%s_degrees_to_proteins.tsv' % source_node, 'a+') as file_write:
			for source in ordered_source_counts_dict :
				file_write.write("%s\t%s\t%s\n" % (source, ordered_source_counts_dict[source], ordered_source_counts_dict[source] / sum(ordered_source_counts_dict.values())))


		# # # # # ##################################################################
		# # # # # ############ SAME ANALYSIS FOR RANDOM NETWORKS ###################
		# # # # # # ##################################################################
		
		with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\Further_Analysis\\full_analysis_for_source_to_proteins_undirected\\Log_degree_analysis_on_random_networks_for_%s_proteins.txt' % source_node, 'a+') as log_write:
			log_write.write("In covid net : len(sources_linked_to_proteins) and len(set(sources_linked_to_proteins)) are : %s and %s \n" % (len(sources_linked_to_proteins), len(set(sources_linked_to_proteins))) ) # # 9208 435
			log_write.write("In covid net : average_degree_of_my_sources_nodes = %s\n" % average_degree_of_my_sources_nodes)
			log_write.write("In covid net : max_degree_of_my_sources_nodes = %s\n" % max_degree_of_my_sources_nodes)

			nodefiles_numbers = []
			edgefiles_numbers = []
			for filename in os.listdir("A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\random_networks\\Mock_networks_21Apr"):
				if 'nodes_' in filename :
					filename = filename.split(".csv")[0]
					nodefiles_numbers.append(int(filename[19:]))
				if 'edges_' in filename :
					filename = filename.split(".csv")[0]
					edgefiles_numbers.append(int(filename[19:]))

			# log_write.write("len(nodefiles_numbers) = %s\n" % len(nodefiles_numbers)) #2051
			# log_write.write("len(edgefiles_numbers) = %s\n" % len(edgefiles_numbers)) #2051
			# log_write.write("difference between edgefiles_numbers and nodefiles_numbers = %s\n" % len(set(nodefiles_numbers).difference(set(edgefiles_numbers)))) #0

			nodefiles_numbers = sorted(nodefiles_numbers)
			edgefiles_numbers = sorted(edgefiles_numbers)

			networks_to_avoid = [3, 11, 22, 23, 33, 34, 38, 57, 72, 77, 123, 127, 148, 156, 157, 183, 205, 209, 210, 217, 232, 252, 268, 281, 291, 312, 341, 353, 362, 364, 413, 421, 433, 443, 451, 494, 525, 532, 569, 592, 629, 658, 662, 674, 676, 701, 705, 712, 736, 757, 788, 796, 825, 827, 838, 873, 876, 908, 910, 916, 917, 918, 938, 958, 965, 1088, 1125, 1128, 1145, 1146, 1165, 1171, 1183, 1185, 1192, 1195, 1211, 1244, 1254, 1257, 1300, 1353, 1354, 1364, 1372, 1394, 1395, 1401, 1410, 1427, 1454, 1477, 1497, 1517, 1532, 1549, 1551, 1575, 1578, 1611, 1624, 1629, 1637, 1650, 1655, 1697, 1700, 1713, 1736, 1745, 1762, 1790, 1800, 1815, 1824, 1827, 1866, 1872, 1886, 1887, 2053, 2072, 2106, 2115, 2128, 2160, 2217, 2343, 2383, 2388, 2413, 2458, 2465]
			nodefiles_numbers = [item for item in nodefiles_numbers if item not in networks_to_avoid]

			for iteration in nodefiles_numbers[1:] :
				print("\n########################## ITERATION %s -- %s ###########################\n" % (iteration, source_node))
				log_write.write("\n########################## ITERATION %s ###########################\n" % iteration)
			
				with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\random_networks\\Mock_networks_21Apr\\COVID19_GDDS_nodes_%s.csv' % str(iteration), 'r') as nodecsv: # Open the file                   
					nodereader = csv.reader(nodecsv) # Read the csv  
					nodes = [n for n in nodereader][1:]                     
					node_names = [n[0] for n in nodes] # Get a list of only the node names 

				with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\random_networks\\Mock_networks_21Apr\\COVID19_GDDS_edges_%s.csv' % str(iteration), 'r') as edgecsv: # Open the file
					edgereader = csv.reader(edgecsv) # Read the csv     
					edges = [tuple(e) for e in edgereader][1:] # Retrieve the data

					log_write.write("len(node_names) : %s\n" % len(node_names)) #16007 --> 13448
					log_write.write("len(edges) : %s\n" % len(edges)) #53917 --> 37835

				G_random = nx.Graph()
				G_random.add_nodes_from(node_names)
				G_random.add_edges_from(edges)

				descr_dict_random = {}
				description_random_set=set()
				for node in nodes: 
					descr_dict_random[node[0]] = node[2]
					description_random_set.add(node[2])

				# print(nx.info(G_random))
				log_write.write("\nRANDOM GRAPH %s\n" % iteration)
				log_write.write("%s\n" % nx.info(G_random))

				#####################################################################################################
				################################# GET EDGES CONNECTING sources TO proteins ################
				#####################################################################################################

				sources_linked_to_proteins = []
				for u,v,c in G_random.edges(data=True) :	
					if (descr_dict_random[u] == source_node and descr_dict_random[v] == 'Human PPI (target)') :
						sources_linked_to_proteins.append(u)
					if (descr_dict_random[v] == source_node and descr_dict_random[u] == 'Human PPI (target)') :
						sources_linked_to_proteins.append(v)		
				
				source_counts = sorted(sources_linked_to_proteins)

				source_counts_dict = defaultdict( int )
				for source in  source_counts:
				    source_counts_dict[source] += 1

				my_sources_nodes = [n for n in G_random.nodes if descr_dict_random[n] == source]
				for source in  my_sources_nodes :
					if source_counts_dict.get(source) == None :
						# print("%s is missing in the linked sources : Adding it with a frequency of 0.\n" % protein)
						source_counts_dict[source]= 0

				average_degree_of_my_sources_nodes = statistics.mean([ source_counts_dict[source] for source in  source_counts_dict ])
				log_write.write("Average degree of %s linked to Human PPI (target) : %s\n" % (source_node, average_degree_of_my_sources_nodes))

				ordered_source_counts_dict = OrderedDict(sorted(source_counts_dict.items()))

				with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\Further_Analysis\\full_analysis_for_source_to_proteins_undirected\\%s-prot\\COVID19_GDDS_%s_degrees_to_proteins_%s.tsv' % (source_node, source_node, str(iteration)), 'a+') as file_write:
					for source in ordered_source_counts_dict :
						if (sum(ordered_source_counts_dict.values()) != 0) :
							file_write.write("%s\t%s\t%s\n" % (source, ordered_source_counts_dict[source], ordered_source_counts_dict[source] / sum(ordered_source_counts_dict.values())))
						else :
							file_write.write("%s\t%s\t%s\n" % (source, ordered_source_counts_dict[source], ordered_source_counts_dict[source]))



		# # # # ##################################################################
		# # # # ############ COMPUTE DEGREES MEANS AND Z-SCORES ##################
		# # # # ######## INCLUDING ABSENCE OF source AS ZERO DEGREES ###########
		# # # # ##################################################################


		# # For covid network, compute the mean and the standard deviation of proteins degrees and normalized degrees
		covid_sources = []
		covid_source_degree_dict = {}
		covid_source_norm_degree_dict = {}

		total_source_degree_dict = defaultdict(list)
		total_source_normalized_degree_dict = defaultdict(list)

		with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\Further_Analysis\\full_analysis_for_source_to_proteins_undirected\\COVID19_GDDS_%s_degrees_to_proteins.tsv' % source_node, 'r') as data :
			for line in data :
				line = line.split('\t')
				source = line[0]
				covid_sources.append(source)
				source_degree = int(line[1])
				source_norm_degree = float(line[2].replace('\n', ''))
				if covid_source_degree_dict.get(source) == None :
					covid_source_degree_dict[source] = source_degree
				else :
					print("warning : %s duplicate ?!" % protein)
				covid_source_norm_degree_dict[source] = source_norm_degree

				if total_source_degree_dict.get(source) == None :
					total_source_degree_dict[source] = [source_degree]
				else :
					total_source_degree_dict[source].append(source_degree)

				if total_source_normalized_degree_dict.get(source) == None :
					total_source_normalized_degree_dict[source] = [source_norm_degree]
				else :
					total_source_normalized_degree_dict[source].append(source_norm_degree)

		mean_source_deg = statistics.mean(covid_source_degree_dict.values())
		mean_norm_source_deg = statistics.mean(covid_source_norm_degree_dict.values())
		sd_source_deg = statistics.pstdev(covid_source_degree_dict.values())
		sd_norm_source_deg = statistics.pstdev(covid_source_norm_degree_dict.values())
		print(mean_source_deg, sd_source_deg, mean_norm_source_deg, sd_norm_source_deg) #
		print(len(covid_sources),len(covid_source_degree_dict),len(covid_source_norm_degree_dict)) #
		# 20.14879649890591 18.9696222570052 0.002188183807439825 0.00206012405050013
		# 457 457 457






		nodefiles_numbers = []
		edgefiles_numbers = []
		for filename in os.listdir("A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\random_networks\\Mock_networks_21Apr"):
			if 'nodes_' in filename :
				filename = filename.split(".csv")[0]
				nodefiles_numbers.append(int(filename[19:]))
			if 'edges_' in filename :
				filename = filename.split(".csv")[0]
				edgefiles_numbers.append(int(filename[19:]))

		# log_write.write("len(nodefiles_numbers) = %s\n" % len(nodefiles_numbers)) #2051
		# log_write.write("len(edgefiles_numbers) = %s\n" % len(edgefiles_numbers)) #2051
		# log_write.write("difference between edgefiles_numbers and nodefiles_numbers = %s\n" % len(set(nodefiles_numbers).difference(set(edgefiles_numbers)))) #0
		print("len(nodefiles_numbers) = %s\n" % len(nodefiles_numbers)) #2051
		print("len(edgefiles_numbers) = %s\n" % len(edgefiles_numbers)) #2051
		print("difference between edgefiles_numbers and nodefiles_numbers = %s\n" % len(set(nodefiles_numbers).difference(set(edgefiles_numbers)))) #0

		# nodefiles_numbers = sorted(nodefiles_numbers)
		# edgefiles_numbers = sorted(edgefiles_numbers)

		networks_to_avoid = [3, 11, 22, 23, 33, 34, 38, 57, 72, 77, 123, 127, 148, 156, 157, 183, 205, 209, 210, 217, 232, 252, 268, 281, 291, 312, 341, 353, 362, 364, 413, 421, 433, 443, 451, 494, 525, 532, 569, 592, 629, 658, 662, 674, 676, 701, 705, 712, 736, 757, 788, 796, 825, 827, 838, 873, 876, 908, 910, 916, 917, 918, 938, 958, 965, 1088, 1125, 1128, 1145, 1146, 1165, 1171, 1183, 1185, 1192, 1195, 1211, 1244, 1254, 1257, 1300, 1353, 1354, 1364, 1372, 1394, 1395, 1401, 1410, 1427, 1454, 1477, 1497, 1517, 1532, 1549, 1551, 1575, 1578, 1611, 1624, 1629, 1637, 1650, 1655, 1697, 1700, 1713, 1736, 1745, 1762, 1790, 1800, 1815, 1824, 1827, 1866, 1872, 1886, 1887, 2053, 2072, 2106, 2115, 2128, 2160, 2217, 2343, 2383, 2388, 2413, 2458, 2465]
		nodefiles_numbers = [item for item in nodefiles_numbers if item not in networks_to_avoid]

		for iteration in nodefiles_numbers[1:] :
			# log_write.write("\n########################## ITERATION %s ###########################\n" % iteration)
			print("\n########################## ITERATION %s ###########################\n" % iteration)

			with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\Further_Analysis\\full_analysis_for_source_to_proteins_undirected\\%s-prot\\COVID19_GDDS_%s_degrees_to_proteins_%s.tsv' % (source_node, source_node, str(iteration)), 'r') as data: # Open the file  
				for line in data :
					line = line.split('\t')
					source = line[0]
					degree = float(line[1])
					normalized_degree = float(line[2].replace('\n', ''))

					if total_source_degree_dict.get(source) == None :
						total_source_degree_dict[source] = [degree]
					else :
						total_source_degree_dict[source].append(degree)

					if total_source_normalized_degree_dict.get(source) == None :
						total_source_normalized_degree_dict[source] = [normalized_degree]
					else :
						total_source_normalized_degree_dict[source].append(normalized_degree)

		print(len(total_source_degree_dict)) # 19928

		print(len(total_source_normalized_degree_dict)) # 19928

		### make the files with all calculations
		with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\Further_Analysis\\full_analysis_for_source_to_proteins_undirected\\COVID19_GDDS_Results_%s_to_Proteins-withZeros.tsv' % source_node, 'a+') as file_write:
			file_write.write("%s\t# of nets it is found in\tmin-degree\tmax-degree\tmean-degree\tmedian-degree\tsd-degree\n" % source_node)
			for source in  total_source_degree_dict :
				nets_number = len(total_source_degree_dict[source])
				min_deg = min(total_source_degree_dict[source])
				max_deg = max(total_source_degree_dict[source])
				mean_deg = statistics.mean(total_source_degree_dict[source])
				median_deg = statistics.median(total_source_degree_dict[source])
				sd_deg = statistics.pstdev(total_source_degree_dict[source])
				file_write.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (source, nets_number, min_deg, max_deg, mean_deg, median_deg, sd_deg))


		with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\Further_Analysis\\full_analysis_for_source_to_proteins_undirected\\COVID19_GDDS_Results_%s_to_Proteins_Relative-withZeros.tsv' % source_node, 'a+') as file_write:
			file_write.write("%s\t# of nets it is found in\tmin-degree-rel\tmax-degree-rel\tmean-degree-rel\tmedian-degree-rel\tsd-degree-rel\n" % source_node)
			for source in  total_source_normalized_degree_dict :
				nets_number = len(total_source_normalized_degree_dict[source])
				min_deg = min(total_source_normalized_degree_dict[source])
				max_deg = max(total_source_normalized_degree_dict[source])
				mean_deg = statistics.mean(total_source_normalized_degree_dict[source])
				median_deg = statistics.median(total_source_normalized_degree_dict[source])
				sd_deg = statistics.pstdev(total_source_normalized_degree_dict[source])
				file_write.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (source, nets_number, min_deg, max_deg, mean_deg, median_deg, sd_deg))

		### compute Z-score for the covid sources
		# For each source found in covid network, compute 2 Z-scores.

		with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\Further_Analysis\\full_analysis_for_source_to_proteins_undirected\\COVID19_GDDS_Results_%s_to_CovidProteins-withZeros.tsv' % source_node, 'a+') as file_write:
			file_write.write("covid %s\t# of nets it is found in\tmin-degree\tmax-degree\tmean-degree\tmedian-degree\tsd-degree\tcovid-degree\tZ-score\n" % source_node)
			for source in  covid_sources :
				nets_number = len(total_source_degree_dict[source])
				min_deg = min(total_source_degree_dict[source])
				max_deg = max(total_source_degree_dict[source])
				mean_deg = statistics.mean(total_source_degree_dict[source])
				median_deg = statistics.median(total_source_degree_dict[source])
				sd_deg = statistics.pstdev(total_source_degree_dict[source])
				covid_deg = covid_source_degree_dict[source]
				if sd_deg != 0 :
					z_score = ( covid_deg - mean_deg) / sd_deg
				else : 
					z_score = 'NaN'
				file_write.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (source, nets_number, min_deg, max_deg, mean_deg, median_deg, sd_deg, covid_deg, z_score))


		with open('A:\\Downloads\\Projects\\workFromHome\\Projects\\Covid\\Further_Analysis\\full_analysis_for_source_to_proteins_undirected\\COVID19_GDDS_Results_%s_to_CovidProteins_Relative-withZeros.tsv' % source_node, 'a+') as file_write:
			file_write.write("covid %s\t# of nets it is found in\tmin-degree-rel\tmax-degree-rel\tmean-degree-rel\tmedian-degree-rel\tsd-degree-rel\tcovid-degree-rel\tZ-score-rel\n" % source_node)
			for source in  covid_sources :
				nets_number = len(total_source_normalized_degree_dict[source])
				min_deg = min(total_source_normalized_degree_dict[source])
				max_deg = max(total_source_normalized_degree_dict[source])
				mean_deg = statistics.mean(total_source_normalized_degree_dict[source])
				median_deg = statistics.median(total_source_normalized_degree_dict[source])
				sd_deg = statistics.pstdev(total_source_normalized_degree_dict[source])
				covid_deg = covid_source_norm_degree_dict[source]
				if sd_deg != 0 :
					z_score = ( covid_deg - mean_deg) / sd_deg
				else : 
					z_score = 'NaN'
				file_write.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (source, nets_number, min_deg, max_deg, mean_deg, median_deg, sd_deg, covid_deg, z_score))










stop = timeit.default_timer()
print(stop - start)
  
test_normal_distributions.py,
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

transform_Zscore_in_pvalue.py,
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

validation.py,
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
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
from collections import Counter
from collections import defaultdict
from collections import OrderedDict
import statistics 
import re
import gzip
import math


# covid_disease_list = []
# with open('COVID19_GDDS_Results_Diseases_to_Proteins_rel.tsv', 'r') as file_read :
# 	data = file_read.readlines()
# 	for line in data[1:] :
# 		line = line.split('\t')
# 		covid_disease = line[0]
# 		covid_disease_list.append(covid_disease)
# print(len(set(covid_disease_list)))	#4176	
# covid_disease_set = set(covid_disease_list)


# top_disease_list = []
# with open('ranked_diseases.tsv', 'r') as file_read :
# 	data = file_read.readlines()
# 	for line in data[1:359] :
# 		line = line.split('\t')
# 		covid_disease = line[0]
# 		top_disease_list.append(covid_disease)
# print(len(set(top_disease_list)))	#358	
# top_disease_set = set(top_disease_list)


# zhou_disease_list = []
# with open('ZhouSymptoms_41467_2014_BFncomms5212_MOESM1045_ESM.txt', 'r') as file_read :
# 	data = file_read.readlines()
# 	for line in data[1:] :
# 		line = line.split('\t')
# 		zhou_disease = line[1]
# 		zhou_disease_list.append(zhou_disease)
# print(len(set(zhou_disease_list)))	#4219
# zhou_disease_set = set(zhou_disease_list)

# overlap_diseases = covid_disease_set.intersection(zhou_disease_set)
# print(len(overlap_diseases)) #940

# overlap_top_diseases = zhou_disease_set.intersection(top_disease_set)
# print(len(overlap_top_diseases)) #84
# print(overlap_top_diseases) #84






# covid_symptom_list = []
# with open('COVID19_GDDS_Results_Symptoms_to_Proteins_rel.tsv', 'r') as file_read :
# 	data = file_read.readlines()
# 	for line in data[1:] :
# 		line = line.split('\t')
# 		covid_symptom = line[0]
# 		covid_symptom_list.append(covid_symptom)
# print(len(set(covid_symptom_list)))	#2157	
# covid_symptom_set = set(covid_symptom_list)

# zhou_symptom_list = []
# with open('ZhouSymptoms_41467_2014_BFncomms5212_MOESM1045_ESM.txt', 'r') as file_read :
# 	data = file_read.readlines()
# 	for line in data[1:] :
# 		line = line.split('\t')
# 		zhou_symptom = line[0]
# 		zhou_symptom_list.append(zhou_symptom)
# print(len(set(zhou_symptom_list)))	#322
# zhou_symptom_set = set(zhou_symptom_list)

# overlap_symptoms = covid_symptom_set.intersection(zhou_symptom_set)
# print(len(overlap_symptoms)) #60

# overlap_symptoms = zhou_symptom_set.intersection(covid_symptom_set)
# print(len(overlap_symptoms)) #60





covid_disease_list = []
covid_disease_scores = {}
with open('COVID19_GDDS_Results_Diseases_to_Proteins_rel.tsv', 'r') as file_read :
	data = file_read.readlines()
	for line in data[1:] :
		line = line.split('\t')
		covid_disease = line[0]
		score = float(line[-1].replace('\n',''))
		covid_disease_list.append(covid_disease)
		covid_disease_scores[covid_disease] = score
print(len(set(covid_disease_list)))	#4176	
covid_disease_set = set(covid_disease_list)


# covid_disease_scores_sorted = {k: v for k, v in sorted(covid_disease_scores.items(), key=lambda item: item[1])}
# print(covid_disease_scores_sorted)

covid_disease_scores_filtered = {k: v for k, v in sorted(covid_disease_scores.items(), key=lambda item: item[1]) if v > 2}
# print(covid_disease_scores_filtered)
print(len(covid_disease_scores_filtered)) #358 > 1 and 86 > 2, lets go with this one.
top_covid_diseases = [k for k in covid_disease_scores_filtered]
print(top_covid_diseases)
top_covid_diseases_set = set(top_covid_diseases)

with open("top_disease.tsv", 'a+') as file_write :
	for disease in covid_disease_scores_filtered :
		file_write.write("%s\t%s\n" % (disease,covid_disease_scores_filtered[disease])) 










# zhou_disease_list = []
# zhou_disease_symptom_dict = {}
# with open('ZhouSymptoms_41467_2014_BFncomms5212_MOESM1045_ESM.txt', 'r') as file_read :
# 	data = file_read.readlines()
# 	for line in data[1:] :
# 		line = line.split('\t')
# 		zhou_symptom = line[0]
# 		zhou_disease = line[1]
# 		zhou_disease_list.append(zhou_disease)
# 		if zhou_disease not in zhou_disease_symptom_dict :
# 			zhou_disease_symptom_dict[zhou_disease] = []
# 		zhou_disease_symptom_dict[zhou_disease].append(zhou_symptom)

# print(len(set(zhou_disease_list)))	#4219
# zhou_disease_set = set(zhou_disease_list)

# overlap_diseases = top_covid_diseases_set.intersection(zhou_disease_set)
# print(len(overlap_diseases)) #21
# print(overlap_diseases)
# #{'Brain Edema', 'Muscle Spasticity', 'Cerebral Hemorrhage', 'Nijmegen Breakage Syndrome', 'Pseudoxanthoma Elasticum', 'Contracture', 'Blindness, Cortical', 'Airway Obstruction', 'Learning Disorders', 'Hyperbilirubinemia, Neonatal', 'Mitochondrial Encephalomyopathies', 'Hepatitis C', 'Lymphohistiocytosis, Hemophagocytic', 'Neoplasm, Residual', 'Ehlers-Danlos Syndrome', 'Infection', 'Tricuspid Valve Prolapse', 'Severe Acute Respiratory Syndrome', 'Liver Failure', 'Retroviridae Infections', 'Adrenal Gland Neoplasms'}
# disease_overlap_dict = {}
# for disease in overlap_diseases :
# 	disease_overlap_dict[disease] = covid_disease_scores[disease]
# disease_overlap_dict_sorted = {k: v for k, v in sorted(disease_overlap_dict.items(), key=lambda item: item[1])}

# symptoms_occurrences = []
# for disease in disease_overlap_dict_sorted :
# 	#print(disease, disease_overlap_dict_sorted[disease], zhou_disease_symptom_dict[disease])
# 	symptoms_occurrences.append(zhou_disease_symptom_dict[disease])

# # print(symptoms_occurrences)

# symptoms_counts = []
# for symptoms_list in symptoms_occurrences :
# 	for symptom in symptoms_list :
# 		symptoms_counts.append(symptom)
# # print(symptoms_counts)

# symptoms_counts = sorted(symptoms_counts)
# # print(symptoms_counts)

# symptoms_counts_dict = defaultdict( int )
# for symptom in  symptoms_counts:
#     symptoms_counts_dict[symptom] += 1
# symptoms_counts_dict_sorted = {k: v for k, v in sorted(symptoms_counts_dict.items(), key=lambda item: item[1])}
# for symptom in symptoms_counts_dict_sorted :
# 	print(symptom, symptoms_counts_dict_sorted[symptom])


stop = timeit.default_timer()
print(stop - start)
  
