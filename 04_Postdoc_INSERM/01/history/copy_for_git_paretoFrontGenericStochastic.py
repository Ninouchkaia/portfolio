
import math
import sys
import pylab as plt
import numpy as np

inputFile = sys.argv[1]
outputFile = sys.argv[2]

with open(inputFile, 'r') as file_read :
# with open("A:\\Downloads\\Projects\\workFromHome\\Projects\\ABM2021\\explorations\\outputs_ABM_9_2.txt", 'r') as file_read :
# with open("A:\\Downloads\\Projects\\workFromHome\\Projects\\ABM2021\\explorations\\ABM_8_4\\ABM_8_4\\population20000.csv", 'r') as file_read :

	my_points = []
	data = file_read.readlines()
	for line in data[1:] :
		line = line.replace("\n","").split(",")
		apo = int(line[0])
		needSig = int(line[1])
		layers = int(line[2])
		alpha = int(line[3])
		monoPhago = int(line[4])
		NLCPhago = int(line[5])
		M2Phago = int(line[6])
		M2Kill = int(line[7])
		cllDist = int(line[8])
		MonoDist = int(line[9])
		nlcDist = int(line[10])
		macroDist = int(line[11])
		nlcThreshold = int(line[12])
		signalInitMean = int(line[13])
		signalInitStd = int(line[14])
		diffTime = int(line[15])
		diffInitStd = int(line[16])
		gammaLifeInit = int(line[17])
		alphaDistrib = float(line[18])
		delta_fitness_via = float(line[19])
		delta_fitness_conc = float(line[20])
		euclid = math.sqrt(delta_fitness_via*delta_fitness_via + delta_fitness_conc * delta_fitness_conc)
		my_points.append([delta_fitness_via,delta_fitness_conc, euclid, 1, 
			apo, needSig, layers, alpha, 
			monoPhago, NLCPhago, M2Phago, M2Kill, 
			cllDist, MonoDist, nlcDist, macroDist, nlcThreshold, signalInitMean, signalInitStd,
			diffTime, diffInitStd, gammaLifeInit, alphaDistrib])
# print(len(my_points))

b_set = set(tuple(x) for x in my_points)
my_points_singles = [ list(x) for x in b_set ]
# print(len(my_points_singles))


my_points_sorted = sorted(my_points_singles, key=lambda x: x[2])

for i in range(0, len(my_points_sorted)) :
	x1 = my_points_sorted[i][0]
	y1 = my_points_sorted[i][1]

	for j in range(0, len(my_points_sorted)) :
		x2 = my_points_sorted[j][0]
		y2 = my_points_sorted[j][1]
		if not ((x1 == x2) and (y1 == y2)) :
			if (x2 <= x1) and (y2 <= y1) :
				my_points_sorted[i][3] = 0
				break

pareto_front = []
for point in my_points_sorted :
	if point[3] == 1 :
		pareto_front.append((point[0],point[1], point[4:]))


fig, ax = plt.subplots()

x_val = [x[0] for x in my_points_sorted]
y_val = [x[1] for x in my_points_sorted]
x_val_pareto = [x[0] for x in pareto_front]
y_val_pareto = [x[1] for x in pareto_front]

# plt.scatter(x_val, y_val,color='black',s=0.75)
# plt.scatter(x_val_pareto, y_val_pareto,color='red',s=0.75)

plt.scatter(x_val, y_val,color='black',s=3)
plt.scatter(x_val_pareto, y_val_pareto,color='red',s=3)
ax.set_xlabel(r'$\Delta Viability Fitness$', fontsize=15)
ax.set_ylabel(r'$\Delta ConcentrationFitness$', fontsize=15)
ax.set_title('Pareto front')
ax.grid(True)
plt.tight_layout()
plt.savefig('%s.png' % outputFile, bbox_inches='tight')

# plt.show()

print("len(pareto_front)",len(pareto_front))

pareto_front_sorted = sorted(pareto_front, key=lambda x: x[0])
with open("%s.txt" % outputFile, 'w') as file_write :
	file_write.write("delta_fitness_via,delta_fitness_conc, apo, needSig, layers, alpha, monoPhago, NLCPhago, M2Phago, M2Kill, cllDist, MonoDist, nlcDist, macroDist, nlcThreshold, signalInitMean, signalInitStd, diffTime, diffInitStd, LifeInitGamma, alphaDistrib\n")
	for sets in pareto_front_sorted :
		line = (",".join(str(x) for x in sets)).replace("[","").replace("]","")
		file_write.write(line)
		file_write.write("\n")
 


