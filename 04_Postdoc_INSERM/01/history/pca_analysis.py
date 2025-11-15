#!/usr/bin/python
# -*- coding: utf-8 -*-
# jupyter notebook --browser='C:/Program Files (x86)/Google/Chrome/Application/chrome.exe %s'


import timeit
start = timeit.default_timer()

import sys
import matplotlib.pyplot as plt

# import os

#### usage : python scripts/pca_analysis.py best_via_sets_all_patients.tsv

my_file = sys.argv[1] # best_via_sets_all_patients.tsv
# output = sys.argv[2]

import pandas as pd

df = pd.read_csv(my_file, sep='\t', header=0, index_col=0)
print(df)

my_data = df.T
my_data.reset_index(inplace=True)
my_data.rename(columns = {'index':'patient'}, inplace = True)
my_data.drop(["delta_fitness_via", "delta_fitness_conc"], axis = 1, inplace = True)

print(my_data)


from sklearn.preprocessing import StandardScaler
my_features = list(my_data)
my_features = my_features[1:]
print(my_features)

features = my_features

# features = []
# for feature in my_features :
# 	feature = feature.replace(" ", "")
# 	features.append(feature)

# Separating out the features
x = my_data.loc[:, features].values

# Separating out the target
y = my_data.loc[:,['patient']].values

# Standardizing the features
x = StandardScaler().fit_transform(x)


from sklearn.decomposition import PCA
pca = PCA(n_components=2)
principalComponents = pca.fit_transform(x)
principalDf = pd.DataFrame(data = principalComponents, columns = ['principal component 1', 'principal component 2'])


finalDf = pd.concat([principalDf, my_data[['patient']]], axis = 1)




fig = plt.figure(figsize = (8,8))
ax = fig.add_subplot(1,1,1) 
ax.set_xlabel('Principal Component 1', fontsize = 15)
ax.set_ylabel('Principal Component 2', fontsize = 15)


patients = my_data.patient.values.tolist()

colors = ['#00FFFF','#458B74','#000000','#0000FF','#8A2BE2','#FF4040', '#0FE721', '#FF0000', '#FC00FF', '#00FFB6', '#8B2F00']

for patient, color in zip(patients,colors):
    indicesToKeep = finalDf['patient'] == patient
    ax.scatter(finalDf.loc[indicesToKeep, 'principal component 1']
               , finalDf.loc[indicesToKeep, 'principal component 2']
               , c = color
               , s = 50)
ax.legend(patients)
ax.grid()

# plt.show()
ax.set_title('2-component PCA - %s\nExplained Variance = %s' % (my_file[:-17],pca.explained_variance_ratio_), fontsize = 20)

plt.savefig('pca_analysis_%s.pdf' % my_file[:-4])
plt.savefig('pca_analysis_%s.png' % my_file[:-4])

stop = timeit.default_timer()
print(stop - start)  


