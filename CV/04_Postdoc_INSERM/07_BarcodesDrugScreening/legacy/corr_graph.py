import timeit
start = timeit.default_timer()

import pandas as pd
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import seaborn as sns
 
df = pd.read_csv("result6_fillna_control_renamed_filtered6.csv", sep=';', header=0, index_col=0)
print (df) #

matrix = df.corr(
    method = 'pearson',  # The method of correlation
    min_periods = 1).round(2)     # Min number of observations required)

# # Transform it in a links data frame (3 columns only):
links = matrix.stack().reset_index()
print(links)
links.columns = ['var1', 'var2', 'value']

# Keep only correlation over a threshold and remove self correlation (cor(A,A)=1)
links_filtered=links.loc[ (links['value'] > 0.6) & (links['var1'] != links['var2']) ]
 
# Build your graph
G=nx.from_pandas_edgelist(links_filtered, 'var1', 'var2')
 
# Plot the network:
nx.draw(G, with_labels=True, node_color='orange', node_size=50, edge_color='black', linewidths=1, font_size=5)
# nx.draw(G, with_labels=False, node_color='orange', node_size=50, edge_color='black', linewidths=1, font_size=5)
plt.show() 

stop = timeit.default_timer()
print(stop - start)  