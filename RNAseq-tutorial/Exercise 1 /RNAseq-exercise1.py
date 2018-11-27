#!/usr/bin/env python3

# """"
# Exercise 1 from Mike Schatz tutorial on analyzing RNAseq data
# This python 3 script is to Visualize the provided e.coli differential gene expression
# from 100 genes across 10 timepoints
#
# http://schatz-lab.org/teaching/exercises/rnaseq/
#
# https://matplotlib.org/gallery/images_contours_and_fields/image_annotated_heatmap.html
#
# """"

import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA

df = pd.read_table(sys.argv[1], index_col=0)

col_names = df.columns.values.tolist()
row_names = df.index.values.tolist()

pca = PCA( n_components =2)
fit = pca.fit( df )
x = fit.transform( df )
# print( fit.explained_variance_ratio_ )
print( fit.components_ )
print( fit.components_.shape )
print( x )
print( x.shape )

fig, ax = plt.subplots(figsize=(6,6))

rng = np.random.RandomState(0)
colors = rng.rand(len(row_names))
ax.scatter( x[:, 0], x[:, 1], c=colors , cmap="RdYlGn" )
ax.set_title( "PCA of FPKM" )
ax.set_xlabel( "PCA1")
plt.xlim(-100,110)
ax.set_ylabel( "PCA2" )
ax.axhline(y=0, color='r', dashes=(1.0,1.0))

# some = int(len( row_names ))


# for i in range( some ):
#
#      ax.annotate( row_names[i], ( x[i,0], x[i,1] ))

plt.tight_layout()
fig.savefig( "pca_ecoliTimeExperiments.png" )
plt.close()



for i in range(len(row_names)):
    xaxis = col_names
    data = df.T.iloc[:,i]
    plt.plot(data)
    plt.ylim(0, 110)
    plt.xlabel('Experimental Time Point')
    plt.ylabel('FPKM?')
    plt.title("Time vs E.coli mRNA expression")
    continue

plt.plot(color=colors)
plt.savefig("PolyEColiLines.png")
plt.close()
    
#PCA
# poly lines 





