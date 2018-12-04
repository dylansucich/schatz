#!/Users/dylansucich/miniconda3/bin/python3

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


import os
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from scipy.cluster.hierarchy import dendrogram, linkage, leaves_list

# df = pd.read_table(sys.argv[1], index_col=0)

# col_names = df.columns.values.tolist()
# row_names = df.index.values.tolist()

class RNAseq:

	self = pd.read_table(sys.argv[1], index_col=0)
	df = pd.read_table(sys.argv[1], index_col=0)
	timepoints = self.columns.values.tolist()
	genes = self.index.values.tolist()
	

	def __init__(self):
		
		self = pd.read_table(sys.argv[1], index_col=0)
		df = pd.read_table(sys.argv[1], index_col=0)
		timepoints = self.columns.values.tolist()
		genes = self.index.values.tolist()

	def pca(self):

		pca = PCA( n_components =2)
		fit = pca.fit(self)
		x = fit.transform(self)
		# print( fit.explained_variance_ratio_ )
		print("PCA Components: ", fit.components_ )
		print("PCA Shape: ", fit.components_.shape )
		print("PCA Fit?: ", x )
		print("PCA Fit Shape: ", x.shape )
		
		fig, ax = plt.subplots(figsize=(6,6))
		ax.scatter( x[:, 0], x[:, 1], c=colors , cmap="RdYlGn" )
		ax.set_title( "PCA of FPKM" )
		ax.set_xlabel( "PCA1")
		# plt.xlim(-100,110)
		ax.set_ylabel( "PCA2" )
		ax.axhline(y=0, color='r', dashes=(1.0,1.0))
		plt.tight_layout()
		fig.savefig( "pca_timecourse_rnaseq.png" )
		plt.close()

	def pt(self):

		for i in range(len(genes)):
		    xaxis = timepoints
		    data = self.T.iloc[:,i]
		    plt.plot(data)
		    plt.ylim(0, 110)
		    continue

		plt.xlabel('Experimental Time Point')
		plt.ylabel('FPKM')
		title = input("Please input Parallel lines plot title: ")
		plt.title(title)
		plt.plot(color=colors)
		plt.savefig("plot_parallel_lines.png")
		plt.close()

	def heatmap(self):
		print(genes)
		print(timepoints)
		matrix = self.values
		Z = linkage(matrix, "ward")
		ZT = linkage(matrix.T, "ward")

		fig, ax = plt.subplots()

		title = input("Please input Clustered Dendrogram plot title: ")
		plt.title(title)

		plt.xlabel("sample")
		plt.ylabel("distance")
		dendrogram( ZT, labels= timepoints, show_leaf_counts=False, leaf_rotation=90., leaf_font_size=12., show_contracted=True)

		fig.savefig("dendrogram.png")
		plt.close(fig)

		idx_rows = leaves_list(Z)
		data = matrix[idx_rows, :]
		# idx_columns = leaves_list(ZT)
		idx_columns = range(10)
		data = data[:, idx_columns]

		X = (data-np.average(data,axis=0))/np.std(data,axis=0)

		m = np.max(np.abs(X))
		y = np.arange(0,100)

		fig, ax = plt.subplots()
		
		title = input("Please input Heatmap title: ")

		ax.set_title(title)
		im = ax.pcolor(X, cmap="viridis", vmin = -m, vmax = m- 0.5)
		ax.grid(False)
		ax.set_xticks(np.arange(0.5, X.shape[1]+0.5),)
		ax.set_xticklabels(timepoints, rotation=50)
		ax.set_yticks(np.arange(0.5, len(genes), 5), genes)
		ax.set_yticklabels(genes)
		cbar = fig.colorbar(im, ax=ax)
		fig.subplots_adjust(left = 0.05, bottom = 0.15, right = 1.0, top = 0.95)
		fig.savefig("Heatmap_clustered.png" ) # Save the image
		plt.close(fig) # Close the canvas
		print("Clustered Heatmap for Complete...")


d = "output/"
if not os.path.exists(d):
    os.makedirs(d)
os.chdir(d) 

print(RNAseq.df)
df = RNAseq.df
timepoints = RNAseq.timepoints
genes = RNAseq.genes
rng = np.random.RandomState(len(genes))
colors = rng.rand(len(genes))
pca = RNAseq.pca(df)
plot = RNAseq.pt(df)
heatmap = RNAseq.heatmap(df)
os.chdir("..") 















