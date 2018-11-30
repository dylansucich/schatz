#!/usr/bin/env python3

# Extract Genes that are differentially expressed in Medial subnucleus 
# of Central Amygdala SOM+ gabaergic neurons

# First run plot log2(p-value) on y and log2FoldChange on x
# Highlight above same threshold seen in paper for y and on x
 
# See if extracting by p-value or adjusted p-value results in different 
# lists

import os
import sys
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("agg")

import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from scipy.cluster.hierarchy import dendrogram, linkage, leaves_list

df = pd.read_table(sys.argv[1], index_col=0)

row_names = df.index.tolist()
col_names = df.columns.tolist()

print(col_names)
# print(len(row_names))

total_genes1 = len(row_names)


pval1 = df.loc[:,"pvalue"]
pval2 = pval1.dropna()

droppedpval = len(pval1) - len(pval2)
# print(len(pval1))
# print(len(pval2))

covered = len(pval1)

pval3 = df.loc[:,"padj"]
pval4 = pval3.dropna()

droppedpadj = len(pval3) - len(pval4)
# print(len(pval3))
# print(len(pval4))

sizes = [len(pval2), droppedpval]
sizes2 = [len(pval4), droppedpadj]

pvalpercent_covered = (len(pval2) / total_genes1) * 100  
padjpercent_covered = (len(pval4) / total_genes1) * 100  

labels = ["Covered by pval, Total"]

print("\nTotal number of Mus Musculus Coding genes: " + str(len(row_names)))
print("Genes detected through Ribo-tag RNA-seq: " + str(len(pval2)))
print("Percent: " + str(pvalpercent_covered))
print("Genes detected through Ribo-tag RNA-seq with adjusted pval: " + str(len(pval4)))
print("Percent: " + str(padjpercent_covered) + "\n")


df_pval = pval2
df_padj = pval4

df_log2_pval = np.log2(pval2)
df_log2_padj = np.log2(pval4)
df_log10_pval = np.log(pval2)
df_log10_padj = np.log(pval4)

df_base = df.loc[:, "baseMean"]
df_log2 = df.loc[:,"log2FoldChange"]
df_log2 = df_log2.dropna()
df_lfcSE = df.loc[:, "lfcSE"]
df_lfcSE = df_lfcSE.dropna()
df_stat = df.loc[:, "stat"]
df_stat = df_stat.dropna()

mask_pval = (0.05 < df_pval.iloc[:])
mask_padj = (0.05 < df_padj.iloc[:])
mask_log2_pval = (0.05 < df_log2_pval.iloc[:])
mask_log2_padj = (0.05 < df_log2_padj.iloc[:])
mask_log10_pval = (0.05 < df_log10_pval.iloc[:])
mask_log10_padj = (0.05 < df_log10_padj.iloc[:])

df_log2_na1 = df.loc[:,["log2FoldChange", "pvalue"]]
print("df_log2_filtered:", len(df_log2_na1))
df_log2_na1 = df_log2_na1.dropna()
print("df_log2_filtered_pval_na: ", len(df_log2_na1))
print("df_log2_pval: ", len(df_log2_pval))

df_log2_na2 = df.loc[:,["log2FoldChange", "padj"]]
print("\ndf_log2_filtered:", len(df_log2_na2))
df_log2_na2 = df_log2_na2.dropna()
print("df_log2_filtered_padj_na:", len(df_log2_na2))
print("df_log2_padj: ", len(df_log2_padj))

df_mpval = df_log2_na1.isin(df_pval)
print("\ndf_mpval: ", len(df_mpval))
print("df_pval: ", len(df_pval))
df_mpadj= df_log2_na2.isin(df_padj)
print("df_mpadj: ", len(df_mpadj))
print("df_padj: ", len(df_padj))

d = "output/"
if not os.path.exists(d):
    os.makedirs(d)

os.chdir(d)   

plt.pie(sizes)
plt.title("Coverage by pval of all mouse genes")
plt.savefig("Coveragebypval.png")
plt.close()

plt.pie(sizes2)
plt.title("Coverage by padj of all mouse genes")
plt.savefig("Coveragebypadj.png")
plt.close()

plt.hist(df_base, bins = 20, color = "r")
plt.title("Distribution of Mean Gene Expression")
plt.savefig("histogram_basemean.png")
plt.close()

plt.hist(df_log2, bins = 100, color = "r")
plt.title("Distribution of log2 Fold change in Gene Expression")
plt.savefig("histogram_log2FoldChange.png")
plt.close()

plt.hist(df_lfcSE, bins = 100, color = "r")
plt.title("Distribution of 'lfcSE' in Gene Expression")
plt.savefig("histogram_lfcSE.png")
plt.close()

plt.hist(df_log2, bins = 100, color = "r")
plt.title("Distribution of 'Stat' in Gene Expression")
plt.savefig("histogram_stat.png")
plt.close()

plt.hist(df_pval, bins = 100, color = "r", alpha = 0.5)
plt.hist(df_pval[mask_pval], bins = 100, color = "b")
plt.title("Distribution of Genes P-values")
plt.savefig("histogram_pval.png")
plt.close()

plt.hist(df_padj, bins = 100, color = "r", alpha = 0.5)
plt.hist(df_padj[mask_padj], bins = 100, color = "b")
plt.title("Distribution of Genes Adjusted P-values")
plt.savefig("histogram_padj.png")
plt.close()

plt.hist(df_log2_pval, bins = 100, color = "r", alpha = 0.5)
plt.hist(df_log2_pval[mask_log2_pval], bins = 100, color = "b")
plt.title("Distribution of Genes log2 P-values")
plt.savefig("histogram_log2_log2pval.png")
plt.close()

plt.hist(df_log2_padj, bins = 100, color = "r", alpha = 0.5)
plt.hist(df_log2_padj[mask_log2_padj], bins = 100, color = "b")
plt.title("Distribution of Genes Adjusted log2 P-values")
plt.savefig("histogram_log2_padj.png")
plt.close()

plt.hist(df_log10_pval, bins = 100, color = "r", alpha = 0.5)
plt.hist(df_log10_pval[mask_log10_pval], bins = 100, color = "b")
plt.title("Distribution of Genes log10 P-values")
plt.savefig("histogram_log10_pval.png")
plt.close()

plt.hist(df_log10_padj, bins = 100, color = "r", alpha = 0.5)
plt.hist(df_log10_padj[mask_log10_padj], bins = 100, color = "b")
plt.title("Distribution of Genes Adjusted log10 P-values")
plt.savefig("histogram_log10_padj.png")
plt.close()

plt.scatter(df_mpval.loc[:,"log2FoldChange"], df_pval, color = "r", alpha = 0.5)
plt.title("Scatter test pval")
plt.savefig("scatter_test_pval.png")
plt.close()

plt.scatter(df_mpadj.loc[:,"log2FoldChange"], df_padj, color = "r", alpha = 0.5)
plt.title("Scatter test padj")
plt.savefig("scatter_test_padj.png")
plt.close()

plt.scatter(df_mpval.loc[:,"log2FoldChange"], df_pval, color = "r", alpha = 0.5)
plt.title("Scatter test i pval")
plt.gca().invert_yaxis()
plt.savefig("scatter_test_i_pval.png")
plt.close()

plt.scatter(df_mpadj.loc[:,"log2FoldChange"], df_padj, color = "r", alpha = 0.5)
plt.title("Scatter test i padj")
plt.gca().invert_yaxis()
plt.savefig("scatter_test_i_padj.png")
plt.close()

plt.scatter(df_mpval.loc[:,"log2FoldChange"], df_log2_pval, color = "r", alpha = 0.5)
plt.title("Scatter test log2 i pval")
plt.gca().invert_yaxis()
plt.savefig("scatter_test_log2_i_pval.png")
plt.close()

plt.scatter(df_mpadj.loc[:,"log2FoldChange"], df_log2_padj, color = "r", alpha = 0.5)
plt.title("Scatter test log2 i padj")
plt.gca().invert_yaxis()
plt.savefig("scatter_test_log2_i_padj.png")
plt.close()

plt.scatter(df_mpval.loc[:,"log2FoldChange"], df_log10_pval, color = "r", alpha = 0.5)
plt.title("Scatter test log10 pval")
plt.gca().invert_yaxis()
plt.savefig("scatter_test_log10_i_pval.png")
plt.close()

plt.scatter(df_mpadj.loc[:,"log2FoldChange"], df_log10_padj, color = "r", alpha = 0.5)
plt.title("Scatter test log10 i padj")
plt.gca().invert_yaxis()
plt.savefig("scatter_test_log10_i_padj.png")
plt.close()


os.chdir("..")







# for df.loc[:,""]