#!/usr/bin/env python3

#    compute average depth per exon from the reference gene list
#      depth_2nd.py depth.out > GOI_expressionchanges.out

import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
from scipy.cluster.hierarchy import dendrogram, linkage, leaves_list

df = pd.read_table(sys.argv[1], header=None)



timepoints = []


df1 = df.fillna(0)


timepoints = set(df1[0].values.tolist())
timepoints =list(timepoints)
        

# class TimecourseAnalysis():

#     def __init__(self)


# for i in range(len(timepoints)):
#     df_t1 = df1[(df1[0] == timepoints[i] )]
#     df_i = df_t1[[1, 3]]
#
#     df_i.columns = ["Gene", timepoints[i] + " Depth"]
#
#     continue
#


df_t1 = df1[(df1[0] == timepoints[0] )]
df_t1 = df_t1[[1, 3]]
df_t1.columns = ["Gene", timepoints[0]]
df_t1.index = df_t1["Gene"]
df_t1 = df_t1[timepoints[0]]
df_t1.columns = [timepoints[0]]

df_t2 = df1[(df1[0] == timepoints[1] )]
df_t2 = df_t2[[1, 3]]
df_t2.columns = ["Gene", timepoints[1]]
df_t2.index = df_t2["Gene"]
df_t2 = df_t2[timepoints[1]]
df_t2.columns = [timepoints[1]]

df_t3 = df1[(df1[0] == timepoints[2] )]
df_t3 = df_t3[[1, 3]]
df_t3.columns = ["Gene", timepoints[2]]
df_t3.index = df_t3["Gene"]
df_t3 = df_t3[timepoints[2]]
df_t3.columns = [timepoints[2]]

df_t4 = df1[(df1[0] == timepoints[3] )]
df_t4 = df_t4[[1, 3]]
df_t4.columns = ["Gene", timepoints[3]]
df_t4.index = df_t4["Gene"]
df_t4 = df_t4[timepoints[3]]
df_t4.columns = [timepoints[3]]

df_t5 = df1[(df1[0] == timepoints[4] )]
df_t5 = df_t5[[1, 3]]
df_t5.columns = ["Gene", timepoints[4]]
df_t5.index = df_t5["Gene"]
df_t5 = df_t5[timepoints[4]]
df_t5.columns = [timepoints[4]]

df_t6 = df1[(df1[0] == timepoints[5] )]
df_t6 = df_t6[[1, 3]]
df_t6.columns = ["Gene", timepoints[5]]
df_t6.index = df_t6["Gene"]
df_t6 = df_t6[timepoints[5]]
df_t6.columns = [timepoints[5]]

df_t7 = df1[(df1[0] == timepoints[6] )]
df_t7 = df_t7[[1, 3]]
df_t7.columns = ["Gene", timepoints[6]]
df_t7.index = df_t7["Gene"]
df_t7 = df_t7[timepoints[6]]
df_t7.columns = [timepoints[6]]

df_t8 = df1[(df1[0] == timepoints[7] )]
df_t8 = df_t8[[1, 3]]
df_t8.columns = ["Gene", timepoints[7]]
df_t8.index = df_t8["Gene"]
df_t8 = df_t8[timepoints[7]]
df_t8.columns = [timepoints[7]]

df_t9 = df1[(df1[0] == timepoints[8] )]
df_t9 = df_t9[[1,3]]
df_t9.columns = ["Gene", timepoints[8]]
df_t9.index = df_t9["Gene"]
df_t9 = df_t9[timepoints[8]]
df_t9.columns = [timepoints[8]]


df_t10 = df1[(df1[0] == timepoints[9] )]
df_t10 = df_t10[[1, 3]]
df_t10.columns = ["Gene", timepoints[9]]
df_t10.index = df_t10["Gene"]
df_t10 = df_t10[timepoints[9]]
df_t10.columns = [timepoints[9]]

# df_t1 = df_t6.loc["t6 Depth"]



frames = [df_t1, df_t2, df_t3, df_t4, df_t5, df_t6, df_t7, df_t8, df_t9, df_t10]
result = pd.concat(frames, axis=1)
# result = result.reset_index()


dfRegrex = result.reindex(columns=sorted(result.columns, key=lambda x: int(x[1:]) if x!='p' else 0))
result.columns = ["t1 Depth", "t2 Depth", "t3 Depth", "t4 Depth", "t5 Depth", "t6 Depth", "t7 Depth", "t8 Depth", "t9 Depth", "t10 Depth"]

print(dfRegrex)

# with open("dfRegrex.out", 'w') as file:
#     for item in dfRegrex.T:
#         file.write(str(item) + "\n")

col_names = dfRegrex.columns.values.tolist()
row_names = dfRegrex.index.values.tolist()

rng = np.random.RandomState(0)
colors = rng.rand(len(row_names))

d = "Invidividual GOI Depth graphs/"
if not os.path.exists(d):
    os.makedirs(d)

os.chdir(d)   


m = np.max(np.abs(dfRegrex))


for i in range(len(row_names)):
    xaxis = col_names
    data = dfRegrex.T.iloc[:,i]
    figure, ax = plt.subplots(figsize=(10, 6))
    plt.plot(data)
    # plt.ylim()
    plt.xlabel('Experimental Time Point')
    plt.ylabel('Depth')
    plt.title("Timepoint vs Depth of mRNA coverage for Gene:" + row_names[i] + " in E.coli")
    # plt.plot(color=colors)
    # plt.legend(row_names[i])
    plt.savefig("Parallel_lines_" + row_names[i] + ".png")
    plt.close()
    print("Completed Mean Individual Gene expression graph for Gene: " + row_names[i] + " ..")
    continue



d = ".."
os.chdir(d)  


xaxis = col_names
data = dfRegrex.T.iloc[:,:]
figure, ax = plt.subplots(figsize=(10, 6))
plt.plot(data)
plt.ylim(0, 110)
plt.xlabel('Experimental Time Point')
plt.ylabel('Depth?')
plt.title("Time vs Depth mRNA Coverage of Genes of Interest in E.coli")
plt.plot(color=colors)
# plt.legend(row_names)
plt.show()
plt.savefig("ParallelLines_EColi_GOI.png")
plt.close() 


# dfRegrexlog2 = np.log2(dfRegrex + 1)
# print(dfRegrexlog2)

# d = "Invidividual log2 Mean GOI Depth graphs/"
# if not os.path.exists(d):
#     os.makedirs(d)

# os.chdir(d)   


# m = np.max(np.abs(dfRegrex))


# for i in range(len(row_names)):
#     xaxis = col_names
#     data = dfRegrexlog2.T.iloc[:,i]
#     figure, ax = plt.subplots(figsize=(10, 6))
#     plt.plot(data)
#     # plt.ylim()
#     plt.xlabel('Experimental Time Point')
#     plt.ylabel('Log2 Depth')
#     plt.title("Timepoint vs Log2 Mean Depth of mRNA coverage for Gene: " + row_names[i] + " in E.coli")
#     plt.plot(color=colors)
#     plt.legend(row_names[i])
#     plt.savefig("Parallel_lines_" + row_names[i] + "_Mean_Log2.png")
#     plt.close()
#     print("Completed Log2 Mean Individual Gene expression graph for Gene:" + row_names[i] + " ..")
#     continue



# d = ".."
# os.chdir(d)  


# xaxis = col_names
# data = dfRegrexlog2.T.iloc[:,:]
# figure, ax = plt.subplots(figsize=(10, 6))
# plt.plot(data)
# # plt.ylim(0, 110)
# plt.xlabel('Experimental Time Point')
# plt.ylabel('Log2 Depth?')
# plt.title("Time vs Mean Log2 Depth mRNA Coverage of Genes of Interest in E.coli")
# plt.plot(color=colors)
# plt.legend(row_names)
# plt.show()
# plt.savefig("PolyEColiLines_GOI.png")
# plt.close() 

# exp_num = dfRegrex.columns.values.tolist()
# gene_names = dfRegrex.index.tolist()
# data = dfRegrex.values 
# matrix = result.values

# Z = linkage(data, "ward")
# ZT = linkage(data.T, "ward")


# idx_rows = leaves_list(Z)
# data = matrix[idx_rows, :]
# # idx_columns = leaves_list(ZT)
# # idx_columns = range(10)
# # data = data[:, idx_columns]

# # fig, ax = plt.subplots(figsize=(15, 3))
# # plt.plot(dfRegrex.T)
# # # ax.set_xlabel( "Gene", rotation=45)
# # ax.set_ylabel( "Depth" )
# # plt.ylim(0, 500)
# # plt.tight_layout()
# # fig.savefig( "DepthvsPositionEcoliTimeSeries_verts" + arg + ".png" )
# # plt.close()


# X = (data-np.average(data, axis=0))/np.std(data,axis=0)

# m = np.max(np.abs(X))

# fig, ax = plt.subplots(figsize=(15, 6))
# ax.set_title("Heatmap of Genes of Interest by Clustered Experiment Timepoints")
# im = ax.pcolor(X, cmap="cool", vmin = -m, vmax = m- 0.5)
# ax.grid(False)
# ax.set_xticks(np.arange(0.5, X.shape[1]+0.5),)
# ax.set_xticklabels(exp_num, rotation=50)
# ax.set_yticklabels(idx_rows, rotation=50)
# cbar = fig.colorbar(im, ax=ax)
# fig.subplots_adjust(left = 0.05, bottom = 0.15, right = 1.0, top = 0.95)
# fig.savefig("Depth_all_heatmap.png") # Save the image
# plt.close(fig) # Close the canvas



# fig, ax = plt.subplots(figsize=(15, 3))
# ax.set_title("Depth of All genes over time")
# plt.plot(dfRegrex.T)
# ax.set_xlabel( "Gene" )
# ax.set_ylabel( "Depth" )
# plt.ylim(0, 500)
# plt.tight_layout()
# fig.savefig( "Depth_ParallelCoordinates_All.png" )
# plt.close()

# increased = dfRegrex[(dfRegrex["t10"] > dfRegrex["t1"])] 
# list_increased = increased.index.tolist()
# # print(list_increased)
# fig, ax = plt.subplots(figsize=(15, 3))
# ax.set_title("Depth of Increased genes over time")
# plt.plot(increased.T)
# ax.set_xlabel( "Gene", rotation=45)
# ax.set_ylabel( "Depth" )
# plt.ylim(0, 500)
# plt.tight_layout()
# fig.savefig( "Depth_ParallelCoodinates_increased.png" )
# plt.close()


# decreased = dfRegrex[(dfRegrex["t1"] > dfRegrex["t10"])]
# list_decreased = decreased.index.tolist()
# # print(list_decreased)
# fig, ax = plt.subplots(figsize=(15, 3))
# ax.set_title("Depth of Decreased genes over time")
# plt.plot(decreased.T)
# ax.set_xlabel( "Gene" )
# ax.set_ylabel( "Depth" )
# plt.ylim(0, 500)
# plt.tight_layout()
# fig.savefig( "Depth_ParallelCoordinates_decreased.png" )
# plt.close()

# # neutral =  dfRegrex[(dfRegrex["t1"] < dfRegrex["t10"]) & (dfRegrex["t10"] > dfRegrex["t1"]) ]
# # list_neutral = neutral.index.tolist()
# # print(list_neutral)




# # filtered = result[(result < 400)]
# # sliced = data[(data.Position > start) & (data.Position < end)]
# # print(filtered)


# # plt.figure()
# # plt.plot(filtered)
# # plt.show()

# # pd.plotting.parallel_coordinates(result, result.index)
    
# # plt.show()

# # print(filtered)
# # data = filtered.values
# # filtered.dropna()


# data = increased.values
# Z = linkage(data, "ward")
# idx_rows = leaves_list(Z)
# data = data[idx_rows, :]
# X = (data-np.average(data, axis=0))/np.std(data,axis=0)

# m = np.max(np.abs(X))

# fig, ax = plt.subplots(figsize=(15, 6))
# ax.set_title("Heatmap of Increased Genes by Clustered Experiment Timepoints")
# im = ax.pcolor(X, cmap="viridis", vmin = -m, vmax = m- 0.5)
# ax.grid(False)
# ax.set_xticks(np.arange(0.5, X.shape[1]+0.5),)
# ax.set_xticklabels(exp_num, rotation=50)
# ax.set_yticklabels(idx_rows, rotation=50)
# cbar = fig.colorbar(im, ax=ax)
# fig.subplots_adjust(left = 0.05, bottom = 0.15, right = 1.0, top = 0.95)
# fig.savefig("Depth_increased_heatmap.png") # Save the image
# plt.close(fig) # Close the canvas

# data = decreased.values
# Z = linkage(data, "ward")
# idx_rows = leaves_list(Z)
# data = data[idx_rows, :]
# X = (data-np.average(data, axis=0))/np.std(data,axis=0)

# m = np.max(np.abs(X))

# fig, ax = plt.subplots(figsize=(15, 6))


# ax.set_title("Heatmap of Decreased Genes by Experiment Timepoints")
# im = ax.pcolor(X, cmap="viridis", vmin = -m, vmax = m- 0.5)
# ax.grid(False)
# ax.set_xticks(np.arange(0.5, X.shape[1]+0.5),)
# ax.set_xticklabels(exp_num, rotation=50)
# ax.set_yticklabels(idx_rows, rotation=50)
# cbar = fig.colorbar(im, ax=ax)
# fig.subplots_adjust(left = 0.05, bottom = 0.15, right = 1.0, top = 0.95)
# fig.savefig("Depth_decreased_heatmap.png") # Save the image
# plt.close(fig) # Close the canvas


# interesting_genes = {}
# gene_list = []
# dfRegrex_index = dfRegrex.index.tolist()

# for gene in dfRegrex_index:
#     if gene in list_increased:
#         interesting_genes[gene] = "increased"
#         continue
#     if gene in list_decreased:
#         interesting_genes[gene] = "decreased"
#         continue
#     else:
#         interesting_genes[gene] = "neutral"
#         continue
        


# for gene in dfRegrex_index:
#     if interesting_genes[gene] = "increased":
#         continue
#     if interesting_genes[gene] = "decreased":
#         continue
#     else:
#         interesting_genes[gene] = "neutral"
#         continue

# d = "GeneOntology/"
# if not os.path.exists(d):
#     os.makedirs(d)

# os.chdir(d)     

# with open("increased.out", 'w') as file:
#     for item in increasedNP:
#         file.write(str(item).upper() + "\n")

# with open("decreased.out", 'w') as file:
#     for item in decreasedNP:
#         file.write(str(item).upper() + "\n")

# with open("Meanslog2_inc.out", 'w') as file:
#     for item in Meanslog2_inc:
#         file.write(str(item).upper() + "\n")

# with open("Meanslog2_dec.out", 'w') as file:
#     for item in Meanslog2_dec:
#         file.write(str(item).upper() + "\n")

# with open("Means_inc.out", 'w') as file:
#     for item in Means_inc:
#         file.write(str(item).upper() + "\n")

# with open("Means_dec.out", 'w') as file:
#     for item in Means_dec:
#         file.write(str(item).upper() + "\n")
# for gene in interesting_genes:
#     print("Gene: " + gene + "\tExpression: " + interesting_genes[gene] )

# print(interesting_genes, interesting_genes[gene])















