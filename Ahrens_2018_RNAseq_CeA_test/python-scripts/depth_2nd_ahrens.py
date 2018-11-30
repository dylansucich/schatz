#!/usr/bin/env python3

#    compute average depth per exon from the reference gene list
#      depth_2nd.py depth.out > GOI_expressionchanges.out

import os
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from scipy.cluster.hierarchy import dendrogram, linkage, leaves_list

df = pd.read_table(sys.argv[1], index_col=0)

print(".......................................")
print("\nSample data source:" + sys.argv[1])
print("\nSource paper:\n\n\tAhrens et al. 2018\n\tJ.Neuroscience\n\tCentral Amygdala circuit modulates anxiety\n")
samples_list = []
samples_list = df.columns.tolist()
samples_list = samples_list[0:]
samples_list_count = len(samples_list)

IP_total_count = 0
IP_M_count = 0 
IP_F_count = 0

Control_total_count = 0
Control_M_count = 0
Control_F_count = 0

GOI_list = []

sample_types = df.columns.values.tolist()
gene_names = df.index.tolist()


print(".......................................\n" + "\nList of Samples: \n" + str(samples_list) + "\n")
print("\n\tNumber of Genes Covered: " + str(len(gene_names)) +"\n")
print("\tNumber of Samples: " + str(samples_list_count))

for i in sample_types:
    if i.startswith("IP"):
        IP_total_count += 1 
        
        if i.startswith("IP_Male"):
            IP_M_count += 1
            
        else:
            IP_F_count +=1
        continue
        
    if i.startswith("Input"):
        Control_total_count += 1
        
        if i.startswith("Input_Male"):
            Control_M_count += 1
            
        else:
            Control_F_count +=1
        continue
        
    continue
    
print("\n\tIP Samples: " + str(IP_total_count) + "\n\tMale,Female: " + str(IP_M_count) + ", " + str(IP_F_count) + "\n\n\tInput Samples: " + str(Control_total_count) + "\n\tMale, Female: " + str(Control_M_count) + ", " + str(Control_F_count) + "\n\n")
print(".......................................")
print("IP Samples:\n" + str(sample_types[0:10]) + "\n")
print("INPUT Samples:\n" + str(sample_types[10:]) + "\n.......................................")

pd.set_option("mode.chained_assignment", None)
# print("Ignore the Error below about '.loc', it is alright\n")


GOI_list = ["Slc17a1", "Slc17a2", "Slc17a3", "Slc17a4", "Slc17a5", "Slc17a6", "Slc17a7", "Slc17a8", "Slc17a9", "Grin1", "Grin2b", "Gls", "Gls2", "Gad1", "Gad2", "Slc32a1", "Hspa8"]


# GOI_list = ["Sst", "Penk", "Pdyn", "Tac2", "Nts", "Vgf", "Tac1", "Nxph1", "Oprk1", "Cck", "Pnoc", "Crh", "Pomc", "Ucn3", "Uts2", "Oxt", "Ucn", "Avp", "Ucn2", "Nms", "Nmu"]
print("Neuropeptides of Interest: ")
print(str(GOI_list[0:5]) +"\n" + str(GOI_list[5:10]) +"\n" + str(GOI_list[10:15]) +"\n" + str(GOI_list[15:21]))

print(".......................................")

print("Differential expression matrix for Input and IP constructing...")

matrix1 = df.iloc[:,0:10]
matrix1.loc[:, "All IP mean"] = (matrix1.mean(axis=1) + 1)
dfALL_IP = pd.DataFrame(matrix1.loc[:, "All IP mean"])

matrix2 = df.iloc[:,0:10]
matrix2.loc[:, "All IP log2 mean"] = np.log2((matrix2.mean(axis=1) +1))
dfALL_IP_log2 = pd.DataFrame(matrix2.loc[:, "All IP log2 mean"])

matrix3 = df.iloc[:,10:]
matrix3.loc[:, "All INPUT mean"] = (matrix3.mean(axis=1) + 1)
dfALL_INPUT = pd.DataFrame(matrix3.loc[:, "All INPUT mean"])

matrix4 = df.iloc[:,10:]
matrix4.loc[:, "All INPUT log2 mean"] = np.log2((matrix4.mean(axis=1)) + 1 )
dfALL_INPUT_log2 = pd.DataFrame(matrix4.loc[:, "All INPUT log2 mean"])

ALL_samples = dfALL_IP.join(dfALL_INPUT)
ALL_samples_log2 = dfALL_IP_log2.join(dfALL_INPUT_log2)
# print(ALL_samples)
# print(ALL_samples_log2)

ALL_samples_genes = ALL_samples.index.tolist()
ALL_samples_groups = ALL_samples.columns.tolist()

ALL_samples_matrix = ALL_samples.values

dfRegrexMeanslog2 = ALL_samples_log2
dfRegrexMeanslog2_columns = dfRegrexMeanslog2.columns.values.tolist()
dfRegrexMeanslog2_gene_names = dfRegrexMeanslog2.index.tolist()

dfRegrexMeans = ALL_samples
dfRegrexMeans_columns = dfRegrexMeans.columns.values.tolist()
dfRegrexMeans_gene_names = dfRegrexMeans.index.tolist()

increasedMeanslog2 = dfRegrexMeanslog2[(dfRegrexMeanslog2["All IP log2 mean"] > dfRegrexMeanslog2["All INPUT log2 mean"])]
list_increasedMeanslog2 = increasedMeanslog2.index.tolist()
# print(list_increasedMeanslog2)

increasedMeans = dfRegrexMeans[(dfRegrexMeans["All IP mean"] > dfRegrexMeans["All INPUT mean"])]
list_increasedMeans = increasedMeans.index.tolist()
# print(list_increasedMeans)

decreasedMeanslog2 = dfRegrexMeanslog2[(dfRegrexMeanslog2["All IP log2 mean"] < dfRegrexMeanslog2["All INPUT log2 mean"])]
list_decreasedMeanslog2 = decreasedMeanslog2.index.tolist()
# print(list_decreasedMeanslog2)

decreasedMeans = dfRegrexMeans[(dfRegrexMeans["All IP mean"] < dfRegrexMeans["All INPUT mean"])]
list_decreasedMeans = decreasedMeans.index.tolist()
# print(list_decreasedMeans)


Meanslog2_interesting_genes = {}
Meanslog2_gene_list = []
dfRegrexMeanslog2_index = dfRegrexMeanslog2.index.tolist()
Meanslog2_inc = [] 
Meanslog2_dec = []
Meanslog2_neu = []
Meanslog2_inc_count = 0 
Meanslog2_dec_count = 0 
Meanslog2_neu_count = 0 

for gene in dfRegrexMeanslog2_index:
    if gene in list_increasedMeanslog2:
        Meanslog2_interesting_genes[gene] = "increased"
        Meanslog2_inc.append(gene)
        Meanslog2_inc_count += 1
        continue
    if gene in list_decreasedMeanslog2:
        Meanslog2_interesting_genes[gene] = "decreased"
        Meanslog2_dec.append(gene)
        Meanslog2_dec_count += 1
        continue
    else:
        Meanslog2_interesting_genes[gene] = "neutral"
        Meanslog2_neu.append(gene)
        Meanslog2_neu_count += 1
        continue
        
        
Means_interesting_genes = {}
Means_gene_list = []
dfRegrexMeans_index = dfRegrexMeans.index.tolist()
Means_inc = [] 
Means_dec = []
Means_neu = []
Means_inc_count = 0 
Means_dec_count = 0 
Means_neu_count = 0 
Means_total = 0 

for gene in dfRegrexMeans_index:
    if gene in list_increasedMeans:
        Means_interesting_genes[gene] = "increased"
        Means_inc.append(gene) 
        Means_inc_count += 1 
        continue
    if gene in list_decreasedMeans:
        Means_interesting_genes[gene] = "decreased"
        Means_dec.append(gene)
        Means_dec_count += 1
        continue
    else:
        Means_interesting_genes[gene] = "neutral"
        Means_neu.append(gene)
        Means_neu_count += 1
        continue
        
increasedNP = []
incNPcount = 0 
decreasedNP = []  
decNPcount = 0   

for gene in GOI_list:
    if gene in  Meanslog2_interesting_genes:
        # print("log2means_Neuropeptide: " + gene + "\tIP vs Input: " + Meanslog2_interesting_genes[gene])
        if Meanslog2_interesting_genes[gene] == "increased":
            increasedNP.append(gene)
            incNPcount += 1
            continue
        if Meanslog2_interesting_genes[gene] == "decreased":
            decreasedNP.append(gene)
            decNPcount += 1
            continue
        continue
    if gene in  Means_interesting_genes:
        # print("means_Neuropeptide: " + gene + "\tIP vs Input: " + Means_interesting_genes[gene])
        if Means_interesting_genes[gene] == "increased":
            increasedNP.append(gene)
            continue
        if Means_interesting_genes[gene] == "decreased":
            decreasedNP.append(gene)
            continue
        continue



# fileincIPvsInput = increasedNP.out
# fileincIPvsInput.writelines( "%s\n" % item for item in increasedNP )
# filedecIPvsInput = decreasedNP.out
# filedecIPvsInput.writelines( "%s\n" % item for item in decreasedNP )


# dir_path = os.path.join(self.feed, self.address)
#
# print(dir_path) # will return 'feed/address'
# os.makedirs(dir_path, exist_ok=True)  # create directory [current_path]/feed/address
# output = open(os.path.join(dir_path, file_name), 'wb')

print("\nNeuropeptides enriched in IP vs. Input...")
print("Increased: "+ str(incNPcount) + "\n" + str(increasedNP) + "\nDecreased: " + str(decNPcount) + "\n" + str(decreasedNP))
   
print("\nMeans Increased: " + str(Means_inc_count) + "\nMeans Decreased: " + str(Means_dec_count) + "\nMeans Neutral: " + str(Means_neu_count)) 
print("Means Total: " + str(Means_inc_count + Means_dec_count + Means_neu_count))
print("\nLog2 Means Increased: " + str(Meanslog2_inc_count) + "\nLog2 Means Decreased: " + str(Meanslog2_dec_count) + "\nLog2 Means Neutral: " + str(Meanslog2_neu_count)) 
print("Log2 Means Total: " + str(Meanslog2_inc_count + Meanslog2_dec_count + Meanslog2_neu_count))

print("\nComplete...")

print(".......................................")



ALL_samples = dfALL_IP.join(dfALL_INPUT)
ALL_samples_log2 = dfALL_IP_log2.join(dfALL_INPUT_log2)

dfRegrexALL = ALL_samples
dfRegrexALLlog2 = ALL_samples_log2

# Parallel Lines !!!!!!!
col_names = dfRegrexALL.columns.values.tolist()
row_names = dfRegrexALL.index.values.tolist()

rng = np.random.RandomState(0)
colors = rng.rand(len(row_names))

# d = "Invidividual NP Depth graphs for Mean Depth/"
# if not os.path.exists(d):
#     os.makedirs(d)

# os.chdir(d)   


# m = np.max(np.abs(dfRegrex))


# for i in range(len(row_names)):
#     xaxis = col_names
#     data = dfRegrexALL.T.loc[:,row_names[i]]
#     figure, ax = plt.subplots(figsize=(10, 6))
#     plt.plot(data)
#     # plt.ylim()
#     plt.xlabel('Average Means for IP vs. Input')
#     plt.ylabel('Depth?')
#     plt.title("Differential Mean Depth of mRNA coverage for Gene:" + row_names[i] + " in CeA")
#     plt.plot(color=colors)
#     plt.savefig("Parallel_lines_" + GOI_list[i] + ".png")
#     plt.close()
#     print("Completed Individual Gene expression graph for Gene:" + row_names[i] + " ..")
#     continue



# d = ".."
# os.chdir(d)  


# xaxis = col_names
# data = dfRegrexALL.T.loc[:, GOI_list[:]]

# figure, ax = plt.subplots(figsize=(10, 6))
# plt.plot(data)
# # plt.ylim(0, 110)
# plt.xlabel('Average Means for IP vs. Input')
# plt.ylabel('Depth?')
# plt.title("Differential Depth of mRNA coverage for all Genes in Sst+ vs. Sst- GABA neurons in the CeA")
# plt.plot(color=colors)
# plt.legend(GOI_list)
# plt.savefig("PolylinesIPvsInput.png")
# plt.close() 
# print("Completed Gene expression graph Mean Polylines..\n")

# GABA = [""]

d = "Invidividual GABA genes VIP-CCK vs. Sst Depth graphs for Log2 Mean Depth/"
if not os.path.exists(d):
    os.makedirs(d)

os.chdir(d)   


# m = np.max(np.abs(dfRegrex))


for i in range(len(GOI_list)):
    xaxis = col_names
    data = dfRegrexALLlog2.T.loc[:, GOI_list[i]]
    figure, ax = plt.subplots(figsize=(10, 6))
    plt.plot(data)
    # plt.ylim()
    plt.xlabel('Average Means for IP vs. Input')
    plt.ylabel('Log2 Depth?')
    plt.title("Differential log2 Mean Depth of mRNA coverage for Gene:" +  GOI_list[i] + " in CeA")
    plt.plot(color=colors)
    plt.savefig("Parallel_lines_" + GOI_list[i] + "_Means_log2.png")
    plt.close()
    print("Completed Individual Gene expression graph for Gene:" + GOI_list[i] + " ..")
    continue



d = ".."
os.chdir(d)  


xaxis = GOI_list
data = dfRegrexALLlog2.T.loc[:,GOI_list]
figure, ax = plt.subplots(figsize=(10, 6))
plt.plot(data)
# plt.ylim(0, 110)
plt.xlabel('Average Means for IP vs. Input')
plt.ylabel('Depth?')
plt.title("Differential Depth of mRNA coverage for Glut/GABA markers in Sst+ vs. Sst- GABA neurons in the CeA")
plt.plot(color=colors)
# plt.legend(row_names)
plt.show()
plt.savefig("PolylinesGlut-GABA_IPvsInput_Means_log2.png")
plt.close() 

print("Completed Gene expression graph log2 Mean Polylines...\n")


# print(ALL_samples)
# print(ALL_samples_log2)

# ALL_samples_genes = ALL_samples.index.tolist()
# ALL_samples_groups = ALL_samples.columns.tolist()

# ALL_samples_matrix = ALL_samples.values

# print(ALL_samples_matrix)



# for i in ALL_samples:
#     fig, ax = plt.subplots()
#     plt.plot(ALL_samples_genes, ALL_samples )
#     ax.set_title("Plot Clustered X and Y for " + i)
#     ax.set_xlabel( ALL_samples_genes, rotation=45)
#     ax.set_ylabel( ALL_samples_groups )
#     plt.ylim(0, 500)
#     plt.tight_layout()
#     fig.savefig("BasicPlot" + i + ".png" )
#     plt.close(fig)
#     print("Basic Plot for " + i + " Complete...")
#     continue


    
# All_samples = pd.merge(dfALL_IP, dfALL_INPUT, on= 0 , how= "inner")



# if str(sys.argv[2]) or str(sys.argv[3]) == "ALL":
#     print("\n\n\nIndividual Sample Analysis= True")
#
#     # For All samples
#     for i in range(len(sample_types)):
#         print(".......................................\n" + sample_types[i])
#         matrix = df.iloc[:,i]
#
#         idx_rows = matrix.index
#         data = matrix[idx_rows, :]
#         idx_columns = leaves_list(ZT)
#         # print("Clustering Complete...")
#         # # idx_columns = range(10)
#         data = matrix.values
#
#         X = (data-np.average(data, axis=0))/np.std(data,axis=0)
#
#         m = np.max(np.abs(X))
#         print("Max of normalized data: " + str(m))
#         Xmin = np.min(np.abs(X))
#         print("Min of normalized data: " + str(Xmin) +"\n")
#
#         fig, ax = plt.subplots()
#         plt.plot(X)
#         ax.set_title("Plot Clustered X and Y for " + str(sample_types[i]))
#         ax.set_xlabel( "?", rotation=45)
#         ax.set_ylabel( "Expression" )
#         plt.ylim(0, 500)
#         plt.tight_layout()
#         fig.savefig("BasicPlot" + str(sample_types[i]) + ".png" )
#         plt.close(fig)
#         print("Basic Plot for " + str(sample_types[i]) + " Complete...")
#
#
#         fig, ax = plt.subplots()
#         plt.plot(X.T)
#         ax.set_title("Plot Transposed Clustered X and Y for " + str(sample_types[i]))
#         ax.set_xlabel( "?", rotation=45)
#         ax.set_ylabel( "Expression" )
#         plt.ylim(0, 500)
#         plt.tight_layout()
#         fig.savefig("BasicPlot_Transposed"+ str(sample_types[i]) + ".png" )
#         plt.close(fig)
#         print("Transposed Basic Plot for " + str(sample_types[i]) + " Complete...")
#         print(".......................................\n")
#
#
#         fig, ax = plt.subplots()
#         ax.set_title("Heatmap of Clustered X and Y for " + str(sample_types[i]))
#         im = ax.pcolor(X, cmap="viridis", vmin = -m, vmax = m- 0.5)
#         ax.grid(False)
#         ax.set_xticks(np.arange(0.5, X.shape[1]+0.5),)
#         ax.set_xticklabels(idx_columns, rotation=50)
#         ax.set_yticklabels(idx_rows, rotation=50)
#         cbar = fig.colorbar(im, ax=ax)
#         fig.subplots_adjust(left = 0.05, bottom = 0.15, right = 1.0, top = 0.95)
#         fig.savefig("Clustered_Heatmap_" + str(sample_types[i]) +".png" ) # Save the image
#         plt.close(fig) # Close the canvas
#         print("Clustered Heatmap for " + str(sample_types[i]) + " Complete...")
#
#
#         fig, ax = plt.subplots()
#         ax.set_title("Transposed Heatmap of Clustered X and Y for " + str(sample_types[i]))
#         im = ax.pcolor(X.T, cmap="viridis", vmin = -m, vmax = m- 0.5)
#         ax.grid(False)
#         ax.set_xticks(np.arange(0.5, X.shape[1]+0.5),)
#         ax.set_xticklabels(idx_columns, rotation=50)
#         ax.set_yticklabels(idx_rows, rotation=50)
#         cbar = fig.colorbar(im, ax=ax)
#         fig.subplots_adjust(left = 0.05, bottom = 0.15, right = 1.0, top = 0.95)
#         fig.savefig("Clustered_Heatmap_Transposed_"+ str(sample_types[i]) + ".png" ) 
# Save the image
#         plt.close(fig) # Close the canvas
        # print("Transposed Clustered Heatmap for " + str(sample_types[i]) + " Complete..." + "\n\n")
#         print(".......................................\n")
#
#     else:
#
#         print(".......................................\n" +
#         "No Individiual Sample Analysis.......................................\n")


#######################
# for all IP
# matrix = ALL_samples.values

# Z = linkage(ALL_samples, "ward")
# ZT = linkage(ALL_samples.T, "ward")

# fig, ax = plt.subplots()
# plt.title("Dendrogram of Input vs IP" )
# plt.xlabel("sample")
# plt.ylabel("distance")
# dendrogram(
#     ZT,
#     show_leaf_counts=False,
#     leaf_rotation=90.,
#     leaf_font_size=12.,
#     show_contracted=True)

# fig.savefig("ClusteredDendrogram_Transposed_InputvsIP.png")
# plt.close(fig)
# print("Dendrogram for Input vs IP Complete...")

# idx_rows = leaves_list(Z)
# data = matrix[idx_rows, :]
# idx_columns = leaves_list(ZT)
# print("Clustering Complete...")
# print(".......................................")
# # idx_columns = range(10)
# data = data[:, idx_columns]

# # X = (data-np.average(data, axis=0))/np.std(data,axis=0)
# X = (data-np.average(data, axis=0))/np.std(data,axis=0)
# m = np.max(np.abs(X))
# print("Max of normalized data: " + str(m))
# Xmin = np.min(np.abs(X))
# print("Min of normalized data: " + str(Xmin))

# fig, ax = plt.subplots(figsize=(20, 3))
# plt.plot(X)
# ax.set_title("Plot Clustered X and Y for Input vs IP")
# ax.set_xlabel(idx_columns, rotation=45)
# ax.set_ylabel( idx_rows)
# plt.ylim(Xmin, 50)
# plt.tight_layout()
# fig.savefig("BasicPlotInputvsIP.png" )
# plt.close(fig)
# print("Basic Plot for Input vs IP Complete...")
# print(".......................................")

##########################
# fig, ax = plt.subplots(figsize=(20,3))
# plt.plot(X.T)
# ax.set_title("Plot Transposed Clustered X and Y for " + "InputvsIP")
# ax.set_xlabel( idx_rows, rotation=45)
# ax.set_ylabel( idx_columns )
# plt.ylim(Xmin, m)
# plt.tight_layout()
# fig.savefig("BasicPlot_Transposed"+ "InputvsIP.png" )
# plt.close(fig)
# print("Transposed Basic Plot for " + "InputvsIPComplete...")

# fig, ax = plt.subplots(figsize=(6, 8))
# ax.set_title("Log Transformed Heatmap of Clustered X and Y for " + "Input vs IP")
# im = ax.pcolor(X, cmap="cool", vmin = -m, vmax = m- 0.5)
# ax.grid(False)
# ax.set_xticks(np.arange(0.5, X.shape[1]+0.5),)
# ax.set_xticklabels(idx_columns, rotation=50)
# ax.set_yticklabels(idx_rows, rotation=50)
# cbar = fig.colorbar(im, ax=ax)
# fig.subplots_adjust(left = 0.05, bottom = 0.15, right = 1.0, top = 0.95)
# # plt.show()
# fig.savefig("TransformedClustered_Heatmap_" + "InputvsIP.png" ) # Save the image
# plt.close(fig) # Close the canvas
# print("Clustered Heatmap for " + "InputvsIP Complete...")
#
#
# trans = X.T
# fig, ax = plt.subplots(figsize=(6, 8))
# ax.set_title("Transposed Heatmap of Clustered X and Y for " + "InputvsIP")
# im = ax.pcolor(X.T, cmap="cool", vmin = -m, vmax = m- 0.5)
# ax.grid(False)
# ax.set_xticks(np.arange(0.5, trans.shape[1]+0.5),)
# ax.set_xticklabels(idx_rows, rotation=50)
# ax.set_yticklabels(idx_columns, rotation=50)
# cbar = fig.colorbar(im, ax=ax)
# fig.subplots_adjust(left = 0.05, bottom = 0.15, right = 1.0, top = 0.95)
# fig.savefig("Clustered_Heatmap_Transposed_"+ "InputvsIP.png" ) # Save the image
# plt.close(fig) # Close the canvas
# print("Transposed Clustered Heatmap for " + "Input v sIP Complete..." + "\n")

################
# matrix = ALL_samples_log2.values

# Z = linkage(ALL_samples_log2, "ward")
# ZT = linkage(ALL_samples_log2.T, "ward")

# fig, ax = plt.subplots()
# plt.title("Dendrogram of Input vs IP" )
# plt.xlabel("sample")
# plt.ylabel("distance")
# dendrogram(
#     ZT,
#     show_leaf_counts=False,
#     leaf_rotation=90.,
#     leaf_font_size=12.,
#     show_contracted=True)

# fig.savefig("log2ClusteredDendrogram_Transposed_InputvsIP.png")
# plt.close(fig)
# print("Log2 Dendrogram for Input vs IP Complete...")

# idx_rows = leaves_list(Z)
# data = matrix[idx_rows, :]
# idx_columns = leaves_list(ZT)
# print("Clustering Complete...")
# print(".......................................")
# # idx_columns = range(10)
# data = data[:, idx_columns]

# # X = (data-np.average(data, axis=0))/np.std(data,axis=0)
# X = (data-np.average(data, axis=0))/np.std(data,axis=0)
# m = np.max(np.abs(X))
# print("Max of normalized data: " + str(m))
# Xmin = np.min(np.abs(X))
# print("Min of normalized data: " + str(Xmin))

# fig, ax = plt.subplots(figsize=(8, 4))
# plt.plot(X)
# ax.set_title("Log2 Plot Clustered X and Y for Input vs IP")
# ax.set_xlabel("Genes")
# ax.set_ylabel( "Expression" )
# plt.ylim(Xmin, 6)
# # ax.set_xticks(np.arange(0.5, X.shape[1]+0.5),)
# ax.set_xticklabels(idx_rows, rotation=50)
# ax.set_yticklabels(idx_columns, rotation=50)
# plt.tight_layout()
# fig.savefig("Log2_BasicPlot_InputvsIP.png" )
# plt.close(fig)
# print("Log2 Basic Plot for Input vs IP Complete...")
# print(".......................................")


# # fig, ax = plt.subplots(figsize=(20,3))
# # plt.plot(X.T)
# # ax.set_title("Log 2 Plot Transposed Clustered X and Y for " + "InputvsIP")
# # ax.set_xlabel( "?", rotation=45)
# # ax.set_ylabel( "?" )
# # plt.ylim(Xmin, m)
# # plt.tight_layout()
# # fig.savefig("Log2_BasicPlot_TransposedInputvsIP.png" )
# # plt.close(fig)
# # print("Log2 Transposed Basic Plot for " + "InputvsIPComplete...")

# fig, ax = plt.subplots(figsize=(6, 8))
# ax.set_title("Log 2 Heatmap of Clustered X and Y for " + "InputvsIP")
# im = ax.pcolor(X, cmap="viridis", vmin = -m, vmax = m- 0.5)
# ax.grid(False)
# ax.set_xticks(np.arange(0.5, X.shape[1]+0.5),)
# ax.set_xticklabels(idx_rows, rotation=50)
# ax.set_yticklabels(idx_columns, rotation=50)
# cbar = fig.colorbar(im, ax=ax)
# fig.subplots_adjust(left = 0.05, bottom = 0.15, right = 1.0, top = 0.95)
# fig.savefig("Log2_Clustered_Heatmap_"+ "InputvsIP.png" ) # Save the image
# plt.close(fig) # Close the canvas
# print("Creating HEATMAP NOW!!\nLog2 Clustered Heatmap for " + "Input vs IP Complete...")
# print(".......................................")

# df = pd.read_table(sys.argv[1], index_col=0)

# col_names = df.columns.values.tolist()
# row_names = df.index.values.tolist()

# pca = PCA( n_components =2)
# fit = pca.fit( df )
# x = fit.transform( df )

# fig, ax = plt.subplots(figsize=(6,6))

# rng = np.random.RandomState(0)
# colors = rng.rand(len(row_names))
# ax.scatter( x[:, 0], x[:, 1], c=colors , cmap="RdYlGn" )
# ax.set_title( "IP vs. Input PCA" )
# ax.set_xlabel( "PCA1")
# # plt.xlim(-10000,200000)
# # # plt.ylim(-1000, 1000)
# ax.set_ylabel( "PCA2" )
# ax.axhline(y=0, color='r', dashes=(1.0,1.0))
# plt.tight_layout()
# # plt.show()
# fig.savefig( "pca_Ahrens.png" )
# plt.close()
# print("\nPCA for Regular data for " + "Input vs IP Complete...")
# print("\nExplained Variance Ratio: ")
# print( fit.explained_variance_ratio_ )
# print("\nComponents: ")
# print( fit.components_ )
# print("\nShape of Components: ")
# print( fit.components_.shape )
# print("\n PCA data fit for all samples: ")
# print( x )
# print("\nShape PCA data fit: ")
# print( x.shape )
# print(".......................................")

# pca = PCA( n_components =2)
# fit = pca.fit( df.T )
# x = fit.transform( df.T )

# fig, ax = plt.subplots(figsize=(6,6))

# rng = np.random.RandomState(0)
# colors = rng.rand(len(col_names))
# ax.scatter( x[:, 0], x[:, 1], c=colors , cmap="RdYlGn" )
# ax.set_title( "Transposed IP vs. Input PCA" )
# ax.set_xlabel( "PCA1")
# # plt.xlim(-100,110)
# ax.set_ylabel( "PCA2" )
# ax.axhline(y=0, color='r', dashes=(1.0,1.0))
# plt.tight_layout()
# # plt.show()
# fig.savefig( "pca_Ahrens_transposed.png" )
# plt.close()
# print("\nPCA for Transposed data for " + "Input vs IP Complete...\n")
# print("\nExplained Variance Ratio")
# print( fit.explained_variance_ratio_ )
# print("\nComponents")
# print( fit.components_ )
# print("\nShape of Components")
# print( fit.components_.shape )
# print("\n PCA data fit for all samples ")
# print( x )
# print("\nShape PCA data fit")
# print( x.shape )
# print(".......................................")


# ALL_samples_log2
# col_names = ALL_samples_log2.columns.values.tolist()
# row_names = ALL_samples_log2.index.values.tolist()

# pca = PCA( n_components =2)
# fit = pca.fit( ALL_samples_log2 )
# x = fit.transform( ALL_samples_log2 )
# # print( "Log2 Transformed data \n" + fit.explained_variance_ratio_ )

# fig, ax = plt.subplots(figsize=(6,6))

# rng = np.random.RandomState(0)
# colors = rng.rand(len(row_names))
# ax.scatter( x[:, 0], x[:, 1], c=colors , cmap="RdYlGn" )
# ax.set_title( "IP vs. InputLog2 Transformed data PCA" )
# ax.set_xlabel( "PCA1")
# # plt.xlim(-10000,200000)
# # # plt.ylim(-1000, 1000)
# ax.set_ylabel( "PCA2" )
# ax.axhline(y=0, color='r', dashes=(1.0,1.0))
# plt.tight_layout()
# # plt.show()
# fig.savefig( "pca_Ahrens_log2.png" )
# plt.close()
# print("\nPCA for Log2 Transformed data " + "Input vs IP Complete...\n")
# print("\nExplained Variance Ratio")
# print( fit.explained_variance_ratio_ )
# print("\nComponents")
# print( fit.components_ )
# print("\nShape of Components")
# print( fit.components_.shape )
# print("\n PCA data fit for all samples ")
# print( x )
# print("\nShape PCA data fit")
# print( x.shape )
# print(".......................................")



# pca = PCA( n_components =2)
# fit = pca.fit( ALL_samples_log2.T )
# x = fit.transform( ALL_samples_log2.T )
# # print( "Log2 Transformed and Transposed data \n " + fit.explained_variance_ratio_ )

# fig, ax = plt.subplots(figsize=(6,6))

# rng = np.random.RandomState(0)
# colors = rng.rand(len(col_names))
# ax.scatter( x[:, 0], x[:, 1], c=colors , cmap="RdYlGn" )
# ax.set_title( "Log2 Transformed and Transposed PCA" )
# ax.set_xlabel( "PCA1")
# # plt.xlim(-100,110)
# ax.set_ylabel( "PCA2" )
# ax.axhline(y=0, color='r', dashes=(1.0,1.0))
# plt.tight_layout()
# # plt.show()
# fig.savefig( "pca_Ahrens_transposed_log2.png" )
# plt.close()
# print("\nPCA for Log2 Transformed and Transposed data " + "Input vs IP Complete...\n")
# print("\nExplained Variance Ratio")
# print( fit.explained_variance_ratio_ )
# print("\nComponents")
# print( fit.components_ )
# print("\nShape of Components")
# print( fit.components_.shape )
# print("\n PCA data fit for all samples ")
# print( x )
# print("\nShape PCA data fit")
# print( x.shape )
# print(".......................................")

d = "GeneOntology/"
if not os.path.exists(d):
    os.makedirs(d)

os.chdir(d)     

with open("increasedglutGABA.out", 'w') as file:
    for item in increasedNP:
        file.write(str(item).upper() + "\n")

with open("decreasedglutGABA.out", 'w') as file:
    for item in decreasedNP:
        file.write(str(item).upper() + "\n")

with open("Meanslog2_inc_glutGABA.out", 'w') as file:
    for item in Meanslog2_inc:
        file.write(str(item).upper() + "\n")

with open("Meanslog2_dec_glutGABA.out", 'w') as file:
    for item in Meanslog2_dec:
        file.write(str(item).upper() + "\n")

with open("Means_inc_glutGABA.out", 'w') as file:
    for item in Means_inc:
        file.write(str(item).upper() + "\n")

with open("Means_dec_glutGABA.out", 'w') as file:
    for item in Means_dec:
        file.write(str(item).upper() + "\n")


#######################
# some = int(len( row_names ))


# for i in range( some ):
#
#      ax.annotate( row_names[i], ( x[i,0], x[i,1] ))


# for i in range(len(row_names)):
#     xaxis = col_names
#     data = df.T.iloc[:,i]
#     plt.plot(data)
#     plt.ylim(0, 150)
#     plt.xlabel('InputvsIP')
#     plt.ylabel('FPKM?')
#     plt.title("InputvsIP vs mRNA expression")
#     plt.plot(color=colors)
#     plt.savefig("InputvsIP" + str(i) +".png")
#     plt.close()
#     continue

# fig, ax = plt.subplots(figsize=(6, 8))
# ax.set_title("Log 2 Transposed Heatmap of Clustered X and Y for " + "InputvsIP")
# im = ax.pcolor(X.T, cmap="cool", vmin = -m, vmax = m- 0.5)
# ax.grid(False)
# ax.set_xticks(np.arange(0.5, X.shape[1]+0.5),)
# ax.set_xticklabels(idx_rows, rotation=50)
# ax.set_yticklabels(idx_columns, rotation=50)
# cbar = fig.colorbar(im, ax=ax)
# fig.subplots_adjust(left = 0.05, bottom = 0.15, right = 1.0, top = 0.95)
# fig.savefig("Log2_Clustered_Heatmap_Transposed_"+ "InputvsIP.png" ) # Save the image
# plt.close(fig) # Close the canvas
# print("Log 2 Transposed Clustered Heatmap for " + "InputvsIP Complete..." + "\n\n")

#
# # for all Input
# for i in range(len(sample_types)):
#
#     print(sample_types[10:])
#     matrix = df.iloc[:,10:]
#     matrix = matrix.values
#     i = " input "
#
#     Z = linkage(matrix, "ward")
#     ZT = linkage(matrix.T, "ward")
#
#     fig, ax = plt.subplots()
#     plt.title("Dendrogram of " + str(i) )
#     plt.xlabel("sample")
#     plt.ylabel("distance")
#     dendrogram(
#         ZT,
#         show_leaf_counts=False,
#         leaf_rotation=90.,
#         leaf_font_size=12.,
#         show_contracted=True)
#
#     fig.savefig("ClusteredDendrogram_Transposed" + str(i) + ".png")
#     plt.close(fig)
#     print("Dendrogram for " + str(i) + " Complete...")
#
#
#     idx_rows = leaves_list(Z)
#     data = matrix[idx_rows, :]
#     idx_columns = leaves_list(ZT)
#     print("Clustering Complete...")
#     # idx_columns = range(10)
#     data = data[:, idx_columns]
#
#     X = (data-np.average(data, axis=0))/np.std(data,axis=0)
#
#     m = np.max(np.abs(X))
#     print("Max of normalized data: " + str(m))
#     Xmin = np.min(np.abs(X))
#     print("Min of normalized data: " + str(Xmin) +"\n")
#
#     idx_rows = leaves_list(Z)
#     data = matrix[idx_rows, :]
#     idx_columns = leaves_list(ZT)
#     # idx_columns = range(10)
#     data = data[:, idx_columns]
#
#     fig, ax = plt.subplots()
#     plt.plot(X)
#     ax.set_title("Plot Clustered X and Y for " + str(i))
#     ax.set_xlabel( "?", rotation=45)
#     ax.set_ylabel( "Expression" )
#     plt.ylim(0, 500)
#     plt.tight_layout()
#     fig.savefig("BasicPlot" + str(i) + ".png" )
#     plt.close(fig)
#     print("Basic Plot for " + str(i) + " Complete...")
#
#
#     fig, ax = plt.subplots()
#     plt.plot(X.T)
#     ax.set_title("Plot Transposed Clustered X and Y for " + str(i))
#     ax.set_xlabel( "?", rotation=45)
#     ax.set_ylabel( "Expression" )
#     plt.ylim(0, 500)
#     plt.tight_layout()
#     fig.savefig("BasicPlot_Transposed"+ str(i) + ".png" )
#     plt.close(fig)
#     print("Transposed Basic Plot for " + str(i) + " Complete...")
#
#
#     fig, ax = plt.subplots()
#     ax.set_title("Heatmap of Clustered X and Y for " + str(i))
#     im = ax.pcolor(X, cmap="viridis", vmin = -m, vmax = m- 0.5)
#     ax.grid(False)
#     ax.set_xticks(np.arange(0.5, X.shape[1]+0.5),)
#     ax.set_xticklabels(idx_columns, rotation=50)
#     ax.set_yticklabels(idx_rows, rotation=50)
#     cbar = fig.colorbar(im, ax=ax)
#     fig.subplots_adjust(left = 0.05, bottom = 0.15, right = 1.0, top = 0.95)
#     fig.savefig("Clustered_Heatmap_" + str(i) + ".png" ) # Save the image
#     plt.close(fig) # Close the canvas
#     print("Clustered Heatmap for " + str(i) + " Complete...")
#
#
#     fig, ax = plt.subplots()
#     ax.set_title("Transposed Heatmap of Clustered X and Y for " + str(i))
#     im = ax.pcolor(X.T, cmap="viridis", vmin = -m, vmax = m- 0.5)
#     ax.grid(False)
#     ax.set_xticks(np.arange(0.5, X.shape[1]+0.5),)
#     ax.set_xticklabels(idx_columns, rotation=50)
#     ax.set_yticklabels(idx_rows, rotation=50)
#     cbar = fig.colorbar(im, ax=ax)
#     fig.subplots_adjust(left = 0.05, bottom = 0.15, right = 1.0, top = 0.95)
#     fig.savefig("Clustered_Heatmap_Transposed_"+ str(i) + ".png" ) # Save the image
#     plt.close(fig) # Close the canvas
#     print("Transposed Clustered Heatmap for " + str(i) + " Complete..." + "\n\n")
#     break
#
#
#
# for i in range(len(sample_types)):
#
#     print(sample_types[i:10])
#     matrix = df.iloc[:,i:10]
#     matrix["All IP mean"] = matrix.mean(axis=1)
#     print(matrix["All IP mean"])
#     matrix = matrix.values
#     i = " ALL IP "
#
#     Z = linkage(matrix, "ward")
#     ZT = linkage(matrix.T, "ward")
#
#     idx_rows = leaves_list(Z)
#     data = matrix[idx_rows, :]
#     idx_columns = leaves_list(ZT)
#     print("Clustering Complete...")
#     # idx_columns = range(10)
#     data = data[:, idx_columns]
#
#     X = np.log2((data-np.average(data, axis=0))/np.std(data,axis=0))
#     # X = (data-np.average(data, axis=0))/np.std(data,axis=0)
#     print(X)
#     m = np.max(np.abs(X))
#     print("Max of normalized data: " + str(m))
#     Xmin = np.min(np.abs(X))
#     print("Min of normalized data: " + str(Xmin) +"\n")
#     break
#
# All_IP = X
#
# for i in range(len(sample_types)):
#
#     print(sample_types[10:])
#     matrix = df.iloc[:,10:]
#     matrix["All INPUT mean"] = matrix.mean(axis=1)
#     print(matrix["All INPUT mean"])
#     matrix = matrix.values
#     i = " input "
#     break
#
#     Z = linkage(matrix, "ward")
#     ZT = linkage(matrix.T, "ward")
#
#     idx_rows = leaves_list(Z)
#     data = matrix[idx_rows, :]
#     idx_columns = leaves_list(ZT)
#     print("Clustering Complete...")
#     # idx_columns = range(10)
#     data = data[:, idx_columns]
#
#     X = np.log2((data-np.average(data, axis=0))/np.std(data,axis=0))
#     # X = (data-np.average(data, axis=0))/np.std(data,axis=0)
#     print(X)
#     m = np.max(np.abs(X))
#     print("Max of normalized data: " + str(m))
#     Xmin = np.min(np.abs(X))
#     print("Min of normalized data: " + str(Xmin) +"\n")
#     break
#
# All_Input = X
#
# print(All_IP)
# print(All_Input)









# # For input_male
# for i in range(len(sample_types)):
#
#     print(sample_types[i])
#     matrix = df.iloc[:,11,12]
#     print(matrix)
#     matrix = matrix.values
#     print(matrix)
#
#     Z = linkage(matrix, "ward")
#     ZT = linkage(matrix.T, "ward")
#
#     fig, ax = plt.subplots()
#     plt.title("Dendrogram of " + str(i) )
#     plt.xlabel("sample")
#     plt.ylabel("distance")
#     dendrogram(
#         ZT,
#         show_leaf_counts=False,
#         leaf_rotation=90.,
#         leaf_font_size=12.,
#         show_contracted=True)
#
#     fig.savefig("ClusteredDendrogram_Transposed" + str(i) + ".png")
#     plt.close(fig)
#     print("Dendrogram for " + str(i) + "Complete...")
#
#
#     idx_rows = leaves_list(Z)
#     data = matrix[idx_rows, :]
#     idx_columns = leaves_list(ZT)
#     print("Clustering Complete...")
#     # idx_columns = range(10)
#     data = data[:, idx_columns]
#
#     X = (data-np.average(data, axis=0))/np.std(data,axis=0)
#
#     m = np.max(np.abs(X))
#     print("Max of normalized data: " + str(m))
#     Xmin = np.min(np.abs(X))
#     print("Min of normalized data: " + str(Xmin) +"\n")
#
#     idx_rows = leaves_list(Z)
#     data = matrix[idx_rows, :]
#     idx_columns = leaves_list(ZT)
#     # idx_columns = range(10)
#     data = data[:, idx_columns]
#
#     fig, ax = plt.subplots()
#     plt.plot(X)
#     ax.set_title("Plot Clustered X and Y for " + str(i))
#     ax.set_xlabel( "?", rotation=45)
#     ax.set_ylabel( "Expression" )
#     # plt.ylim(0, 500)
#     plt.tight_layout()
#     fig.savefig("BasicPlot" + str(i) + ".png" )
#     plt.close(fig)
#     print("Basic Plot for " + str(i) + "Complete...")
#
#
#     fig, ax = plt.subplots()
#     plt.plot(X.T)
#     ax.set_title("Plot Transposed Clustered X and Y for " + str(i))
#     ax.set_xlabel( "?", rotation=45)
#     ax.set_ylabel( "Expression" )
#     # plt.ylim(0, 500)
#     plt.tight_layout()
#     fig.savefig("BasicPlot_Transposed"+ str(i) + ".png" )
#     plt.close(fig)
#     print("Transposed Basic Plot for " + str(i) + "Complete...")
#
#
#     fig, ax = plt.subplots()
#     ax.set_title("Heatmap of Clustered X and Y for " + str(i))
#     im = ax.pcolor(X, cmap="viridis", vmin = -m, vmax = m- 0.5)
#     ax.grid(False)
#     ax.set_xticks(np.arange(0.5, X.shape[1]+0.5),)
#     ax.set_xticklabels(idx_columns, rotation=50)
#     ax.set_yticklabels(idx_rows, rotation=50)
#     cbar = fig.colorbar(im, ax=ax)
#     fig.subplots_adjust(left = 0.05, bottom = 0.15, right = 1.0, top = 0.95)
#     fig.savefig("Clustered_Heatmap_" + str(i) ".png" ) # Save the image
#     plt.close(fig) # Close the canvas
#     print("Clustered Heatmap for " + str(i) + "Complete...")
#
#
#     fig, ax = plt.subplots()
#     ax.set_title("Transposed Heatmap of Clustered X and Y for " + str(i))
#     im = ax.pcolor(X.T, cmap="viridis", vmin = -m, vmax = m- 0.5)
#     ax.grid(False)
#     ax.set_xticks(np.arange(0.5, X.shape[1]+0.5),)
#     ax.set_xticklabels(idx_columns, rotation=50)
#     ax.set_yticklabels(idx_rows, rotation=50)
#     cbar = fig.colorbar(im, ax=ax)
#     fig.subplots_adjust(left = 0.05, bottom = 0.15, right = 1.0, top = 0.95)
#     fig.savefig("Clustered_Heatmap_Transposed_"+ str(i) + ".png" ) # Save the image
#     plt.close(fig) # Close the canvas
#     print("Transposed Clustered Heatmap for " + str(i) + "Complete..." + "\n\n")
#     break
    
#
#
# matrix = df.values
#
# Z = linkage(matrix, "ward")
# ZT = linkage(matrix.T, "ward")
#
# fig, ax = plt.subplots()
# plt.title("Dendrogram of " + sys.argv[1] )
# plt.xlabel("sample")
# plt.ylabel("distance")
# dendrogram(
#     ZT,
#     show_leaf_counts=False,
#     leaf_rotation=90.,
#     leaf_font_size=12.,
#     show_contracted=True)
#
# fig.savefig("ClusteredDendrogram_Transposed" + sys.argv[1] + ".png")
# plt.close(fig)
#
#
# idx_rows = leaves_list(Z)
# data = matrix[idx_rows, :]
# idx_columns = leaves_list(ZT)
# print("Clustering Complete...")
# # idx_columns = range(10)
# data = data[:, idx_columns]
#
# X = (data-np.average(data, axis=0))/np.std(data,axis=0)
#
# m = np.max(np.abs(X))
# print("Max of normalized data: " + str(m))
# Xmin = np.min(np.abs(X))
# print("Min of normalized data: " + str(Xmin) +"\n")
#
# idx_rows = leaves_list(Z)
# data = matrix[idx_rows, :]
# idx_columns = leaves_list(ZT)
# # idx_columns = range(10)
# data = data[:, idx_columns]
#
# fig, ax = plt.subplots()
# plt.plot(X)
# ax.set_title("Plot Clustered X and Y")
# ax.set_xlabel( "?", rotation=45)
# ax.set_ylabel( "Expression" )
# # plt.ylim(0, 500)
# plt.tight_layout()
# fig.savefig("BasicPlot" + sys.argv[1] + ".png" )
# plt.close(fig)
# print("Basic Plot Complete...")
#
#
# fig, ax = plt.subplots()
# plt.plot(X.T)
# ax.set_title("Plot Transposed Clustered X and Y")
# ax.set_xlabel( "?", rotation=45)
# ax.set_ylabel( "Expression" )
# # plt.ylim(0, 500)
# plt.tight_layout()
# fig.savefig("BasicPlot_Transposed"+ sys.argv[1] + ".png" )
# plt.close(fig)
# print("Transposed Basic Plot Complete...")
#
#
# fig, ax = plt.subplots()
# ax.set_title("Heatmap of Clustered X and Y")
# im = ax.pcolor(X, cmap="viridis", vmin = -m, vmax = m- 0.5)
# ax.grid(False)
# ax.set_xticks(np.arange(0.5, X.shape[1]+0.5),)
# ax.set_xticklabels(idx_columns, rotation=50)
# ax.set_yticklabels(idx_rows, rotation=50)
# cbar = fig.colorbar(im, ax=ax)
# fig.subplots_adjust(left = 0.05, bottom = 0.15, right = 1.0, top = 0.95)
# fig.savefig("Clustered_Heatmap_"+ sys.argv[1] + ".png" ) # Save the image
# plt.close(fig) # Close the canvas
# print("Clustered Heatmap Complete...")
#
#
# fig, ax = plt.subplots()
# ax.set_title("Transposed Heatmap of Clustered X and Y")
# im = ax.pcolor(X.T, cmap="viridis", vmin = -m, vmax = m- 0.5)
# ax.grid(False)
# ax.set_xticks(np.arange(0.5, X.shape[1]+0.5),)
# ax.set_xticklabels(idx_columns, rotation=50)
# ax.set_yticklabels(idx_rows, rotation=50)
# cbar = fig.colorbar(im, ax=ax)
# fig.subplots_adjust(left = 0.05, bottom = 0.15, right = 1.0, top = 0.95)
# fig.savefig("Clustered_Heatmap_Transposed_"+ sys.argv[1] + ".png" ) # Save the image
# plt.close(fig) # Close the canvas
# print("Transposed Clustered Heatmap Complete..." + "\n\n")

    
    












# timepoints = set(df1[0].values.tolist())
# timepoints =list(timepoints)
#
#
#
#
# # for i in range(len(timepoints)):
# #     df_t1 = df1[(df1[0] == timepoints[i] )]
# #     df_i = df_t1[[1, 3]]
# #
# #     df_i.columns = ["Gene", timepoints[i] + " Depth"]
# #
# #     continue
# #
#
#
# df_t1 = df1[(df1[0] == timepoints[0] )]
# df_t1 = df_t1[[1, 3]]
# df_t1.columns = ["Gene", timepoints[0]]
# df_t1.index = df_t1["Gene"]
# df_t1 = df_t1[timepoints[0]]
# df_t1.columns = [timepoints[0]]
#
# df_t2 = df1[(df1[0] == timepoints[1] )]
# df_t2 = df_t2[[1, 3]]
# df_t2.columns = ["Gene", timepoints[1]]
# df_t2.index = df_t2["Gene"]
# df_t2 = df_t2[timepoints[1]]
# df_t2.columns = [timepoints[1]]
#
# df_t3 = df1[(df1[0] == timepoints[2] )]
# df_t3 = df_t3[[1, 3]]
# df_t3.columns = ["Gene", timepoints[2]]
# df_t3.index = df_t3["Gene"]
# df_t3 = df_t3[timepoints[2]]
# df_t3.columns = [timepoints[2]]
#
# df_t4 = df1[(df1[0] == timepoints[3] )]
# df_t4 = df_t4[[1, 3]]
# df_t4.columns = ["Gene", timepoints[3]]
# df_t4.index = df_t4["Gene"]
# df_t4 = df_t4[timepoints[3]]
# df_t4.columns = [timepoints[3]]
#
# df_t5 = df1[(df1[0] == timepoints[4] )]
# df_t5 = df_t5[[1, 3]]
# df_t5.columns = ["Gene", timepoints[4]]
# df_t5.index = df_t5["Gene"]
# df_t5 = df_t5[timepoints[4]]
# df_t5.columns = [timepoints[4]]
#
# df_t6 = df1[(df1[0] == timepoints[5] )]
# df_t6 = df_t6[[1, 3]]
# df_t6.columns = ["Gene", timepoints[5]]
# df_t6.index = df_t6["Gene"]
# df_t6 = df_t6[timepoints[5]]
# df_t6.columns = [timepoints[5]]
#
# df_t7 = df1[(df1[0] == timepoints[6] )]
# df_t7 = df_t7[[1, 3]]
# df_t7.columns = ["Gene", timepoints[6]]
# df_t7.index = df_t7["Gene"]
# df_t7 = df_t7[timepoints[6]]
# df_t7.columns = [timepoints[6]]
#
# df_t8 = df1[(df1[0] == timepoints[7] )]
# df_t8 = df_t8[[1, 3]]
# df_t8.columns = ["Gene", timepoints[7]]
# df_t8.index = df_t8["Gene"]
# df_t8 = df_t8[timepoints[7]]
# df_t8.columns = [timepoints[7]]
#
# df_t9 = df1[(df1[0] == timepoints[8] )]
# df_t9 = df_t9[[1,3]]
# df_t9.columns = ["Gene", timepoints[8]]
# df_t9.index = df_t9["Gene"]
# df_t9 = df_t9[timepoints[8]]
# df_t9.columns = [timepoints[8]]
#
#
# df_t10 = df1[(df1[0] == timepoints[9] )]
# df_t10 = df_t10[[1, 3]]
# df_t10.columns = ["Gene", timepoints[9]]
# df_t10.index = df_t10["Gene"]
# df_t10 = df_t10[timepoints[9]]
# df_t10.columns = [timepoints[9]]
#
# # df_t1 = df_t6.loc["t6 Depth"]
#
#
#
# frames = [df_t1, df_t2, df_t3, df_t4, df_t5, df_t6, df_t7, df_t8, df_t9, df_t10]

#
# data = dfRegrex.values
# matrix = result.values
#
# Z = linkage(data, "ward")
# ZT = linkage(data.T, "ward")
#
#
# idx_rows = leaves_list(Z)
# data = matrix[idx_rows, :]
# # idx_columns = leaves_list(ZT)
# # idx_columns = range(10)
# # data = data[:, idx_columns]
#
# # fig, ax = plt.subplots(figsize=(15, 3))
# # plt.plot(dfRegrex.T)
# # # ax.set_xlabel( "Gene", rotation=45)
# # ax.set_ylabel( "Depth" )
# # plt.ylim(0, 500)
# # plt.tight_layout()
# # fig.savefig( "DepthvsPositionEcoliTimeSeries_verts" + arg + ".png" )
# # plt.close()
#
#
# X = (data-np.average(data, axis=0))/np.std(data,axis=0)
#
# m = np.max(np.abs(X))
#
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
#
#
#
# fig, ax = plt.subplots(figsize=(15, 3))
# ax.set_title("Depth of All genes over time")
# plt.plot(dfRegrex)
# ax.set_xlabel( "Gene" )
# ax.set_ylabel( "Depth" )
# plt.ylim(0, 500)
# plt.tight_layout()
# fig.savefig( "Depth_ParallelCoordinates_All.png" )
# plt.close()
#
# increased = dfALL_samples[(df_ALL_["t10"] > dfRegrex["t1"])]
# list_increased = increased.index.tolist()
# print(list_increased)
# fig, ax = plt.subplots(figsize=(15, 3))
# ax.set_title("Depth of Increased genes over time")
# plt.plot(increased)
# ax.set_xlabel( "Gene", rotation=45)
# ax.set_ylabel( "Depth" )
# plt.ylim(0, 500)
# plt.tight_layout()
# fig.savefig( "Depth_ParallelCoodinates_increased.png" )
# plt.close()
#
#
# decreased = dfRegrex[(dfRegrex["t1"] > dfRegrex["t10"])]
# list_decreased = decreased.index.tolist()
# print(list_decreased)
# fig, ax = plt.subplots(figsize=(15, 3))
# ax.set_title("Depth of Decreased genes over time")
# plt.plot(decreased)
# ax.set_xlabel( "Gene" )
# ax.set_ylabel( "Depth" )
# plt.ylim(0, 500)
# plt.tight_layout()
# fig.savefig( "Depth_ParallelCoordinates_decreased.png" )
# plt.close()
#
# neutral =  dfRegrex[(dfRegrex["t1"] < dfRegrex["t10"]) & (dfRegrex["t10"] > dfRegrex["t1"]) ]
# list_neutral = neutral.index.tolist()
# print(list_neutral)
#
#
#
#
# # filtered = result[(result < 400)]
# # sliced = data[(data.Position > start) & (data.Position < end)]
# # print(filtered)
#
#
# # plt.figure()
# # plt.plot(filtered)
# # plt.show()
#
# # pd.plotting.parallel_coordinates(result, result.index)
#
# # plt.show()
#
# # print(filtered)
# # data = filtered.values
# # filtered.dropna()
#
#
# data = increased.values
# Z = linkage(data, "ward")
# idx_rows = leaves_list(Z)
# data = data[idx_rows, :]
# X = (data-np.average(data, axis=0))/np.std(data,axis=0)
#
# m = np.max(np.abs(X))
#
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
#
# data = decreased.values
# Z = linkage(data, "ward")
# idx_rows = leaves_list(Z)
# data = data[idx_rows, :]
# X = (data-np.average(data, axis=0))/np.std(data,axis=0)
#
# m = np.max(np.abs(X))
#
# fig, ax = plt.subplots(figsize=(15, 6))
#
#
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
#
#
#
