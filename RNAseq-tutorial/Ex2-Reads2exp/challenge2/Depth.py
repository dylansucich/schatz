#!/usr/bin/env python3

#    compute average depth per exon from the reference gene list
#      Depth.py refgenes.ptt t* > depth.out

import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

geneExon ={}
start_end = {}
gene_locations = {}
gene_list = []
start_list = []
end_list = []
exon_length = []

reference_genes = pd.read_table(sys.argv[1], header=None)

# print(reference_genes)

for line in open( sys.argv[1] ):
    info = line.strip("\r\n").split()
    start1 = int(info[0].split("..")[0])
    end1 = int(info[0].split("..")[1])
    exon_length1 = int(end1 - start1)
    start = start_list.append(int(info[0].split("..")[0]))
    end = end_list.append(int(info[0].split("..")[1]))
    gene = info[4]
    # gene = gene_list.append(info[4])

    start_end[gene] = start1, end1
    geneExon[gene] = exon_length1

    continue
    # for key, value in start_end.items():
print(geneExon)   
        
        
        
    
    # gene dict(zip(list_with_keys, list_with_values))
    
    
    #
    # exon_length.append(exon_length1)
    # gene_locations_startend[start] = end
    # gene_locations[gene_list] = gene_locations_startend



# people = {1: {"gene" : gene_list, "start" : start_list,"end": end_list,
#     "exon length": exon_length}}

# people[3] = {}
#
# people[3]['name'] = 'Luna'
# people[3]['age'] = '24'
# people[3]['sex'] = 'Female'
# people[3]['married'] = 'No'



# startdf = pd.DataFrame(people[1]["start"])
#
# fig, ax = plt.subplots(figsize=(15,3))
# ax.plot(data1)
# x1position = people[1]["start"]
# for xc in x1position:
#     plt.axvline(x=xc, color='b', linestyle='-')
#
# x1position = people[1]["end"]
# for xc in x1position:
#     plt.axvline(x=xc, color='r', linestyle='--')
#
# plt.show()
#
#
#
# exonsdf = pd.DataFrame({"gene" : gene_list, "start" : start_list,"end": end_list,
#     "exon length": exon_length})



time_point_depth = {}
gene_depth = []
avg = []
genes = []
exp = []
for arg in sys.argv[2:]:
    df = pd.read_table(arg, header= None)
    df.columns = ["Species", "Position", "Depth"]
    df.set_index("Position")
    data = df[["Position", "Depth"]]
    data.set_index("Position")
    
    for gene, (start, end) in start_end.items():
        data[(data.Position > start) & (data.Position < end)]
        sliced = data[(data.Position > start) & (data.Position < end)]
        val = arg.split("/")[-1]
        print(val.split("_")[0]+ "\t" + gene + "\t" + "\t" + str((np.sum(sliced["Depth"],axis=0)/float(geneExon[gene]))))
        # avg = avg.append(sliced["Depth"].mean())
#         genes = genes.append(gene)
#         exp = exp.append(arg.split("_")[0])
#
#         # gene_depth[int(sliced["Depth"].mean())] = gene, arg.split("_")[0]
#
#         filtered = pd.DataFrame({"exp" : exp, "gene" : genes, "Mean Depth": avg})
        
        continue
    

    
# #     time_points.append(arg.split("_")[0])
# #
# #     fields = line.rstrip("\r\n").split("\t")
# #     fields[2]
# #     gene_start = int(fields[3])
# #     gene_end = int(fields[4])
# #     gene_id = fields[8]
#
#
#
#
#
#
#

# #     # print(data.iloc[0:,0])
# #
# #
# #     if data[0] > 190:
# #         print(data[0])
# #
# #
# #
# #
# #
# #
# #     if exonsdf["start"].isin(data["Position"]) == True:
# #         filtered = data.iloc[(data["Position"] >= exonsdf["start"]) & (data["Position"] <= exonsdf["end"])]
# #         print(filtered)
# #     #
# #     # int(exonsdf["start"][gene]
# #         continue
#
#     # for gene in enumerate(exonsdf):
# #         print(gene)
# #         if Position in data <= exonsdf["start"]:
# #             print(data["Position"])
# #
#
#
#
#     fig, ax = plt.subplots(figsize=(12,6))
#     ax.plot(data, c="r")
# # #     ax.plot(], c="r")
# # #
# # #     # for i in range(0, ,exonsdf[["start"]]), 2):
# # # #         plt.plot(exonsdf[["start"]][i:i+2], exonsdf[["end"]][i:i+2], 'ro-')
# # #
#     ax.set_title( "Depth vs. Position for " + arg)
#     ax.set_xlabel( "Position")
#     plt.xlim(-1000,110000)
#     ax.set_ylabel( "Depth" )
#     plt.ylim(0,1500)
#     # ax.autoscale_view()
#     plt.tight_layout()
#     fig.savefig( "DepthvsPositionEcoliTimeSeries_verts" + arg + ".png" )
#     plt.close()
#     continue
#
# # starts = exonsdf[["gene", "start"]]
# # ax.plot(exonsdf[["start", "end"]], c="r")
#
#
# # rng = np.random.RandomState(0)
# # colors = rng.rand(101273)
# # plt.plot(color=colors, cmap='spring')
# # fig.savefig( "DepthvsPositionEcoliTimeSeries_loop.png" )
# # plt.close()
#
#     # for gene, start, end in exons:
#
#
#
#
#
#
#
#
#
# # for line in enumerate(gene_start_end):
# #     print(line)
#
#
#
#
#
#
# # for line in enumerate(reference_genes):
# #     fields = line.rstrip("\r\n")
#
#
#
#     #dictionary gene
#
# # depth_file = pd.read_table(sys.argv[2])
#
#     # depth and number of bases covered
#     # in any given tissue, only about half
#     # coverage sum of ind depths / length of the
#     # some genes no genes