#!/usr/bin/env python3

#    compute average depth per exon from the reference gene list
#      DepthperExon.py refgenes.ptt t#.

import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


gene_locations_startend = {}
gene_locations = {}
gene_list = []
start_list = []
end_list = []
exon_length = []


for line in open( sys.argv[1] ):
    info = line.strip("\r\n").split()
    start1 = int(info[0].split("..")[0])
    end1 = int(info[0].split("..")[1])
    exon_length1 = int(end1 - start1)
    start = start_list.append(int(info[0].split("..")[0]))
    end = end_list.append(int(info[0].split("..")[1]))
    exon_length.append(exon_length1)
    gene = gene_list.append(info[4])
    gene_locations_startend[start] = end

people = {1: {"gene" : gene_list, "start" : start_list,"end": end_list,
    "exon length": exon_length}}
    
startdf = pd.DataFrame(people[1]["start"])


time_points = list(sys.argv[2:])
print(time_points)

df1 = pd.read_table(sys.argv[2], header=None)
df2 = pd.read_table(sys.argv[3], header=None)
df3 = pd.read_table(sys.argv[4], header=None)
df4 = pd.read_table(sys.argv[5], header=None)
df5 = pd.read_table(sys.argv[6], header=None)
df6 = pd.read_table(sys.argv[7], header=None)
df7 = pd.read_table(sys.argv[8], header=None)
df8 = pd.read_table(sys.argv[9], header=None)
df9 = pd.read_table(sys.argv[10], header=None)
df10 = pd.read_table(sys.argv[11], header=None)

df1.columns = ["Species", "Position", "Depth"]
df2.columns = ["Species", "Position", "Depth"]
df3.columns = ["Species", "Position", "Depth"]
df4.columns = ["Species", "Position", "Depth"]
df5.columns = ["Species", "Position", "Depth"]
df6.columns = ["Species", "Position", "Depth"]
df7.columns = ["Species", "Position", "Depth"]
df8.columns = ["Species", "Position", "Depth"]
df9.columns = ["Species", "Position", "Depth"]
df10.columns = ["Species", "Position", "Depth"]

data1 = df1[["Position", "Depth"]]
data2 = df2[["Position", "Depth"]]
data3 = df3[["Position", "Depth"]]
data4 = df4[["Position", "Depth"]]
data5 = df5[["Position", "Depth"]]
data6 = df6[["Position", "Depth"]]
data7 = df7[["Position", "Depth"]]
data8 = df8[["Position", "Depth"]]
data9 = df9[["Position", "Depth"]]
data10 = df10[["Position", "Depth"]]

fig, ax = plt.subplots(figsize=(15,3))

ax.plot(data1)
ax.plot(data2)
ax.plot(data3)
ax.plot(data4)
ax.plot(data5)
ax.plot(data6)
ax.plot(data7)
ax.plot(data8)
ax.plot(data9)
ax.plot(data10)

x1position = people[1]["start"]
for xc in x1position:
    plt.axvline(x=xc, color='b', linestyle='-')
    
x2position = people[1]["end"]
for xc in x2position:
    plt.axvline(x=xc, color='r', linestyle='--')

ax.set_title( "Depth vs. Position for E. coli RNAseq Time-series")
ax.set_xlabel( "Position")
plt.xlim(0,300)
ax.set_ylabel( "Depth" )
plt.ylim(0,180)
# ax.autoscale_view()
plt.tight_layout()

rng = np.random.RandomState(0)
colors = rng.rand(101273)
plt.plot(color=colors, cmap='spring')
fig.savefig( "Zoomed_DepthvsPositionEcoliTimeSeries_Vert.png" )
plt.show()
plt.close()

# plt.subplot(10, 1, 1)
# plt.plot(data1)
# plt.xlim(13200,17000)
# plt.ylim(0,180)
# plt.plot(color=colors, cmap='spring')
# plt.title('Timepoint 1')

# plt.subplot(10, 1, 2)
# plt.plot(data2)
# plt.xlim(13200,17000)
# plt.ylim(0,180)
# plt.plot(color=colors, cmap='spring')
# plt.title('Timepoint 2')

# plt.subplot(10, 1, 3)
# plt.plot(data3)
# plt.xlim(13200,17000)
# plt.ylim(0,180)
# plt.plot(color=colors, cmap='spring')
# plt.title('Timepoint 3')

# plt.subplot(10, 1, 4)
# plt.plot(data4)
# plt.xlim(13200,17000)
# plt.ylim(0,180)
# plt.plot(color=colors, cmap='spring')
# plt.title('Timepoint 4')

# plt.subplot(10, 1, 5)
# plt.plot(data5)
# plt.xlim(13200,17000)
# plt.ylim(0,180)
# plt.plot(color=colors, cmap='spring')
# plt.title('Timepoint 5')

# plt.subplot(10, 2, 1)
# plt.plot(data6)
# plt.xlim(13200,17000)
# plt.ylim(0,180)
# plt.plot(color=colors, cmap='spring')
# plt.title('Timepoint 6')

# plt.subplot(10, 2, 2)
# plt.plot(data7)
# plt.xlim(13200,17000)
# plt.ylim(0,180)
# plt.plot(color=colors, cmap='spring')
# plt.title('Timepoint 7')

# plt.subplot(10, 2, 3)
# plt.plot(data8)
# plt.xlim(13200,17000)
# plt.ylim(0,180)
# plt.plot(color=colors, cmap='spring')
# plt.title('Timepoint 8')

# plt.subplot(10, 2, 4)
# plt.plot(data9)
# plt.xlim(13200,17000)
# plt.ylim(0,180)
# plt.plot(color=colors, cmap='spring')
# plt.title('Timepoint 9')

# plt.subplot(10, 2, 5)
# plt.plot(data10)
# plt.xlim(13200,17000)
# plt.ylim(0,180)
# plt.plot(color=colors, cmap='spring')
# plt.title('Timepoint 10')

# plt.savefig( "TimeSeriesSubplots.png" )
# plt.close()


for arg in sys.argv[2:]:
	df = pd.read_table(arg, header=None)
	df.columns = ["Species", "Position", "Depth"]

	data = df[["Position", "Depth"]]
	fig, ax = plt.subplots(figsize=(12,6))
	ax.plot(data, c="r")
	ax.set_title( "Depth vs. Position for " + arg)
	ax.set_xlabel( "Position")
	plt.xlim(-1000,120000)
	ax.set_ylabel( "Depth" )
	plt.ylim(0,1500)
    # x1position = people[1]["start"]
#     for xc in x1position:
#         plt.axvline(x=xc, color='b', linestyle='-')
#     x2position = people[1]["end"]
#     for xc in x2position:
#         plt.axvline(x=xc, color='r', linestyle='--')
	# ax.autoscale_view()
	plt.tight_layout()
	continue

fig.savefig( "DepthvsPositionEcoliTimeSeries_loop.png" )
plt.close()


# gene_locations = {}
# gene_locations[]
# reference_genes = pd.read_table(sys.argv[1], header=None)

# for geneloc in reference_genes:
        
# print(reference_genes)

#     #dictionary gene 
    
# depth_file = pd.read_table(sys.argv[2])

#     # depth and number of bases covered
#     # in any given tissue, only about half
#     # coverage sum of ind depths / length of the
#     # some genes no genes 
    

