#!/usr/bin/env python3

#    compute average depth per exon from the reference gene list
#      DepthperExon.py refgenes.ptt t#.

import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
plt.style.use('ggplot')



for filename in sys.argv[1:]:
	df = pd.read_table(sys.argv[sys.argv[2]], header=None)
	df.columns = ["Species", "Position", "Depth"]

	data = df[["Position", "Depth"]]
	fig, ax = plt.subplots(figsize=(12,6))
	ax.plot(data, c="r")
	ax.set_title( "Depth vs. Position for " + sys.argv[filename])
	ax.set_xlabel( "Position")
	plt.xlim(-1000,120000)
	ax.set_ylabel( "Depth" )
	plt.ylim(0,1500)
	# ax.autoscale_view()
	plt.tight_layout()
	continue

fig.savefig( "DepthvsPositionEcoliTimeSeries.png" )
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
    

