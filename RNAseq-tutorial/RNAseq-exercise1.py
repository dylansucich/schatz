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

genes = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,
34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,
69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100]

df = pd.read_csv(sys.argv[1], sep = "\t", index_col=0)

for i in genes:
    genelist = df.loc["gene_" + str(i)][:]
    xaxis = ["exp_1",  "exp_2",  "exp_3",  "exp_4",  "exp_5",  "exp_6",  "exp_7",  "exp_8",  "exp_9",  "exp_10"]
    yaxis = genelist
    plt.figure()
    plt.plot(xaxis, yaxis)
    plt.title("gene_" + str(i) +".png")
    plt.ylim(0, 110)
    plt.savefig("gene_" + str(i) + ".png")
    plt.close()
    continue

