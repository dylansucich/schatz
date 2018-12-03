#!/usr/bin/env python3

import os
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from scipy.cluster.hierarchy import dendrogram, linkage, leaves_list

directoryPath = "/Users/cmdb/schatz/RNAseq-tutorial/Exercise1/challenge2"
os.chdir(directoryPath)
folder_list=os.listdir(directoryPath)
# -5
completed = []
dirs = next(os.walk('.'))[1]
cufflinks = []
timepoints = pd.DataFrame()
for dir in dirs:
    if dir.endswith("clout"):
        cufflinks.append(dir)

print(cufflinks)


data = {}
col_labels = []
for i in cufflinks:
    print(i)
    time = i 
    numberoffiles = os.listdir(i)
    for file in os.listdir(i):
        if file.startswith("trans"):
            print(file)
            df = pd.read_table(file, index_col=0, header=None)
            print(df)
            df = df.loc[df[2] == "transcript"]
            field = df[8].str.split(";", expand=True)
            
            for i, col in enumerate(field.columns):
                field.iloc[:, i] = field.iloc[:, i].str.replace('"', '')
            # df = df[df. != ""]
            field = field.loc[field[1] != "exon_number"]
            # field.rename(col= ["gene_id", "transcript_id", "FPKM", "frac", "conf_lo", "conf_hi", "cov"])
            newfield = pd.DataFrame()
            for i in range(7):
                print(i)
                split = field[i].str.split(expand=True)[0]
                word = split.iloc[0]
                newfield[i] = field[i].str.split(expand=True)[1]
                if word not in col_labels:
                    col_labels.append(word)
                    
            newfield = newfield.iloc[:, 0:7]
            newfield.columns = col_labels
            df = df.iloc[:, 0:7]
            df = pd.concat([df, newfield], axis=1, sort=False)
            print(df)
            d = "bt21/"
            if not os.path.exists(d):
                os.makedirs(d)

            os.chdir(d)     
            out = str(time + ".out")
            df.to_csv(out, sep='\t', index=False, header=True)
            os.chdir("..")
            print(os.getcwd())
            data[time] = df.FPKM.tolist()
            # timepoints = pd.DataFrame()
            completed.append(out)
            
            
            # df[df[“column_z”] < 20].column_x
            # timepoints[time+str(i)] = df[df[“column_z”] < 20].column_x

print(data)
# for key, value in data:
#     print(key, value)
# print(timepoints)
            
# for i in range(len(field)):
#     # xaxis = col_names
#     data = df.T.iloc[-5,i]
#     plt.plot(data)
#     # plt.ylim(0,)
#     plt.xlabel('Experimental Time Point')
#     plt.ylabel('FPKM?')
#     plt.title("Time vs E.coli mRNA expression Cufflinks" + str(time))
#     continue
# plt.savefig("time_cufflinks_" + str(time) + ".png")
# plt.close()

        
        
# Folders=['01_load.TXT', '02_load.TXT', '03_load.TXT'] # These should be called filenames not folders but anyway.
#
# data_frames = {}    # Initialise a dictionary
#
# for filename in Folders:
#     df = pd.read_csv(filename, sep='\t', header=False)
#     data_frames[filename] = df
#
# # Now you can access any of the dataframes by the filename by using the dictionary:
# # Let's say you want the df associated with 02_load.TXT
#
# df = data_frames['02_load.TXT']
# print(df.head())
#
# trans = os.path.join(dir/"transcripts.gtf")
# print(trans)
#
# print(cufflinks)
#
# print("Path at terminal when executing this file")
# print(os.getcwd() + "\n")
#
# print("This file path, relative to os.getcwd()")
# print(__file__ + "\n")
#
# print("This file full path (following symlinks)")
# full_path = os.path.realpath(__file__)
# print(full_path + "\n")
#
# print("This file directory and name")
# path, filename = os.path.split(full_path)
# print(path + ' --> ' + filename + "\n")
#
# print("This file directory only")
# print(os.path.dirname(full_path))

        
        

