#!/usr/bin/env python

#This script calculates the actual retention index 
#for each character (codon)
#The summary file is for the mean retention index of 
#all the codons
#Input:
#1 - input file of codon score table
#2 - output file
#3 - output summary file

import sys
import numpy as np

infile = open(sys.argv[1],"r")
outfile = open(sys.argv[2], "w")
summaryFile = open(sys.argv[3],"a")

header = infile.readline()
RIs = [] #List of retention indices
for line in infile:
    gene,origin,loss,root_loss,n_small,n_total,percent_total,p = line.split("\t")
    if float(n_small) == 1:
        continue
    g = float(n_small)
    s = float(origin) + float(loss) + float(root_loss)
    m = 1
    retention_index = float(g - s) / (g - m)
    RIs.append(retention_index)
    outfile.write(gene + "\t" + str(retention_index) + "\n")

RIs = np.array(RIs)
mean = np.mean(RIs)
sd = np.std(RIs)
print(mean)
print(sd)
summaryFile.write(sys.argv[1] + "," +  str(round(mean,5)) + "," + str(round(sd,5)) + "\n")

infile.close()
outfile.close()
summaryFile.close()
