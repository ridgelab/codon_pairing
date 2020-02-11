#!/usr/bin/env python

#This script calculate the mean retention index and the standard deviation
#For each permutation.
#Input:
#1) input directory of state change files
#2) output file for summary file

import sys
import os
import numpy as np

allFile = open(sys.argv[3],"w")
allFile.write("Permutation number,meanRI,sdRI\n")

for filename in os.listdir(sys.argv[1]):
    infile = open(sys.argv[1] + "/" + filename,"r")
    n = filename.split("_")[-1]
    
    header = infile.readline()
    RIs = []
    for line in infile:
        fields = line.split("\t")
        if len(fields) != 8:
            continue
        gene,origin,loss,root_loss,n_small,n_total,percent_total,p = fields
        if float(n_small) == 1:
            continue
        g = float(n_small)
        s = float(origin) + float(loss) + float(root_loss)
        m = 1
        retention_index = float(g - s) / (g - m)
        RIs.append(retention_index)
    
    RIs = np.array(RIs)
    mean = np.mean(RIs)
    sd = np.std(RIs)
    allFile.write(n + "," + str(mean) + "," + str(sd) + "\n")
    infile.close()
