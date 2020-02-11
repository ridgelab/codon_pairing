#! /usr/bin/env python

#This script makes a csv of number of codon pairs to the gene length.
#input 1 - directory of pairing files
#input 2 - output csv

import sys
import os

outfile = open(sys.argv[2],"w")
outfile.write("Species,Num_Pairs,Gene_length\n")
for filename in os.listdir(sys.argv[1]):
    f = open(sys.argv[1] + "/" + filename,"r")
    numPairs = 0 #Total number of codon pairs
    length = 0 # Total gene length
    gene = 0 
    for line in f:
        if line[0] == ">":
            if gene != 0:
                outfile.write(filename + "," + str(numPairs) + "," + str(length) + "\n")
                numPairs = 0
                length = 0
            gene += 1

        else:
            codon,n,pairs = line.split()
            numPairs += float(pairs)
            length += float(n)

    outfile.write(filename + "," + str(numPairs) + "," + str(length) + "\n")
    f.close()

outfile.close()

