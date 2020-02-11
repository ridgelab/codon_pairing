#! /usr/bin/env python
#Makes a csv of number of codon pairs to the gene length
#input 1 - directory of pairing files
#input 2 - output csv

import sys
import os

outfile = open(sys.argv[2],"w")
outfile.write("Species,Num_Pairs,Gene_length\n")
for filename in os.listdir(sys.argv[1]):
    f = open(sys.argv[1] + "/" + filename,"r")
    numPairs = 0 #Number of codons in the motifs
    length = 0 #Number of codons in the gene
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
            if float(pairs) > 0.0: #As long as the codon appears > 0 times, we count it.
                numPairs += 1
            length += float(n)

    outfile.write(filename + "," + str(numPairs) + "," + str(length) + "\n")
    f.close()

outfile.close()

