#! /usr/bin/env python

#This script calculates the R-squared coefficient of
#the number of codon pairs/number of possible codon pairs in a gene
#to the number of possible codon pairs/gene length (in codons)

#Input:
#1) input directory of pairing files
#2) output csv file

import sys
import os
import numpy as np
import re
from random import sample

outfile = open(sys.argv[2],"a")

allCodonPairing = [] #Stores numCodonPairs/numPossible for each codon in each gene of each species
allPossibleLength = [] #Stores numPossible/totalGeneLength for each codon in each gene of each species

allFiles = []
for filename in os.listdir(sys.argv[1]):
    allFiles.append(filename)

# Subset the data if necessary to decrease computational time/space
if len(allFiles) > 500:
    subset = sample(allFiles,500)
else:
    subset = allFiles

for filename in subset:
    infile = open(sys.argv[1] + "/" + filename, "r")
    gene_length = 0
    codons = {} #Dictionary of codon to its number of pairs and number of possible pairs
    for line in infile:
        if line[0] != ">":
            line = line.strip()
            codon,n_total,n_pairs = line.split("\t")
            #Make sure this is a normal codon
            if len(codon) == 3:
                check = re.sub(r"[ATCG]",r"",codon)
                if len(check) != 0:
                    continue
            n_total = float(n_total) #Total number of times the codon occurs in the gene
            n_possible = float(n_total) - 1 #Number of possible pairs
            n_pairs = float(n_pairs) #Number of actual pairs
            gene_length += float(n_total) #Add this to total number of codons in the gene
            if codon not in codons: 
                codons[codon] = [n_pairs,n_possible]
           
        else:
            if gene_length == 0:
                continue
            for c in codons:
                n_pairs = codons[c][0]
                n_possible = codons[c][1]
                if n_possible == 0:
                    continue
                allCodonPairing.append(float(n_pairs/n_possible))
                allPossibleLength.append(float(n_possible/gene_length))
                    
            gene_length = 0
            codons = {}


#Calculate r squared coefficient
degree = 1 #Indicates a linear model
coeffs = np.polyfit(allCodonPairing, allPossibleLength, degree)
correlation = np.corrcoef(allCodonPairing, allPossibleLength)[0,1]
r_square = correlation**2
print(r_square)
outfile.write(sys.argv[1] + "," + str(round(r_square,4)) + "\n")
outfile.close()
