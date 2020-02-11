#! /usr/bin/env python

#This script makes a tsv file of each codon, its state, and the species that use it.
#Input:
#1 - Input matrix file
#2 - Input character file
#3 - Output file

import sys

matrix = open(sys.argv[1],"r")
chars = open(sys.argv[2],"r")
out = open(sys.argv[3],"w")

one_species = {} #Dictionary of species that have the codon
                #Key is gene, value is list of species that have it
zero_species = {} #Dictionary of species that don't have the codon
characters = [] #Makes a list of the characters (codon in a gene)
for line in chars:
    line = line.strip()
    one_species[line] = []
    zero_species[line] = []
    characters.append(line)

chars.close()

for line in matrix:
    line = line.strip()
    fields = line.split("\t")
    if len(fields) != 2:
        continue
    species = fields[0]
    vals = fields[1]
    index = 0
    for v in vals:
        gene = characters[index]
        if v == "1":
            one_species[gene].append(species)
        elif v == "0":
            zero_species[gene].append(species)
        index += 1

matrix.close()

out.write("Gene\tState\tSpecies\n")
for c in characters:
    out.write(c + "\t0\t")
    s = ",".join(zero_species[c])
    out.write(s + "\n")
    
    out.write(c + "\t1\t")
    s = ",".join(one_species[c])
    out.write(s + "\n")

out.close()
