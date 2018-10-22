#! /usr/bin/env python

import sys
import math

inFile = open(sys.argv[1],"r")
outFile = open(sys.argv[2],"w")

totalGenes = 0
codonFreq = {}

for x in range(1,65):
	codonFreq[x] = 0
header = inFile.readline()
for line in inFile:
        line = line.strip()
        fields = line.split("\t")
	if len(fields) != 3:
                continue
        species = fields[0]
        num_genes = int(fields[2])
	totalGenes += num_genes
        codons = fields[1].replace("(","").replace(")","").replace(",","").replace("'","").split()
       	if len(codons) in codonFreq:
		codonFreq[len(codons)] += num_genes

outFile.write('"Codons_excluded","frequency"\n')
for entry in codonFreq:
        outFile.write(str(entry) + ',')
	logg = 0
	if codonFreq[entry] != 0:
		logg = math.log(codonFreq[entry])
        outFile.write(str(logg) + '\n')

inFile.close()
outFile.close()
