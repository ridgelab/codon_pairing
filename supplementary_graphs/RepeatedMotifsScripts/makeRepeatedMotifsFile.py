#! /usr/bin/env python
import sys
import math

inFile = open(sys.argv[1],"r")
outFile = open(sys.argv[2],"w")

motifToFreq = {}
header = inFile.readline()

for line in inFile:
	fields = line.split("\t")
	if len(fields) != 3:
                continue
	species = fields[0]
	motif = fields[1].replace("(","").replace(")","").replace("'","").replace(",","").split()
	motif = str(motif)	
	if motif in motifToFreq:
		motifToFreq[motif] += int(fields[2])
	else:
		motifToFreq[motif] = int(fields[2])

repeatsToFreq = {}

for entry in motifToFreq:
	if motifToFreq[entry] > 2000:
		print entry
		print motifToFreq[entry]
	if motifToFreq[entry] in repeatsToFreq:
		repeatsToFreq[motifToFreq[entry]] += 1
	else:
		repeatsToFreq[motifToFreq[entry]] = 1

outFile.write('"motif_frequency","count"\n')

for entry in repeatsToFreq:
	logg = 0
        if repeatsToFreq[entry] != 0:
                logg = math.log(repeatsToFreq[entry])
	outFile.write(str(entry) + ',' + str(logg) + '\n')
