#! /usr/bin/env python
import sys
import math

inFile = open(sys.argv[1],"r")
outFile = open(sys.argv[2],"w")

motifToSpecies = {}
header = inFile.readline()

for line in inFile:
	fields = line.split("\t")
	if len(fields) != 3:
                continue
	species = fields[0]
	motif = fields[1].replace("(","").replace(")","").replace("'","").replace(",","").split()
	motif = str(motif)	
	if motif in motifToSpecies:
		motifToSpecies[motif].append(species)
	else:
		myList = []
		myList.append(species)
		motifToSpecies[motif] = myList 

outFile.write('"motifs"\n')
uniqueMotifs = 0
repeatedMotifs = 0
for entry in motifToSpecies:
	if len(motifToSpecies[entry]) > 1:
<<<<<<< HEAD
		repeatedMotifs += 1
	else:
=======
#		print "Motif is in " + str(len(motifToSpecies[entry])) + " species"
		outFile.write("Overlaps\n")
		repeatedMotifs += 1
	else:
		outFile.write("Unique\n")
>>>>>>> 1cb174ba31a54f513c839296f90651f9dc93da67
		uniqueMotifs += 1

outFile.write("Unique Motifs: "+ str(uniqueMotifs) + "\n")
outFile.write("Repeated Motifs: " + str(repeatedMotifs) + "\n")
		
