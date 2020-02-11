#! /usr/bin/env python

import sys

inFile = open(sys.argv[1],"r")
outFile = open(sys.argv[2],"w")

aaDic = {"ATT":"I","ATC":"I","ATA":"I","CTT":"L","CTC":"L","CTA":"L","CTG":"L","TTA":"L","TTG":"L","GTT":"V","GTC":"V","GTA":"V","GTG":"V","TTT":"F","TTC":"F","ATG":"M","TGT":"C","TGC":"C","GCT":"A","GCC":"A","GCA":"A","GCG":"A","GGT":"G","GGC":"G","GGA":"G","GGG":"G","CCT":"P","CCC":"P","CCA":"P","CCG":"P","ACT":"T","ACC":"T","ACA":"T","ACG":"T","TCT":"S","TCC":"S","TCA":"S","TCG":"S","AGT":"S","AGC":"S","TAT":"Y","TAC":"Y","TGG":"W","CAA":"Q","CAG":"Q","AAT":"N","AAC":"N","CAT":"H","CAC":"H","GAA":"E","GAG":"E","GAT":"D","GAC":"D","AAA":"K","AAG":"K","CGT":"R","CGC":"R","CGA":"R","CGG":"R","AGA":"R","AGG":"R","TAA":"*","TAG":"*","TGA":"*"}

def makeCodonFreqSpecies():
	codonFreqSpecies = {'ATT': 0,'ATC': 0,'ATA': 0,'CTT': 0,'CTC': 0,'CTA': 0,'CTG': 0,'TTA': 0,'TTG': 0,'GTT': 0,'GTC': 0,'GTA': 0,'GTG': 0,'TTT': 0,'TTC': 0,'ATG': 0,'TGT': 0,'TGC': 0,'GCT': 0,'GCC': 0,'GCA': 0,'GCG': 0,'GGT': 0,'GGC': 0,'GGA': 0,'GGG': 0,'CCT': 0,'CCC': 0,'CCA': 0,'CCG': 0,'ACT': 0,'ACC': 0,'ACA': 0,'ACG': 0,'TCT': 0,'TCC': 0,'TCA': 0,'TCG': 0,'AGT': 0,'AGC': 0,'TAT': 0,'TAC': 0,'TGG': 0,'CAA': 0,'CAG': 0,'AAT': 0,'AAC': 0,'CAT': 0,'CAC': 0,'GAA': 0,'GAG': 0,'GAT': 0,'GAC': 0,'AAA': 0,'AAG': 0,'CGT': 0,'CGC': 0,'CGA': 0,'CGG': 0,'AGA': 0,'AGG': 0,'TAA': 0,'TAG': 0,'TGA': 0}
	return codonFreqSpecies

def updateCodonFreqClade(codonFreqSpecies,codonFreqClade):
	for entry in codonFreqSpecies:
		if entry in codonFreqClade:
			codonFreqClade[entry].append(codonFreqSpecies[entry] / float(totalGenes))
                else:
			list = []
			list.append(codonFreqSpecies[entry]/ float(totalGenes))
			codonFreqClade[entry] = list
	return codonFreqClade

codonFreqClade = {}
totalGenes = 0
header = inFile.readline()
currentSpecies = ""
totalSpecies = 0

for line in inFile:
	line = line.strip()
	fields = line.split("\t")
	if len(fields) != 3:
		continue
	species = fields[0]
	if species != currentSpecies:
		if currentSpecies != "":
			codonFreqClade = updateCodonFreqClade(codonFreqSpecies,codonFreqClade)

		codonFreqSpecies = makeCodonFreqSpecies()
		currentSpecies = species
		totalSpecies += 1
		totalGenes = 0
	num_unuse = int(fields[2])
	totalGenes += num_unuse
	codons = fields[1].replace("(","").replace(")","").replace(",","").replace("'","").split()
	for codon in codons:
		codonFreqSpecies[codon] += num_unuse

codonFreqClade = updateCodonFreqClade(codonFreqSpecies,codonFreqClade)

outFile.write('"codon","amino_acid","frequency"\n')
for i in range(totalSpecies):
	for entry in codonFreqClade:
		outFile.write('"' + entry + '","' + aaDic[entry] + '",')
		outFile.write('"' + str(codonFreqClade[entry][i]) + '"' + '\n')

inFile.close()
outFile.close()
		
