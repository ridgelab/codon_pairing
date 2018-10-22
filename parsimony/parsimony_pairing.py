#! /usr/bin/env python

import sys
import argparse
from multiprocessing import Process, current_process, freeze_support, Pool
import re
import os


def makeAllPossibleCodons(rna):
	'''
	Input is a flag to specify if the sequence is DNA or RNA.
	Returns a set of all 64 possible codons.
	'''

	from itertools import product
	codons = product("ACGT",repeat=3)
	if rna:
		codons = product("ACGU",repeat=3)
	codonsComb = set()
	for c in codons:
		codonsComb.add("".join(c))
	return codonsComb
def parseArgs():
	'''
	Argument parsing is done.
	Required to have an input file.
	'''
	parser = argparse.ArgumentParser(description='Find Identical and co-tRNA codon pairing.')
	parser.add_argument("-t",help="Number of Cores",action="store",dest="threads",default=0,type=int, required=False)
	parser.add_argument("-i",help="Input Fasta Files",nargs='*',action="store", dest="input", required=False)
	parser.add_argument("-id",help="Input Directory with Fasta Files",action="store", dest="inputDir", required=False)
	parser.add_argument("-o",help="Output Matrix File",action="store",dest="outputMatrix", required=True)
	parser.add_argument("-oc",help="Output Characters File",action="store",dest="outputChars", required=False)
	parser.add_argument("-on",help="Output Species Names File",action="store",dest="outputNames", required=False)
	parser.add_argument("-f",help="Ribosome Footprint",action="store",dest="footprint", type=int, default=9, required=False)
	parser.add_argument("-c",help="Co-tRNA codon pairing",action="store_true",dest="co_trna", required=False)
	parser.add_argument("-comb",help="Combined co-tRNA and identical codon pairing",action="store_true",dest="comb", required=False)
	parser.add_argument("-rna",help="Flag for RNA sequences",action="store_true",dest="rna", required=False)
	args = parser.parse_args()

	if not args.input and not args.inputDir:
		print("You must supply an input file with either -i or -id")
		sys.exit()
	return args


def processOneSpecies(codonPairs):
	'''
	Finds the pairing for each species
	Requires a single input file as a parameter
	Returns a dictionary of the pairing, a list of all codons used, 
	and a dictionary containing the codon and a "0" or "1" for its use.
	'''
	dic = {}
	allCodonsSpecies = set()
	codonsToPairingSpecies = {}
	for header,pairs,unpairs in codonPairs:
		header = header.upper()
		geneName = header.split("GENE=")
		if len(geneName) == 1: #The gene doesn't have a name
			currentGene = ""
			continue
		currentGene = geneName[1].split("] ")[0]
		length = 64
		if args.co_trna or args.comb:
			length = 20
		dic[currentGene] = ['?'] * length
		for codon in pairs:
			if codon not in codonList and codon not in aminoList:
				continue
			geneCodon = currentGene + '__' + codon
			allCodonsSpecies.add(geneCodon)
			codonsToPairingSpecies[geneCodon] = '1'
			if args.co_trna or args.comb:
				index = aminoList.index(codon)
			else:
				index = codonList.index(codon)
			dic[currentGene][index] = '1'
		for codon in unpairs:
			if codon not in codonList and codon not in aminoList:
				continue
			geneCodon = currentGene + '__' + codon
			allCodonsSpecies.add(geneCodon)
			codonsToPairingSpecies[geneCodon] = '0'
			if args.co_trna or args.comb:
				index = aminoList.index(codon)
			else:
				index = codonList.index(codon)
			dic[currentGene][index] = '0'
	return dic,allCodonsSpecies,codonsToPairingSpecies


def getPairs(seq):
	footprint = args.footprint
	pairs = set()
	sequence = []
	codons = []
	if args.comb or args.co_trna:
		from Bio.Seq import Seq
		from Bio.Alphabet import generic_dna
		from Bio.Alphabet import generic_rna
		if args.rna:
			rna = Seq(seq,generic_rna)
			aa = str(rna.translate())
			sequence = re.findall(".",aa)
		else:
			seqaa = Seq(seq,generic_dna)
			aa = str(seqaa.translate())
			sequence =  re.findall(".",aa)
	if  args.co_trna:
		codons = re.findall("...",seq)
	if not args.co_trna and not args.comb:
		sequence = re.findall("...",seq)
	lastFound = dict() #key= codon or aa, value= position of last found codon with pairing
	for x in xrange(len(sequence)):
		curCodon = sequence[x]
		if not curCodon in lastFound or (x - lastFound[curCodon] >= footprint):
			lastFound[curCodon] =x
			continue
		if args.co_trna:
			if codons[x] == codons[lastFound[curCodon]]:
				continue
		pairs.add(curCodon)
		lastFound[curCodon] = x
	codonset = set(sequence)
	unpairs = codonset - pairs
	return pairs,unpairs


def readOneFile(inputFile):
	'''
	Reads one input file that is supplied as a parameter.
	'''
	input = ""
	header = ""
	sequence = ""
	codonPairs = []
	try:
		if inputFile[-3:] =='.gz':
			import gzip
			input = gzip.open(inputFile,'r')
		else:
			input = open(inputFile,'r')
		for line in input:
			if line[0] =='>':
				if sequence !="":
					pairs,unpairs = getPairs(sequence)
					codonPairs.append([header,pairs,unpairs])
				header = line
				sequence = ""
				continue
			sequence +=line.upper().strip()
		if sequence != "":
				pairs,unpairs = getPairs(sequence)
				codonPairs.append([header,pairs,unpairs])
				sequence = ""

	except Exception: #If the input file is malformatted, do not stop the program.
		input.close()
		print('Could not open file ' + str(inputFile))
		return {}
	input.close()


	dic,allCodonsSpecies,codonsToPairingSpecies = processOneSpecies(codonPairs)
	return dic,allCodonsSpecies,inputFile,codonsToPairingSpecies


def combineResults(results,fileNames):
	pairing = {}
	allCodons = set()
	informativeCodons = set()
	codonToFreq = {}
	codonsToPairings = {}
	for result in results:
		dic = result[0]
		allCodonsSpecies = result[1]
		species = result[2]
		codonsToPairingSpecies = result[3] 
		for codon in codonsToPairingSpecies:
			if codon in codonsToPairings:
				codonToFreq[codon] += 1
				if codonsToPairings[codon] != codonsToPairingSpecies[codon]:
					codonsToPairings[codon] = "P" #Found both a 0 and a 1, mark the codon as parsimony informative
				if codonsToPairings[codon] == "P" and codonToFreq[codon] == 4:  #Codon has 0 and 1, and also is found in minimum 4 species
					informativeCodons.add(codon)
			else:
				codonsToPairings[codon] = codonsToPairingSpecies[codon]
				codonToFreq[codon] = 1
		pairing[species] = dic
		allCodons = allCodons | allCodonsSpecies
	speciesToNumChars = {}
	speciesList,informativeCodons,pairing = removeCodons(fileNames,informativeCodons,pairing,codonToFreq,speciesToNumChars)
	return informativeCodons,pairing,speciesList

def removeCodons(speciesList,informativeCodons,pairing,codonToFreq,speciesToNumChars):
	'''
	Removes any codon that isn't found in at least 4 species.
	Removes any species that doesn't have at least 5% of the informative codons.
	Recurses until no more changes are made.
	Returns the list of species, list of informative codons, and pairing dictionary
	of species to codons.
	'''
	if len(informativeCodons) == 0:
		print("No informative codons were found.")
		sys.exit()
	if len(speciesToNumChars) == 0:
		for species in speciesList:
			num = 0
			i = 0
			j = 0
			for geneCodon in informativeCodons:
				fields = geneCodon.split('__')
				gene = fields[0]
				codon = fields[1]
				if codon not in codonList and codon not in aminoList:
					continue
				if args.co_trna or args.comb:
					index = aminoList.index(codon)
				else:
					index = codonList.index(codon)
				if gene in pairing[species]:
					if pairing[species][gene][index] != '?':
						num += 1
					else:
						i += 1
				else:
					j += 1
			speciesToNumChars[species] = num

	needsUpdate = False
	codonsToDelete = set()
	for species in speciesList:
		if speciesToNumChars[species] / float(len(informativeCodons)) < .05: 
			speciesList.remove(species)
			for geneCodon in informativeCodons:
				gene = geneCodon.split('__')[0]
				if gene in pairing[species]:
					if geneCodon in codonToFreq:
						codonToFreq[geneCodon] -= 1
						if codonToFreq[geneCodon] < 4:
							codonsToDelete.add(geneCodon)
							needsUpdate = True

	for geneCodon in codonsToDelete:
		informativeCodons.remove(geneCodon)
		fields = geneCodon.split('__')
		gene = fields[0]
		codon = fields[1]
		for species in speciesList:
			if gene in pairing[species]:
				if args.co_trna or args.comb:
					index = aminoList.index(codon)
				else:
					index = codonList.index(codon)
				if pairing[species][gene][index] != '?':
					speciesToNumChars[species] -= 1

	if needsUpdate:
		removeCodons(speciesList,informativeCodons,pairing,codonToFreq,speciesToNumChars)

	return speciesList,informativeCodons,pairing

def readInputFiles(args):
	'''
	Requires the system arguments to be passed to the function.
	'''
	threads = args.threads
	if threads ==0:
		pool = Pool()
	else:
		pool = Pool(threads)	
	allInputFiles = []
	allSets = set()
	fileToSet = {}
	if args.input:
		allInputFiles = args.input
	elif args.inputDir:
		allFasta = []
		path = args.inputDir
		allFasta = os.listdir(path)
		if path[-1] != '/':
			path += '/'
		allInputFiles = [path +i for i in allFasta]
	if len(allInputFiles) < 1:
		print("At least one input file is required")
		sys.exit()
	results = pool.map(readOneFile,allInputFiles)

	informativeCodons,pairing,speciesList = combineResults(results,allInputFiles)
	if args.outputChars:
		outfileChars = open(args.outputChars,"w")
		for codon in informativeCodons:
			outfileChars.write(codon + "\n")
		outfileChars.close()

	if args.outputNames:
		outfileNames = open(args.outputNames,"w")

	outfile = open(args.outputMatrix,"w")
	outfile.write("xread\n")
	outfile.write(str(len(informativeCodons)) + " " + str(len(speciesList)) + "\n")
	count = 0
	for species in speciesList:
		if args.outputNames:
			outfile.write("Species_" + str(count) + "\t")
			outfileNames.write("Species_" + str(count) + "\t" + species + "\n")
			count += 1
		else:
			speciesName = species.split("/")[1]
			outfile.write(speciesName + "\t")
		for item in informativeCodons:
			fields = item.split('__')
			gene = fields[0]
			codon = fields[1]
			if gene in pairing[species]:
				if args.co_trna or args.comb:
					index = aminoList.index(codon)
				else:
					index = codonList.index(codon)
				outfile.write(pairing[species][gene][index])
			else:
				outfile.write('?')
		outfile.write('\n')
	outfile.write(";\n")
	outfile.write("proc /;\ncomments 0\n;")
	outfile.close()

codonList = ['ATT','ATC','ATA','CTT','CTC','CTA','CTG','TTA','TTG','GTT','GTC','GTA','GTG','TTT','TTC','ATG','TGT','TGC','GCT','GCC','GCA','GCG','GGT','GGC','GGA','GGG','CCT','CCC','CCA','CCG','ACT','ACC','ACA','ACG','TCT','TCC','TCA','TCG','AGT','AGC','TAT','TAC','TGG','CAA','CAG','AAT','AAC','CAT','CAC','GAA','GAG','GAT','GAC','AAA','AAG','CGT','CGC','CGA','CGG','AGA','AGG','TAA','TAG','TGA']
aminoList = ['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V']


if __name__ =='__main__':
	'''
	Main.
	'''

	args = parseArgs()
	codonsComb = makeAllPossibleCodons(args.rna)
	readInputFiles(args)

