#! /usr/bin/env python
'''
Takes a directory of fasta files (or a specified list of fasta files),
calculates codon pairings within each sequence. Adds a tuple of each pair
to a set for a species (the file name), and calculates the distance between
that species and all other species analyzed. Creates a distance matrix to be
used by makeNewick.py (or any other distance algorithm).
'''
import sys
import argparse
from multiprocessing import Process, current_process, freeze_support, Pool
import re
import os


def makeAllPossibleCodons(rna,co_trna):
	'''
	Input: rna is a flag to specify if the sequence is DNA or RNA. co_trna is a flag to
		specify co_tRNA codon pairing (i.e., same amino acid formed)
	Returns a set of all 64 possible codons (DNA or RNA) or all 20 amino acids.
	'''
	if co_trna:
		return set(['A','R','N','D','B','C','E','Q','Z','G','H','I','L','K','M','F','P','S','T','W','Y','V'])
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
	parser.add_argument("-o",help="Output Directory",action="store",dest="output", required=False)
	parser.add_argument("-f",help="Ribosome Footprint",action="store",dest="footprint", type=int, default=9, required=False)
	parser.add_argument("-c",help="Co-tRNA codon pairing",action="store_true",dest="co_trna", required=False)
	parser.add_argument("-rna",help="Flag for RNA sequences",action="store_true",dest="rna", required=False)
	args = parser.parse_args()

	if not args.input and not args.inputDir:
		print "You must supply an input file with either -i or -id"
		sys.exit()
	return args


def getPairs(seq,orderedCodons):
	footprint = args.footprint
	pairs = dict()
	codons = []
	if args.co_trna:
		from Bio.Seq import Seq
		from Bio.Alphabet import generic_dna
		from Bio.Alphabet import generic_rna
		if args.rna:
			rna = Seq(seq,generic_rna)
			aa = str(rna.translate())
			codons = re.findall(".",aa)
		else:
			seq = Seq(seq,generic_dna)
			aa = str(seq.translate())
			codons = re.findall(".",aa)
	else:
		codons = re.findall("...",seq)
	lastFound = dict() #key= codon, value= position of last found codon with pairing
	foundPairing = set()
	for x in xrange(len(codons)):
		curCodon = codons[x]
		if not curCodon in orderedCodons:
			continue
		if not curCodon in lastFound or (x - lastFound[curCodon] >= footprint):
			lastFound[curCodon] =x
			continue
		foundPairing.add(curCodon)
		lastFound[curCodon] = x
	return tuple(sorted(list(foundPairing)))


def readOneFile(inputFile):
	'''
	Reads one input file that is supplied as a parameter.
	Returns a set of codon pairing tuples from a species and the accompanying path to that species file.
	'''
	input = ""
	header = ""
	sequence = ""
	codonPairs = set()
	orderedCodons = makeAllPossibleCodons(args.rna, args.co_trna)
	try:
		if inputFile[-3:] =='.gz':
			import gzip
			input = gzip.open(inputFile,'r')
		else:
			input = open(inputFile,'r')
		for line in input:
			if line[0] =='>':
				if sequence !="":
					codonPairs.add(getPairs(sequence,orderedCodons))
				header = line
				sequence = ""
				continue
			sequence +=line.upper().strip()
		if sequence != "":
			codonPairs.add(getPairs(sequence,orderedCodons))
			sequence = ""
	except Exception: #If the input file is malformatted, do not stop the program.
		input.close()
		return set()

	input.close()
	return (codonPairs, inputFile.split("/")[-1].split(".gz")[0])
	

def readInputFiles(args):
	'''
	Requires the system arguments to be passed to the function.
	Writes the distances between each species baseed on codon pairing to an output file or standard out.
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
		print "At least one input file is required"
		sys.exit()
	n = pool.map(readOneFile,allInputFiles)
	output = sys.stdout
	if args.output:
		output = open(args.output,'w')
	species = {}
	for normalized in n:
		if len(normalized) ==0:
			continue
		species[normalized[1]] =normalized[0]
	output.write("," + ",".join(species.keys()) + "\n")
	for s1,v1 in species.items():
		output.write(s1)
		e = []
		for s2,v2 in species.items():
			combined = v1 & v2
			e.append(len(combined))
			distance = 1-(float(len(combined)) / max(len(v1),len(v2)))
			output.write("," + str(distance))
		output.write("\n")
	if args.output:
		output.close()

if __name__ =='__main__':
	'''
	Main.
	'''

	args = parseArgs()
	readInputFiles(args)

