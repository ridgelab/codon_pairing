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
	parser.add_argument("-o",help="Output File",action="store",dest="output", required=False)
	parser.add_argument("-f",help="Ribosome Footprint",action="store",dest="footprint", type=int, default=9, required=False)
	parser.add_argument("-c",help="Co-tRNA codon pairing",action="store_true",dest="co_trna", required=False)
	parser.add_argument("-b",help="Both Identical and Co-tRNA codon pairing",action="store_true",dest="both", required=False)
	parser.add_argument("-rna",help="Flag for RNA sequences",action="store_true",dest="rna", required=False)
	parser.add_argument("-l",type=str, help="Codon Table. Default: Standard",action="store",dest="codon_table", default="Standard", required=False)
	args = parser.parse_args()

	if not args.input and not args.inputDir:
		sys.stdout.write("You must supply an input file with either -i or -id\n")
		sys.exit()
	if args.co_trna and args.both:
		sys.stdout.write("You cannot use both the co_trna (-c) and both (-b) flags.\n")
		sys.exit()

	return args


def getPairs(seq,orderedCodons):
	footprint = args.footprint
	pairs = dict()
	codons = []
	dna_codons = []
	if args.both or args.co_trna:
		from Bio.Seq import Seq
		from Bio.Alphabet import generic_dna
		from Bio.Alphabet import generic_rna
		if args.rna:
			rna = Seq(seq,generic_rna)
			aa = str(rna.translate(table=args.codon_table))
			codons = re.findall(".",aa)
		else:
			sequence = Seq(seq,generic_dna)
			aa = str(sequence.translate(table=args.codon_table))
			codons = re.findall(".",aa)
	else:
		codons = re.findall("...",seq)
	if args.co_trna: #To ensure that identical codon pairing does not form the amino acid
		dna_codons = re.findall("...",seq)

	lastFound = dict() #key= codon, value= position of last found codon with pairing #For co-trna: key = codon (amino acid), value= dict() where key=dna_codon (codon) and value = last position of it
	foundPairing = set()
	for x in range(len(codons)):
		curCodon = codons[x]
		if not curCodon in orderedCodons:
			continue
		if not args.co_trna:
			if not curCodon in lastFound or (x - lastFound[curCodon] >= footprint): #Must be >= because if footprint is 2 and AAA is found at positions 3 and 4, 4-3 =1, which is 1 less than footprint size.
				lastFound[curCodon] =x
				continue
		else:
			if not curCodon in lastFound:  
				lastFound[curCodon] = dict()
				lastFound[curCodon][dna_codons[x]] =x
				continue
			closestPos = -100
			for key,value in lastFound[curCodon].items():
				if key == dna_codons[x]:
					continue
				if value >closestPos:
					closestPos = value
			if (x - closestPos) >= footprint:
				lastFound[curCodon][dna_codons[x]] =x
				continue

		foundPairing.add(curCodon)
		if args.co_trna:
			lastFound[curCodon][dna_codons[x]] =x
			continue
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
	orderedCodons = makeAllPossibleCodons(args.rna, (args.co_trna | args.both))
	gzipped = False
	try:
		if inputFile.endswith('.gz') or inputFile.endswith(".gzip"):
			gzipped =True
			import gzip
			input = gzip.open(inputFile,'r')
		else:
			input = open(inputFile,'r')
		for line in input:
			if gzipped:
				line = line.decode()
			if line[0] =='>':
				if sequence !="" and (len(sequence)%3==0):
					codonPairs.add(getPairs(sequence,orderedCodons))
				header = line
				sequence = ""
				continue
			sequence +=line.upper().strip()
		if sequence != "" and (len(sequence)%3==0):
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
		sys.stdout.write("At least one input file is required\n")
		sys.exit()
	n = pool.map(readOneFile,allInputFiles)
	species = {}
	for normalized in n:
		if len(normalized) ==0:
			continue
		species[normalized[1]] =normalized[0]
	if len(species) == 0:
		sys.stderr.write("No distances were calculated\n")
		sys.exit()
	output = sys.stdout
	if args.output:
		output = open(args.output,'w')
	sortedSpecies = sorted(species.keys())
	output.write("," + ",".join(sortedSpecies) + "\n")
	for s1 in sortedSpecies:
		v1 =species[s1]
		output.write(s1)
		for s2 in sortedSpecies:
			v2=species[s2]
			combined = v1 & v2
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
	if args.codon_table:
		try:
			if args.both or args.co_trna:
				from Bio.Seq import Seq
				from Bio.Alphabet import generic_dna
				from Bio.Alphabet import generic_rna
				if args.rna:
					rna = Seq("AUG",generic_rna)
					aa = str(rna.translate(args.codon_table))
				else:
					sequence = Seq("ATG",generic_dna)
					aa = str(sequence.translate(table=args.codon_table))
		except Exception:
			sys.stderr.write("Codon Table Does Not Exist.\n")
			sys.exit()


	readInputFiles(args)

