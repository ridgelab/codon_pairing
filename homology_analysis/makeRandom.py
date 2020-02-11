#! /usr/bin/env python

#This script creates random shuffled trees.
#Input is a reference tree and an output file for the shuffled trees

import sys
from random import shuffle
import re

inputF = open(sys.argv[1])
output = open(sys.argv[2],'w')

tree = inputF.readline()
tree = re.sub("\(","\n(\n",tree)
tree = re.sub("\)","\n)\n",tree)
tree = re.sub(",","\n,\n",tree)
tree = re.sub(";","\n;\n",tree)
tree = re.sub("\n\n","\n",tree)
lines = tree.split("\n")
lines = lines[1:-1]
species = []
badChar = ['(',')',',',';']
for line in lines:
	if not line.strip() in badChar:
		species.append(line.strip())

print(species)

for i in range(1000): #Change this number to alter the number of random permutations
    shuffle(species)
    curSpecies = 0
    tree = ""
    for line in lines:
    	if not line.strip() in badChar:
    		tree += species[curSpecies]
    		curSpecies +=1
    	else:
    		tree += line.strip()
    output.write(tree + "\n")

output.close()
inputF.close()

