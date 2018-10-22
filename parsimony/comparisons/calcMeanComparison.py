#! /usr/bin/env python
import sys

infile = open(sys.argv[1])

total = 0
numTrees = 0
for line in infile:
	if line[0] == "(":
		fields = line.split()
		src = float(fields[12])
		ref = float(fields[14])
		numTrees += 1
		if src > ref:
			total += src
		else:
			total += ref

average = float(total)/numTrees
print average
		
