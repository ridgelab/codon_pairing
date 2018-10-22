#! /usr/bin/env python

import sys
import re

newick = ""
tree_line = ""
outfile = open(sys.argv[2],'w')
with open(sys.argv[1],'r') as tnt_output:
	for line_num, line in enumerate(tnt_output,start=1):
		if line.startswith("tree "):
			tree_line = line_num + 1
		if line_num==tree_line:
			newick = re.sub(":","",line)
			newick =re.sub(" ","",newick)
			newick = re.sub("_"," ",newick)
			#newick = re.sub(" (\d+) ",":\g<1> ",line)
			#newick = re.sub(" (\d+), ",":\g<1>, ",newick)
			#newick = re.sub(" (\d+)\)",":\g<1>)",newick)
outfile.write(newick)
outfile.close()
