#! /bin/bash

####################
#Run this script with the following arguments:
#1 - file of tnt trees
#2 Output file for the ete comparisons
#3 - Reference tree for the clade
####################

echo Reading from $1
echo "Ete Comparisons From $1 File" > $2

while read line; do  #Read through lines in the tree file
	if [[ $line == "("* ]];  #Check if the line is a tree
	then
		echo "tree tnt_1 = [&U]" > "tempFile.nex" #Put tnt tree tag in front of the tree
		echo $line >> "tempFile.nex"  #Store the tree in a temporary nexus tree file
		python convertNexusToNewick.py tempFile.nex newickTempFile.nwk #Convert it to newick
		ete3 compare --unrooted -t newickTempFile.nwk -r $3 >> $2  #Do our comparison and store in file
	fi
done <$1

#Clean up time
rm tempFile.nex
rm newickTempFile.nwk

echo "Average branch distance:" >> $2
python calcMeanComparison.py $2 >> $2 #calculate the mean of our comparisons
