#!/bin/sh
#This script is part of the phylo-script pipeline: github.com/reslp/phylo-scripts
#written by Philipp Resl 2015

###################################################################
# Provide the paths to your own files here:
IDFILE="results/IDS_for_tree.txt"
FASTAPATH="/data/scratch/philipp/projects/fungal_phylo/results/alignments" #path has to be absolute!!, all your single locus files need to be in that folder
FASTAFILES=$(ls /data/scratch/philipp/projects/fungal_phylo/results/alignments/*_aligned)
# no need to change anything below this point. But of course you can...
###################################################################

echo "Check if all files are present:"
for file in $FASTAFILES
do
 echo $file
done

echo "Renaming sequence names:"
for file in $FASTAFILES
do
	echo $file
	#sed "s/[.\']/_/g" $file > $file"_replaced"
done


echo "Looking for ID file"
[ -f ./$IDFILE ] && echo "...Found" || exit

for file in $FASTAFILES
do
	~/bin/phylo-scripts/reduce.py $IDFILE $file"_replaced" > $file"_replaced_reduced"
done

rm concat.fas

echo "Create concatenated alignment"
~/bin/phylo-scripts/concat.py $IDFILE results/alignments/*_reduced > concat.fas
echo "Done"
