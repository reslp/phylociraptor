#!/bin/bash

files=$(ls results/alignments/*_trimmed)
cd results
mkdir iqtree_single_gene_trees
cd iqtree_single_gene_trees
for file in $files
do
	cp ../../$file .
	filehere=$(ls *_trimmed)
	#iqtree -s $filehere -m MFP -mset WAG -msub nuclear -cmax 5 -bb 1000 -bnni -nt 8
	iqtree -s $filehere -m Dayhoff -bb 1000 -bnni -nt 8
	rm *_trimmed
done
