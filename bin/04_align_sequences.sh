#!/bin/bash

files=$(ls /data/scratch/philipp/projects/fungal_phylo/results/alignments/*.fas)

cd results
cd alignments
for file in $files
do
	echo $file
	mafft --auto $file > $file"_aligned" 
done
cd ..
mkdir alignments_aligned
mv alignments/*_aligned alignments_aligned/
