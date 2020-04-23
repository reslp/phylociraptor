#!/bin/bash

files=$(ls results/alignments_aligned/*_aligned)

for file in $files
do
	echo $file
	trimal -gappyout -in $file -out $file"_trimmed"
done
cd results
mkdir alignments_trimmed
mv alignments_aligned/*_trimmed alignments_trimmed/
