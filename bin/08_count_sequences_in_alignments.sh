#!/bin/bash

files=$(ls ./results/alignments_trimmed/*trimmed)

for file in $files
do
	printf $file"\t"
	cat $file | grep ">" | wc -l
done
