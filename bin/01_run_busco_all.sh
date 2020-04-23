#!/bin/bash

cd results
mkdir busco
cd busco
files=$(ls ../../data/fungi_incl_some_lichens/*.gz)
for file in $files
do
	gunzip $file
	echo "Preparing to run BUSCO on: "
	echo $file
	uncompressed_file=$(echo $file | sed -e "s/\.gz//g")
	outname=$(echo $uncompressed_file | sed -e "s/\.fna//g")
	#outname=$(echo $outname | sed -e "s/\.\. //g")
	outname=$(basename $outname)
	echo $outname
	echo $file
	python /usr/local/src/busco3/scripts/run_BUSCO.py -i $uncompressed_file -o busco_$outname -l /usr/local/src/busco3/fungi_odb9 -m genome -c 8
	gzip $uncompressed_file
done
