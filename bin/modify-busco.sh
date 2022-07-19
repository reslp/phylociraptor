#!/bin/bash
set -e

echo "Hello from modify-busco.sh! This script modifies the downloaded BUSCO set to maintain only a smaller number of genes."

buscotmp=".busco_tmp/$2"
rm -rf $buscotmp
mkdir -p $buscotmp $buscotmp/info $buscotmp/hmms $buscotmp/prfl

buscodir="results/orthology/busco/busco_set/$2"
if [[ ! -d $buscodir ]]; then
	echo "The BUSCO directory does not exist: $buscodir"
	echo "Make sure you ran phylociraptor setup and that you specified the correct BUSCO set here."
	exit 1
fi


which=$(echo $1 | cut -d "=" -f 1)
genes=$(echo $1 | cut -d "=" -f 2)

if [[ $which == "ngenes" ]]; then
	echo "Will select $genes genes."
	if [[ $3 != "random" ]]; then
		seed=$3
	else
		seed=$RANDOM
	fi
	echo "Used random seed is: $seed"
	buscos="$(shuf --random-source=<(yes $seed) $buscodir/lengths_cutoff | awk '{print $1}' | head -n $genes | tr "\n" " ")"
elif [[ $which == "genes" ]]; then
	buscos=$(echo $genes | tr "," " ")
else
	echo "Something was misspecified for -n or -g"
	exit 1
fi

echo "Selected genes: "
echo $buscos

echo "Copy hmms and prfl files"

for gene in $buscos; do
	cp $buscodir/hmms/$gene.hmm $buscotmp/hmms
	cp $buscodir/prfl/$gene.prfl $buscotmp/prfl
done

echo "Reducing files in info directory"

for gene in $buscos; do
	cat $buscodir/info/ogs.id.info | grep $gene
done > $buscotmp/info/ogs.id.info 

cp $buscodir/info/species.info $buscotmp/info
if [[ -f $buscodir/info/changelog ]]; then cp $buscodir/info/changelog $buscotmp/info/; fi

echo "Reducing lengths and score cutoff files"

for gene in $buscos; do
	cat $buscodir/lengths_cutoff | grep $gene
done > $buscotmp/lengths_cutoff

for gene in $buscos; do
	cat $buscodir/scores_cutoff | grep $gene
done > $buscotmp/scores_cutoff

echo "Copy ancestral sequence files"

cp $buscodir/ancestral $buscotmp/
cp $buscodir/ancestral_variants $buscotmp/

echo "Modify dataset.cfg"

cat $buscodir/dataset.cfg | sed -e "s/number_of_BUSCOs=.*/number_of_BUSCOs=$(for i in $buscos; do echo $i; done | wc -l)/" > $buscotmp/dataset.cfg

echo "Replace busco directory with modified version"
echo "A Backup of the original BUSCO set directory will be kept in results/orthology/busco/busco_set/$2_backup"

mv $buscodir results/orthology/busco/busco_set/$2_backup
mv $buscotmp $buscodir

echo "Modifying BUSCO set done."
