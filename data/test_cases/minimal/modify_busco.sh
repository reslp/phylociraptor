#!/bin/bash


echo "This script modifies the downloaded busco set so that it only includes 10 genes for use with the minimal test_case."

buscos="EOG092C5OAL EOG092C5OPO EOG092C5Q82 EOG092C5S4U EOG092C5T6H EOG092C5U93 EOG092C5UXM EOG092C5V3D EOG092C5Z2K EOG092C608K"

buscotmp="busco_tmp/fungi_odb9"

rm -rf $buscotmp

mkdir -p $buscotmp $buscotmp/info $buscotmp/hmms $buscotmp/prfl

buscodir="results/orthology/busco/busco_set/fungi_odb9"

if [[ ! -d $buscodir ]]; 
then
	echo "Directory $buscodir not found. Did you run phylociraptor setup and is this script in the phylociraptor base directory?"
	exit 1
fi

echo "Copy hmms and prfl files"

for gene in $buscos; do
	cp $buscodir/hmms/$gene.hmm $buscotmp/hmms
	cp $buscodir/prfl/$gene.prfl $buscotmp/prfl
done

echo "Reducing files in info directory"
zcat $buscodir/info/fungi_4751_OrthoDB9_orthogroup_info.txt.gz | head -n 1 > $buscotmp/info/fungi_4751_OrthoDB9_orthogroup_info.txt
for gene in $buscos; do
	zcat $buscodir/info/fungi_4751_OrthoDB9_orthogroup_info.txt | grep $gene
done >> $buscotmp/info/fungi_4751_OrthoDB9_orthogroup_info.txt
gzip $buscotmp/info/fungi_4751_OrthoDB9_orthogroup_info.txt

for gene in $buscos; do
	cat $buscodir/info/ogs.id.info | grep $gene
done > $buscotmp/info/ogs.id.info 

cp $buscodir/info/species.info $buscotmp/info
cp $buscodir/info/changelog $buscotmp/info

echo "Reducing lengths and score cutoff files"

for gene in $buscos; do
	cat $buscodir/lengths_cutoff | grep $gene
done > $buscotmp/lengths_cutoff

for gene in $buscos; do
	cat $buscodir/scores_cutoff | grep $gene
done > $buscotmp/scores_cutoff

echo "Copy ancestral files"

cp $buscodir/ancestral $buscotmp/
cp $buscodir/ancestral_variants $buscotmp/

echo "Modify dataset.cfg"

cat $buscodir/dataset.cfg | sed -e 's/number_of_BUSCOs=.*/number_of_BUSCOs=10/' > $buscotmp/dataset.cfg

echo "Replace busco directory with modified version"
echo "A Backup of the original busco directory will be kept in results/orthology/busco/busco_set_backup"

mv $buscodir results/orthology/busco/busco_set_backup
mv $buscotmp $buscodir

echo "done"
