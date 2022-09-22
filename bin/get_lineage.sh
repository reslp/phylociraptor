#!/bin/bash
set -e

echo "get_lineage.sh: Provided data file: $1"
if [[ $2 != "keep" ]]; then
	echo "get_lineage.sh: Trying to download: ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz to data/lineage_information"
	rm -rf data/lineage_information
	mkdir -p data/lineage_information
	wget -q "ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz" -P data/lineage_information
	echo "get_lineage.sh: Download complete. Will unpack now."
	tar xvfz data/lineage_information/taxdump.tar.gz -C data/lineage_information
	echo "get_lineage.sh: Unpacking complete. Will reformat data with taxonkit:"
	echo "get_lineage.sh: Reformatting taxon data..."
fi
# this below is not yet properly implemented. Not sure how to handle this in the future:
#echo "get_lineage_information.sh: INFO: If unknown species taxa are named in the format Genus_sp the taxid of the genus will be used to increase the amount of available lineage information for unknown taxa."
#tail -n +2 $1 | awk -F"," '{print $1}' | sed 's/_/ /g' | sed 's/sp$//g' | taxonkit name2taxid --data-dir data/lineage_information | taxonkit lineage --data-dir data/lineage_information --taxid-field 2 -R > data/lineage_information/lineage_information.csv
tail -n +2 $1 | awk -F"," '{print $1}' | sed 's/_/ /g' | taxonkit name2taxid --data-dir data/lineage_information | taxonkit lineage --data-dir data/lineage_information --taxid-field 2 -R > data/lineage_information/lineage_information.csv
echo "get_lineage.sh: Reformating is done. NCBI taxnomy data is ready."
