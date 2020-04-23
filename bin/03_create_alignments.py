#!/usr/bin/env python2
# written by Philipp Resl
# this script will create fasta files for all the buscos from all species with >80% of buscos present

import os
import pandas as pd
from Bio import SeqIO

busco_overview = pd.read_csv("/data/scratch/philipp/projects/fungal_phylo/results/busco_table_all_fungi_processed.txt", sep="\t")
genomes = os.listdir("/data/scratch/philipp/projects/fungal_phylo/results/busco/")

species_list = busco_overview["species"].tolist()
print(len(species_list))
#print(species_list)
#first remove species with too low busco coverage
busco_overview = busco_overview.set_index("species")
for sp in species_list:
	if busco_overview.loc[sp, "percent_complete"] < 0.8:
		busco_overview = busco_overview.drop([sp])
		out = sp + "too few BUSCOs, will be removed"
		print out
species_list =  list(busco_overview.index)
print(len(species_list))

#now loop through each busco and extract sequence for each species
buscos = list(busco_overview.columns.values)
for busco in buscos:
	print("Processing: " + busco)
	outfile = open("results/alignments/"+busco+"_all.fas", "w")
	for species in species_list:
		for genome in genomes: # this loop is to get the correct directory name, it is very unelegant
			if species in genome:
				try:
					seqfile = open("results/busco/"+genome+"/single_copy_busco_sequences/"+busco+".faa", "r")
					for seq_record in SeqIO.parse(seqfile, "fasta"):
						name = ">" +species+"\n"
						outfile.write(name)
						outfile.write(str(seq_record.seq)+"\n")
				except: # skip missing buscos
					continue
	outfile.close()
