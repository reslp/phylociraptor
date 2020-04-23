#!/usr/bin/env python2
# written by Philipp Resl
# This script will get all ascomycota busco3 hmms and look in the busco results of 700+ fungal genomes for the presence of the files

import os
import sys
hmms = os.listdir("/usr/local/src/busco3/fungi_odb9/hmms/")
hmms = [hmm.strip(".hmm") for hmm in hmms]

genomes = os.listdir("/data/scratch/philipp/projects/fungal_phylo/results/busco/")

string = "species" + "\t".join(hmms)
print string
for species in genomes:
	outstring = species +"\t"
	try:
		for hmm in hmms:
			buscos = os.listdir("/data/scratch/philipp/projects/fungal_phylo/results/busco/" + species + "/single_copy_busco_sequences/")
			if hmm+".fna" in buscos:
				outstring += "\t"
				outstring += "1"
			else:
				outstring += "\t"
				outstring += "0"
		print outstring
	except:
		out = species + " not found\n"
		sys.stderr.write(out)
		continue
