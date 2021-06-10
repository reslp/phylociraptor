#!/usr/bin/env python
# written by Philipp Resl

import os
import sys
import argparse


if sys.version_info[0] < 3:
    raise Exception("Must be using Python 3")

pars = argparse.ArgumentParser(prog="extract_busco_table.py", description = """This script will get all busco3 hmms and look in the busco results all specified genomes for the presence of the files""", epilog = """written by Philipp Resl""")
pars.add_argument('--hmm', dest="hmm_dir", required=True, help="Directory of the BUSCO hmms.")
pars.add_argument('--busco_results', dest="busco_results", required=True, help="Results directory containing all BUSCO runs.")
pars.add_argument('-o', dest="out", required=True, help="BUSCO table output file.")
args=pars.parse_args()

hmms = os.listdir(args.hmm_dir)
hmms = [hmm.strip(".hmm") for hmm in hmms]

genomes = os.listdir(args.busco_results)

outfile = open(args.out, "w")
header = "species\t"
header += "\t".join(hmms)
header += "\tpercent_complete" 
print(header, file= outfile)
for species in genomes:
	ones = 0
	zeros = 0
	outstring = species
	print("Extracting HMMs for", species, file=sys.stderr)
	try:
		busco_listing_file = open(args.busco_results + species + "/run_busco/single_copy_busco_sequences.txt", "r")
		buscos = []
		for line in busco_listing_file:
			line = line.strip()
			line = line.split(" ")[-1]
			line = line.split("/")[-1]
			if ".faa" in line: # only take aa files, this should be enough
				buscos.append(line)
		buscos = [busco.strip(".faa") for busco in buscos]
		busco_listing_file.close()
		for hmm in hmms:
			if hmm in buscos:
				outstring += "\t"
				outstring += "1"
				ones +=1
			else:
				outstring += "\t"
				outstring += "0"
				zeros +=1
		percent = ones / (ones+zeros)
		outstring += "\t"
		outstring += str(percent)
		print(outstring, file=outfile)
	except:
		out = species + " not found. Skipped.\n"
		print(out, file=sys.stderr)
		continue
outfile.close()
