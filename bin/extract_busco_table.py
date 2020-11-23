#!/usr/bin/env python
# written by Philipp Resl

import os
import sys
import argparse


if sys.version_info[0] < 3:
    raise Exception("Must be using Python 3")

pars = argparse.ArgumentParser(prog="get_busco_table.py", description = """This script will get all busco3 hmms and look in the busco results all specified genomes for the presence of the files""", epilog = """written by Philipp Resl""")
pars.add_argument('--hmm', dest="hmm_dir", required=True, help="Directory of the BUSCO hmms.")
pars.add_argument('--busco_results', dest="busco_results", required=True, help="Results directory containing all BUSCO runs.")
args=pars.parse_args()

hmms = os.listdir(args.hmm_dir)
hmms = [hmm.strip(".hmm") for hmm in hmms]

genomes = os.listdir(args.busco_results)

string = "species\t"
string += "\t".join(hmms)
string += "\tpercent_complete" 
print(string)
ones = 0
zeros = 0
for species in genomes:
	outstring = species
	print("Extracting HMMs for", species, file=sys.stderr)
	try:
		for hmm in hmms:
			buscos = os.listdir(args.busco_results + species + "/run_busco/single_copy_busco_sequences/")
			if hmm+".faa" in buscos:
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
		print(outstring)
	except:
		out = species + " not found. Skipped.\n"
		print(out, file=sys.stderr)
		continue
