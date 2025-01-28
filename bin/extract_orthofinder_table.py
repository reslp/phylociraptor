#!/usr/bin/env python
# written by Philipp Resl

import os
import sys
import argparse
import pandas as pd


if sys.version_info[0] < 3:
    raise Exception("Must be using Python 3")

pars = argparse.ArgumentParser(prog="extract_orthofinder_table.py", description = """This script will parse the Orthogroups.GeneCount.tsv file and transpose it. Optionally it will extract only orthogroups with one or zero genes per species.""", epilog = """written by Philipp Resl""")
pars.add_argument('--of', dest="of_counts", required=True, help="Path to Orthogroups.GeneCount.tsv file.")
pars.add_argument('-o', dest="out", required=True, help="Orthofinder table output file.")
pars.add_argument('--strict', dest="strict", action="store_true", default=False, help="Extract only orthogroups with one or zero number of copies in all taxa")
args=pars.parse_args()

data = pd.read_csv(args.of_counts, sep="\t", header=0)
data = data.drop("Total", axis=1)
data = data.set_index("Orthogroup")
data = data.transpose()

l = len(data.columns)
if args.strict:
	# check if column contains only zero (no ortholog) or one (single copy ortholog):
	allcols = data.columns
	#print(data.head())
	for col in allcols:
		print("Orthogroups matching criteria:", j, "Checking column", i, "of", l, "columns:", col, end="\r") 
		if data[col].isin([0,1]).all() != True:
			data = data.drop(col, axis=1)
		else:
			j += 1
		i += 1

data["percent_complete"] = data[data.columns].sum(axis=1) / l # calculate completeness proportion used for filtering later.
newheader= ["species"] + list(data.columns)
outstring = data.to_csv(index=True, header=False, sep='\t')

with open(args.out, "w") as outfile:
	print("\t".join(newheader), file=outfile)
	print(outstring, file=outfile)
