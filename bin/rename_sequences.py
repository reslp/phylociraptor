#!/usr/bin/env python
# written by Philipp Resl

import os
import sys
import glob
import argparse

if sys.version_info[0] < 3:
    raise Exception("Must be using Python 3")

pars = argparse.ArgumentParser(prog="rename_sequence_files.py", description = """This script will rename sequences so that they are prefixed by a species name""", epilog = """written by Philipp Resl""")
pars.add_argument('--filelist', dest="filelist", required=True, help="List of sequence files (comma separated).")
pars.add_argument('--prefixlist', dest="prefixlist", required=True, help="List of species names (comma separated)")
pars.add_argument('--outdir', dest="outdir", required=True, help="Output directory")

args=pars.parse_args()

files = args.filelist.split(",")
names = args.prefixlist.split(",")
for name, file in zip(names, files):
	with open(file, "r") as f:
		print("Renaming sequences for", name, "in", file)
		counter = 1
		with open(args.outdir + "/" + name + "_proteins.faa", "w") as outfile:
			for line in f:
				line = line.strip()
				if line.startswith(">"):
					line = ">" + line.split(">")[0] + name + "_" + str(counter).zfill(5)
					counter += 1
				print(line, file = outfile)



