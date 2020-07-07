#!/usr/bin/env python
# written by Philipp Resl

import os
import sys
import argparse
from Bio import SeqIO

if sys.version_info[0] < 3:
    raise Exception("Must be using Python 3")

pars = argparse.ArgumentParser(prog="filter_alignments.py", description = """This script will remove alignments with duplicated sequences.""", epilog = """written by Philipp Resl""")
pars.add_argument('--alignments', dest="align", required=True, help="alignment files")
pars.add_argument('--outdir', dest="outdir", required=True, help="output directory")
args=pars.parse_args()

algn_list = args.align.split(" ")

for al in algn_list:
	seqfile = open(al, "r")
	names_list = []
	for seq_record in SeqIO.parse(seqfile, "fasta"):
		names_list.append(seq_record.id)
	if len(names_list) == len(set(names_list)):
		print("Create symlink: ", args.outdir+"/"+al.split("/")[-1])
		os.symlink(al, args.outdir+"/"+al.split("/")[-1])
	else:
		print("File %s contains duplicated sequence IDs!" % al)


