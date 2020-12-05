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
pars.add_argument('--per_sample', dest="perseq", action='store_true', help="if set, only samples with duplicated sequences will be removed. Else the whole alignment")
args=pars.parse_args()

algn_list = args.align.split(" ")

for al in algn_list:
	seqfile = open(al, "r")
	sequences = []
	for seq_record in SeqIO.parse(seqfile, "fasta"):
		sequences.append(seq_record)
	names_list = [seq.id for seq in sequences]
	if len(names_list) == len(set(names_list)):
		print("Create symlink: ", args.outdir+"/"+al.split("/")[-1])
		os.symlink(al, args.outdir+"/"+al.split("/")[-1])
	else:
		print("Warning: File %s contains duplicated sequence IDs!" % al)
		if (args.perseq == True):	
			print("Duplicated sequences will be removed and copy of file will be made.")
			dups = [name for name in names_list if names_list.count(name)>1]
			print("Warning: Will remove %d sequences" % len(dups))
			dups = set(dups)
			sequences = [sequence for sequence in sequences if sequence.id not in dups]
			outfile = open(args.outdir+"/"+al.split("/")[-1], "w")
			for sequence in sequences:
				print(">"+sequence.id, file=outfile)
				print(sequence.seq, file=outfile)					
			outfile.close()
		else:
			print("Warning: %s file will be discarded" % al)	


