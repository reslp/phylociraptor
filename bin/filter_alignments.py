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
pars.add_argument('--min-parsimony', nargs="?", dest="minpars", const=0, help="The minimum number of parsimony informative sites and alignment must have to be kept.")
pars.add_argument('--statistics-file', nargs="?", dest="statfile", const="", help="Statistics file which contains alignment information(incl. parsimony informative sites; this file needs to be created with concat)")
pars.add_argument('--minsp', dest="minsp", action="store", default=3, type=int, help="Number of taxa to be present in alignment at least (default: 3)")
args=pars.parse_args()

algn_list = args.align.split(" ")

stats_dict = {}
if args.statfile:
	with open(args.statfile, "r") as stats:
		stats.readline()
		for stat in stats:
			stats_dict[stat.split("\t")[0]] = int(stat.split("\t")[5])	
for al in algn_list:
	if os.stat(al).st_size == 0:
		print(os.path.basename(al), "\tEMPTY")
		continue
	if stats_dict and stats_dict[os.path.basename(al)] < int(args.minpars):
		print(os.path.basename(al), "\tNOT_INFORMATIVE")
		continue
	
	seqfile = open(al, "r")
	sequences = []
	for seq_record in SeqIO.parse(seqfile, "fasta"):
		sequences.append(seq_record)
	names_list = [seq.id for seq in sequences]
	if len(names_list) < args.minsp:
		print(os.path.basename(al), "\tTOO_FEW_SEQUENCES")
		continue
	if len(names_list) == len(set(names_list)):
		print("Create symlink: ", args.outdir+"/"+al.split("/")[-1], file=sys.stderr)
		os.symlink(al, args.outdir+"/"+al.split("/")[-1])
		print(os.path.basename(al), "\tPASS")
	else:
		print("Warning: File %s contains duplicated sequence IDs!" % al, file=sys.stderr)
		if (args.perseq == True):	
			print("Duplicated sequences will be removed and copy of file will be made.", file=sys.stderr)
			dups = [name for name in names_list if names_list.count(name)>1]
			print("Warning: Will remove %d sequences" % len(dups), file=sys.stderr)
			dups = set(dups)
			sequences = [sequence for sequence in sequences if sequence.id not in dups]
			if len(sequences) < args.minsp:
				print("Warning: %s file will be discarded (too few sequences after duplicate removal)" % al, file=sys.stderr)
				print(os.path.basename(al), "\tTOO_FEW_SEQUENCES")
				continue
			outfile = open(args.outdir+"/"+al.split("/")[-1], "w")
			for sequence in sequences:
				print(">"+sequence.id, file=outfile)
				print(sequence.seq, file=outfile)					
			outfile.close()
			print(os.path.basename(al), "\tREMOVE_DUPLICATED_SEQUENCES")
		else:
			print("Warning: %s file will be discarded" % al, file=sys.stderr)
			print(os.path.basename(al), "\tREMOVED")	


