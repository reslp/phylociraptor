#!/usr/bin/env python
# written by Philipp Resl

import os
import sys
import pandas as pd
from Bio import SeqIO
import argparse
import tarfile
from io import StringIO
from io import TextIOWrapper

if sys.version_info[0] < 3:
    raise Exception("Must be using Python 3")

pars = argparse.ArgumentParser(prog="create_sequence_files.py", description = """This script will create fasta files for all the buscos from all species with >% of buscos present""", epilog = """written by Philipp Resl""")
pars.add_argument('--busco_table', dest="busco_table", required=True, help="Path to BUSCO table.")
pars.add_argument('--busco_results', dest="busco_results", required=True, help="Results directory containing all BUSCO runs.")
pars.add_argument('--cutoff', dest="cutoff", default=0, required=True, help="Percent cutoff %% for BUSCO presence. Species below this threshold will be excluded.")
pars.add_argument('--outdir', dest="outdir", required=True, help="Path to output directory.")
pars.add_argument('--minsp', dest="minsp", required=True, help="Minimum number of species which have to be present to keep the sequences.")
pars.add_argument('--type' , dest="type", required=True, help="Type of sequences (aa or nu).")
pars.add_argument('--genome_statistics' , dest="genome_stats", required=True, help="Path to genome statistics output file.")
pars.add_argument('--gene_statistics' , dest="gene_stats", required=True, help="Path to gene statistics output file.")
pars.add_argument('--batchID' , dest="batchid", default=1, type=int, help="Batch ID (start for subsampling)")
pars.add_argument('--nbatches', dest="nbatches", default=1, type=int, help="Total number of batches (step size for subsampling)")

args=pars.parse_args()

extension=""
if args.type == "nu":
	extension = ".fna"
else:
	extension = ".faa" 


busco_overview = pd.read_csv(args.busco_table, sep="\t")
genomes = os.listdir(args.busco_results)
#print(busco_overview)
print("Settings:")
print("cutoff: ", args.cutoff)
print("minsp: ", args.minsp)
print("type: ", args.type)
print("outdir: ", args.outdir)
print("batchID: %i / %i" %(args.batchid, args.nbatches))

species_list = busco_overview.species.tolist()
print("Original number of species:", len(species_list))
#print(species_list)
#first remove species with too low busco coverage
busco_overview = busco_overview.set_index("species")

genome_file = open(args.genome_stats, "w")
for sp in species_list:
	if busco_overview.loc[sp, "percent_complete"] < float(args.cutoff):
		out = sp + "\tFAILED" + "\t" + str(busco_overview.loc[sp, "percent_complete"]) + "\t" + str(float(args.cutoff))
		print(out, file=genome_file)
		busco_overview = busco_overview.drop([sp])
	else:
		out = sp + "\tOK" + "\t" + str(busco_overview.loc[sp, "percent_complete"]) + "\t" + str(float(args.cutoff)) 
		print(out, file=genome_file)
species_list =  list(busco_overview.index)
print("Species remaining after applying cutoff:", len(species_list))
genome_file.close()

#now loop through each busco and extract sequence for each species
buscos = list(busco_overview.columns.values)
buscos.remove("percent_complete")

target=int(args.batchid)
gene_file = open(args.gene_stats, "w").close()
for i in range(len(buscos)):
	j = i+1
#	print("i: %i; j: %i; target: %i" %(i,j, target))
	if j != target:
#		print("Don't do anything (i: %i)" %i)
		continue
	gene_file = open(args.gene_stats, "a")
	target+=args.nbatches
	busco = buscos[i]
	print("Processing: " + busco)
	numseqs = 0
	outstring = ""
	for species in species_list:
		tar_file_content = open(args.busco_results + "/" + species + "/run_busco/single_copy_busco_sequences.txt", "r")
		for line in tar_file_content:
			line = line.strip()
			if busco+extension in line:
				path_to_busco_file = line.split(" ")[-1]	
				tf = tarfile.open(args.busco_results + "/" + species + "/run_busco/single_copy_busco_sequences.tar", "r")
				#print(path_to_busco_file)
				#print(tf.extractfile(path_to_busco_file))
				tar_file_content = TextIOWrapper(tf.extractfile(path_to_busco_file))
				if tar_file_content:
					with TextIOWrapper(tf.extractfile(path_to_busco_file)) as handle:
						for seq_record in SeqIO.parse(handle, "fasta"):
							name = ">" +species+"\n"
							#outfile.write(name)
							#outfile.write(str(seq_record.seq)+"\n")
							outstring = outstring + name
							outstring = outstring + str(seq_record.seq) + "\n"
				else:
					print(busco, "not found for", species)
					continue
	if outstring.count(">") >= int(args.minsp):	# only keep sequences if total number is larger than specified cutoff above.		
		print(busco + "\t" + "OK" + "\t" + str(outstring.count(">")) +"\t" + str(int(args.minsp)), file=gene_file)
		outfile = open(args.outdir+"/"+busco+"_all.fas", "w")
		outfile.write(outstring)
		outfile.close()
	else:
		print(busco + "\t" + "FAILED" + "\t" + str(outstring.count(">")) +"\t" + str(int(args.minsp)), file=gene_file)
	gene_file.close()
#old working code when busco sequences are not tared.
"""
	for species in species_list:
		for genome in genomes: # this loop is to get the correct directory name, it is very unelegant
			#print(args.busco_results+"/"+genome+"/run_busco/single_copy_busco_sequences/"+busco+extension)
			if species == genome:
				try:
					seqfile = open(args.busco_results + "/" + genome + "/run_busco/single_copy_busco_sequences/" + busco + extension, "r")
					for seq_record in SeqIO.parse(seqfile, "fasta"):
						name = ">" +species+"\n"
						#outfile.write(name)
						#outfile.write(str(seq_record.seq)+"\n")
						outstring = outstring + name
						outstring = outstring + str(seq_record.seq) + "\n"
					seqfile.close()
				except: # skip missing buscos
					continue
	if outstring.count(">") >= int(args.minsp):	# only keep sequences if total number is larger than specified cutoff above.		
		print(busco + "\t" + "OK" + "\t" + str(outstring.count(">")) +"\t" + str(int(args.minsp)), file=gene_file)
		outfile = open(args.outdir+"/"+busco+"_all.fas", "w")
		outfile.write(outstring)
		outfile.close()
	else:
		print(busco + "\t" + "FAILED" + "\t" + str(outstring.count(">")) +"\t" + str(int(args.minsp)), file=gene_file)
gene_file.close()
"""
