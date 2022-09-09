#!/usr/bin/env python3
import sys, time
from Bio import Entrez
import pandas as pd
if len(sys.argv) < 3:
	print("Not enough arguments!")
	print("usage:")
	print("reformat_lineage.py <lineage_file> <outputfile>")
	sys.exit(1)

summary_file = open(sys.argv[1], "r")
outfilename = sys.argv[2]

all_lineages = {}
print("Reformatting all linegaes:")
for line in summary_file:
	line = line.rstrip()
	if len(line.split("\t")) == 4:
		sp = line.split("\t")[0].replace(" ", "_")
		taxid = line.split("\t")[1]
		lineage = line.split("\t")[2].split(";")
		ranks = line.split("\t")[3].split(";")
		d = {}
		for item in zip(lineage, ranks):
			d[item[1]] = item[0]
		all_lineages[sp] = d

print("Creating output CSV file:", outfilename)
df=pd.DataFrame(dict([ (k,pd.Series(v)) for k,v in all_lineages.items() ]))
df1 = df.T
df1.to_csv(outfilename, sep="\t",index=True, index_label="name")
