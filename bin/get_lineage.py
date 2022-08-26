#!/usr/bin/env python3
import sys, time
from Bio import Entrez
import pandas as pd
if len(sys.argv) < 4:
	print("Not enough arguments!")
	print("usage:")
	print("get_lineage.py <download_statistics_file> <outputfile> <email>")
	sys.exit(1)

summary_file = open(sys.argv[1], "r")
outfilename = sys.argv[2]

Entrez.email = sys.argv[3]
handle = Entrez.efetch(db="Taxonomy", id=310279, retmode="xml")
records = Entrez.read(handle)
#import code
#code.interact(local=locals())
i=0
ntries = 10
all_lineage = {}
for line in summary_file:
	if line.startswith("name"):
		continue
	i += 1
	taxid = line.split("\t")[7]
	name = line.split("\t")[0]
#	if name not in ["Pseudois_nayaur", "Merops_nubicus", "Melospiza_melodia", "Gallirallus_okinawae", "Buphagus_erythrorhynchus"]:
#		continue
	print("Retrieving lineage for ",name, taxid)
	d = {}
	ntry = 0
	while ntry <= ntries:
		time.sleep(0.5)	
		handle = Entrez.efetch(db="Taxonomy", id=taxid, retmode="xml")
		records = Entrez.read(handle)
		if len(records) == 0:
			print("Try ", ntry, " to access NCBI database failed. Will try again")
			ntry += 1
			continue
		else:
			break
	if ntry == ntries and len(records) == 0:
		print("Failed to acccess NCBI for taxid: ",taxid, ". Giving up now")
		continue
	for info in records[0]["LineageEx"]:
		d[info["Rank"]] = info["ScientificName"]
	all_lineage[name] = d

df=pd.DataFrame(dict([ (k,pd.Series(v)) for k,v in all_lineage.items() ]))
df1 = df.T
df1.to_csv(outfilename, sep="\t",index=True, index_label="name")
#print(df1)
