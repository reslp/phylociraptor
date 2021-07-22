#!/usr/bin/env python

import sys, os, time
import pandas as pd
from Bio import Entrez
import subprocess
import argparse

pars = argparse.ArgumentParser(prog="genome_download.py", description = """This script will download complete genomes from NCBI genbank""", epilog = """written by Philipp Resl""")
pars.add_argument('--genomes', dest="genomes", required=True, help="List of genomes to be downloaded sperated by commas (,). Requried.")
pars.add_argument('--outdir', dest="outdir",default="", help="output directory. Default: current directory")
pars.add_argument('--entrez_email', dest="email", required=True, help="Email to be used for communication with the NCBI Entrez databases. Required.")
pars.add_argument('--api_key', dest="api_key", help="API key to be used for communication with the NCBI Entrez databases. Optional; not yet implemented.")
args=pars.parse_args()

def now():
	return time.strftime("%Y-%m-%d %H:%M") + " -"


if not args.outdir.endswith("/"):
	args.outdir = args.outdir + "/"
wd = os.getcwd() + "/"
#args.outdir = wd + args.outdir
print(now(), "Using ouput directory:", args.outdir)

if not os.path.isfile(args.outdir+"assembly_summary_genbank.txt"):
	print(now(), "Downloading genome overview...")	
	try:
		outstream = subprocess.run(["wget","-q", "ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_genbank.txt", "-P", args.outdir], stderr=subprocess.STDOUT, stdout = subprocess.PIPE)	
		if outstream.returncode == 0:
			print(now(), "Download finished successfully.")
			output = ""
			with open(args.outdir+"assembly_summary_genbank.txt","r") as handle: # this needs to be done because this file contains two comment head lines.
				handle.readline() # read first line so that pointer is at second line
				for line in handle:
					if line.startswith("#"):
						line = line[2:]
					output += line
			outfile = open(args.outdir+"assembly_summary_genbank.txt","w")
			outfile.write(output) 
			outfile.close()
		else:
			print(now(), "Something went wrong during download")
			sys.exit(1)
	except:
		print(now(), "An error occurred during download")
		sys.exit(1)
else:
	print(now(), "NCBI genome database file is already present. Will not download")

# parse command line arguments
genomes = args.genomes
genomes = genomes.split(",")
Entrez.email = args.email

# read NCBI genome database csv
print(now(), "Reading NCBI genome database file")
data = pd.read_csv(args.outdir + "assembly_summary_genbank.txt", sep="\t", error_bad_lines=False, low_memory=False, na_values="")
data["ftp_path"] = data["ftp_path"].astype(str)


#def check_taxonomy(sp, taxid):
#	print(now(), "Will test if", sp, "has taxid", taxid)
#	entrez_handle = Entrez.esummary(db="taxonomy", id=str(taxid))
#	#entrez_handle = Entrez.esearch(db="Taxonomy", term=sp)
#	record = Entrez.read(entrez_handle)
#	entrez_handle.close()
#	print(record)
#	if record[0]["ScientificName"] == sp:
#		print(now(),  taxid, "seems to be correct for", sp)
#		return True
#	else:
#		print(now(),  taxid, "seems to be wrong for", sp)
#		print(now(), "You can check manually if", sp, "=",record[0]["ScientificName"], ".")
#		return False	

def get_taxid(sp):
	entrez_handle = Entrez.esearch(db="Taxonomy", term=sp)
	record = Entrez.read(entrez_handle)
	entrez_handle.close()
	if len(record["IdList"]) == 0:
		return
	if record["IdList"][0]:
		return str(record["IdList"][0])


def check_already_downloaded(file):
	if os.path.isfile(file):
		return True
	else:
		return False

def download(download_data, genome):
	url = download_data.iloc[0]["ftp_path"]
	genome = genome.replace(" ", "_")
	filename = url.split("/")[-1]
	download_url_genome = url + "/" + filename + "_genomic.fna.gz"
	print(now(), "Will download genome from:", download_url_genome)
	if check_already_downloaded(args.outdir + genome + "_genomic_genbank.fna.gz"):
		print(now(), "A genome has already been downloaded. Will not download again.")
		return "success"
	else:
		try:
			outstream = subprocess.run(["wget", download_url_genome, "-P", args.outdir], stderr=subprocess.STDOUT, stdout=subprocess.PIPE)
			if outstream.returncode == 0:
				print(now(), "Download finished successfully.")
				outstream = subprocess.run(["mv", args.outdir + filename + "_genomic.fna.gz", args.outdir + genome + "_genomic_genbank.fna.gz"])
				print(now(), "Writing meta information file.")
				download_data.to_csv(args.outdir + genome +"_db_genbank.tsv", sep="\t")
			else:
				print(now(), "Something went wrong during download")
				return "failed"
		except:
			print(now(), "An error occurred during download")
			return "failed"
	print(now(), "done")
	return "success"

overview = {}
for genome in genomes:
	print("\n----------------------------------", genome.split("=")[0], "----------------------------------")
	if check_already_downloaded(args.outdir+genome+"_genomic_genbank.fna.gz"):
		print(now(), "A genome has already been downloaded. Will skip")
		overview[genome] = "success"
		continue

	#if specific accession is provided via 'web=accession'
	if "=" in genome:
		accession = genome.split("=")[1]
		print(now(), "A specific accession '"+accession+"' has been provided")
		if data.assembly_accession.str.fullmatch(accession).any():
			species_data = data[data.assembly_accession == accession]
		else:
			print(now(), "The accession wasn't found in the current list of genomes")
			overview[genome.split("=")[0]] = "failed"
			continue
		overview[genome.split("=")[0].replace("_", " ")] = download(species_data, genome.split("=")[0])
		continue

	genome = genome.replace("_", " ")
	taxid = get_taxid(genome)
	if taxid:
		print(now(), "taxid for",genome,"is",taxid)
	else:
		print(now(), "Unable to find taxid for", genome, ". Was the name misspelled?")
		overview[genome] = "failed"
		continue
	print(now(),"Checking number of available genomes...")
	#species_data = data[data["organism_name"] == genome]
	species_data = data[data.organism_name.str.contains(genome)]
	dates = list(species_data["seq_rel_date"])
	dates.sort()
	print(now(), "Found", species_data.shape[0], "assemblies for", genome)
	if (species_data.shape[0] == 0):
		continue
	if (species_data.shape[0] == 1):
		if str(species_data.iloc[0]["taxid"]) == taxid:
			print(now(), "taxids of species and available genome entry match. Will proceed to genome download...")
			overview[genome] = download(species_data, genome)
			continue
		else:
			print(now(), "The genome availabe for", genome, "with taxid:",species_data.iloc[0]["taxid"],"seems to have a different taxid as the species specified:", taxid, ". Nothing will be downloaded, please check manually and adjust the name.")
			overview[genome] = "failed"
			continue
	else:
		print(now(), "Will check if all assemblies have the same taxid...")
		if len(set(list(species_data["taxid"]))) == 1:
			print(now(), "All genomes have the same the same taxid.", list(species_data["taxid"])[0])
			if not str(list(species_data["taxid"])[0]) == taxid:
				print(now(), "The genome availabe for", genome, "with taxid:",species_data.iloc[0]["taxid"],"seems to have a different taxid as the species specified:", taxid, ". Nothing will be downloaded, please check manually and adjust the name. This can happen if the genome is deposited under an infraspecific name such as strains, subspecies etc.")
				overview[genome] = "failed"
				#print(now(), "taxid:", list(species_data["taxid"])[0], "and species name", genome, "do not match. Will skip this species")
				continue	
			else:
				if "reference genome" in list(species_data["refseq_category"]):
					print(now(), "Found reference genome for", genome)
					overview[genome] = download(species_data[species_data.refseq_category == "reference genome"], genome)
					continue
				elif "representative genome" in list(species_data["refseq_category"]):
					print(now(), "Found representative genome for", genome)
					overview[genome] = download(species_data[species_data.refseq_category == "representative genome"], genome)
					continue
				else:
					print(now(), "There is no reference or representative genome available for %s, will download latest genome..." % genome)
					dates = list(species_data["seq_rel_date"])
					latest_date = dates[-1]
					print(now(), "Latest genome of",genome,"was released:", latest_date, ". Will try to download this genome...")
					overview[genome] = download(species_data[species_data.seq_rel_date == "latest_date"], genome)	
					continue
		else:
			print(now(), "Potential genomes have multiple different taxids. Will try to resolve this...")
			taxids = [str(item) for item in set(list(species_data["taxid"]))]
			if taxid not in taxids:
				print(now(), "tax id of species is not in the list of taxon ids with potential genomes. Will not download anything. Please resolve manually.")
				overview[genome] = "failed"
				continue
			taxid_species_data = species_data[species_data["taxid"] == int(taxid)]
			if (taxid_species_data.shape[0] == 1):
				print(now(), "A single genome with the correct taxid remaining.Will try to download this genome now...")
				overview[genome] = download(taxid_species_data[taxid_species_data.taxid == int(taxid)], genome)
				continue
			else:
				print(now(), taxid_species_data.shape[0], "genomes with correct taxid remaining. Checking for reference genome ...")
				if "reference genome" in list(taxid_species_data["refseq_category"]):
					print(now(), "Found reference genome for", genome)
					overview[genome] = download(taxid_species_data[taxid_species_data.refseq_category == "reference genome"], genome)
					continue
				elif "representative genome" in list(taxid_species_data["refseq_category"]):
					print(now(), "Found representative genome for", genome)
					overview[genome] = download(taxid_species_data[taxid_species_data.refseq_category == "representative genome"], genome)
					continue
				else:
					print(now(), "There is now reference or representative genome available for, will download latest genome...", genome)
					dates = list(taxid_species_data["seq_rel_date"])
					latest_date = dates[-1]
					print(now(), "Latest genome of",genome,"was released:", latest_date, ". Will try to download this genome...")
					overview[genome] = download(taxid_species_data[taxid_species_data.seq_rel_date == latest_date], genome)
					continue
	# All possible cases should be covered by the code above, in case something goes wrong nonetheless, treat genome as failed:
	print(now(), "An unspecified problem occured for species", genome, ". Nothing was downloaded")
	overview[genome] = "failed"

print(overview)
print(now(), "Writing overview statistics files")
statsfile = open(args.outdir + "download_overview.txt", "w")
successfile = open(args.outdir + "successfully_downloaded.txt", "w")
failedfile = open(args.outdir + "not_downloaded.txt", "w")
for sp in overview.keys():
	print(sp.replace(" ","_"), ",", overview[sp], sep="", file=statsfile)
	if overview[sp] == "success":
		print(sp.replace(" ","_"), file=successfile)
	if overview[sp] == "failed":
		print(sp.replace(" ","_"), file=failedfile)
statsfile.close()
successfile.close()
failedfile.close()

