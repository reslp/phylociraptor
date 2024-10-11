import os
import glob

outfile_check_dict = {
"setup": ["results/checkpoints/modes/phylogenomics_setup.done"],
"orthology": ["results/checkpoints/modes/orthology.HASH.done", "results/checkpoints/modes/phylogenomics_setup.done"],
"filter-orthology": ["results/checkpoints/modes/filter_orthology.HASH.done", "results/checkpoints/modes/orthology.PREVIOUS.done"],
"align": [ "results/checkpoints/modes/align.HASH.done", "results/checkpoints/modes/filter_orthology.PREVIOUS.done"],
"filter-align": ["results/checkpoints/modes/align.PREVIOUS.done", "results/checkpoints/modes/filter_align.HASH.done"],
"speciestree": ["results/checkpoints/modes/modeltest.PREVIOUS.done", "results/checkpoints/modes/speciestree.HASH.done" ],
"njtree": ["results/checkpoints/modes/njtree.HASH.done", "results/checkpoints/modes/modeltest.PREVIOUS.done"],
"mltree": ["results/checkpoints/modes/modeltest.PREVIOUS.done", "results/checkpoints/modes/mltree.HASH.done"],
"bitree": ["results/checkpoints/modes/modeltest.PREVIOUS.done", "results/checkpoints/modes/bitree.HASH.done"],
"modeltest": ["results/checkpoints/modes/filter_align.PREVIOUS.done", "results/checkpoints/modes/modeltest.HASH.done"]
}

previous_steps = {"setup": "setup", "orthology": "setup", "filter-orthology": "orthology", "align": "filter-orthology", "filter-align": "align", "modeltest": "filter-align", "njtree": "modeltest", "mltree": "modeltest", "speciestree": "modeltest", "bitree": "modeltest"}
steps_to_check = ["setup", "orthology", "filter-orthology", "align", "filter-align", "modeltest", "njtree", "speciestree", "mltree", "bitree"]
step_status = {"setup": 0, "orthology": 0, "filter-orthology": 0, "align": 0, "filter-align": 0, "modeltest": 0, "njtree": 0, "mltree": 0, "speciestree": 0, "bitree": 0}

def has_outfile(mode="", previous_mode = "", hashes={}, previous_hashes={}, debug=False, verbose=False):
	if debug:
		print("Mode:", mode, "hash:", hashes[mode]["global"], "   Previous Mode:", previous_mode, "hash:", previous_hashes[previous_mode]["global"])		
	for f in outfile_check_dict[mode]:
		f = f.replace("HASH", hashes[mode]["global"]).replace("PREVIOUS",previous_hashes[previous_mode]["global"])
		if not os.path.isfile(f):
			print(" This final output file for this step is missing:", f)
			return False
	return True

def return_results_location(mode = "", hashes={}, debug=False, verbose=False):
	if mode == "setup":
		return "  ./results/assemblies"
	if mode == "orthology":
		return "  ./results/orthology/orthology_table." + hashes[mode]["global"] + ".txt"
	if mode == "filter-orthology":
		out = "  ./results/orthology/single-copy-orthologs." + hashes[mode]["global"] + "\n"
		return out.rstrip()
	if mode == "align":
		out = ""
		for aligner in hashes[mode]["per"]:
			out += "  ./results/alignments/full/" + aligner + "." + hashes[mode]["per"][aligner] + "\n"
		return out.rstrip()
	if mode == "filter-align":
		out = ""
		for trimmer in hashes[mode]["per"]:
			for aligner in hashes[mode]["per"][trimmer]:
				out += "  ./results/alignments/filtered/" + aligner + "-" + trimmer + "." + hashes[mode]["per"][trimmer][aligner] + "\n"
		return out.rstrip()
	if mode == "modeltest":
		out = ""
		for genetree in hashes[mode]["per"]:
			for trimmer in hashes[mode]["per"][genetree]:
				for aligner in hashes[mode]["per"][genetree][trimmer]:
					out += "  ./results/modeltest/" + aligner + "-" + trimmer + "." + hashes[mode]["per"]["iqtree"][trimmer][aligner] + "\n"
		return out.rstrip()
	if mode == "njtree":
		out = ""
		if verbose:
			for bscutoff in hashes[mode]["per"]:
				for njtree in hashes[mode]["per"][bscutoff]:
					for genetree in hashes[mode]["per"][bscutoff][njtree]:
						for trimmer in hashes[mode]["per"][bscutoff][njtree][genetree]:
							for aligner in hashes[mode]["per"][bscutoff][njtree][genetree][trimmer]:
								out += "  ./results/phylogeny/" + njtree + "/bootstrap-cutoff-" + bscutoff + "/" +  aligner + "-" + trimmer + "." + hashes[mode]["per"][bscutoff][njtree][genetree][trimmer][aligner] + "\n"
		else:
			for bscutoff in hashes[mode]["per"]:
				for njtree in hashes[mode]["per"][bscutoff]:
					out += "  ./results/phylogeny/" + njtree + "/bootstrap-cutoff-" + bscutoff + "\n"
		return out.rstrip()
	if mode == "mltree":
		out = ""
		if verbose:
			for bscutoff in hashes[mode]["per"]:
				for mltree in hashes[mode]["per"][bscutoff]:
					for genetree in hashes[mode]["per"][bscutoff][mltree]:
						for trimmer in hashes[mode]["per"][bscutoff][mltree][genetree]:
							for aligner in hashes[mode]["per"][bscutoff][mltree][genetree][trimmer]:
								out += "  ./results/phylogeny/" + mltree + "/bootstrap-cutoff-" + bscutoff + "/" +  aligner + "-" + trimmer + "." + hashes[mode]["per"][bscutoff][mltree][genetree][trimmer][aligner] + "\n"
		else:
			for bscutoff in hashes[mode]["per"]:
				for mltree in hashes[mode]["per"][bscutoff]:
					out += "  ./results/phylogeny/" + mltree + "/bootstrap-cutoff-" + bscutoff + "\n"
		return out.rstrip()
	if mode == "bitree":
		out = ""
		if verbose:
			for bscutoff in hashes[mode]["per"]:
				for mltree in hashes[mode]["per"][bscutoff]:
					for genetree in hashes[mode]["per"][bscutoff][bitree]:
						for trimmer in hashes[mode]["per"][bscutoff][bitree][genetree]:
							for aligner in hashes[mode]["per"][bscutoff][bitree][genetree][trimmer]:
								out += "  ./results/phylogeny/" + bitree + "/bootstrap-cutoff-" + bscutoff + "/" +  aligner + "-" + trimmer + "." + hashes[mode]["per"][bscutoff][bitree][genetree][trimmer][aligner] + "\n"
		else:
			for bscutoff in hashes[mode]["per"]:
				for mltree in hashes[mode]["per"][bscutoff]:
					out += "  ./results/phylogeny/" + bitree + "/bootstrap-cutoff-" + bscutoff + "\n"
		return out.rstrip()
	if mode == "speciestree":
		out = ""
		if verbose:
			for bscutoff in hashes[mode]["per"]:
				for sptree in hashes[mode]["per"][bscutoff]:
					for genetree in hashes[mode]["per"][bscutoff][sptree]:
						for trimmer in hashes[mode]["per"][bscutoff][sptree][genetree]:
							for aligner in hashes[mode]["per"][bscutoff][sptree][genetree][trimmer]:
								out += "  ./results/phylogeny/" + sptree + "/bootstrap-cutoff-" + bscutoff + "/" +  aligner + "-" + trimmer + "." + hashes[mode]["per"][bscutoff][sptree][genetree][trimmer][aligner] + "\n"
		else:
			for bscutoff in hashes[mode]["per"]:
				for sptree in hashes[mode]["per"][bscutoff]:
					out += "  ./results/phylogeny/" + sptree + "/bootstrap-cutoff-" + bscutoff + "\n"
		return out.rstrip()

def check_is_running(mode="", previous_mode="", hashes={}, previous_hashes={}, debug=False, verbose=False):
	#print(hashes[mode])
	#print("   ")
	#print(previous_hashes)
	if mode == "align":
		genes = []
		# first get all sequences files after filter-orthology
		if not os.path.isdir("results/alignments"):
			return
		for gene in os.listdir("results/orthology/single-copy-orthologs."+previous_hashes[previous_mode]["global"]+"/"):
			gene = os.path.splitext(gene)[0]
			genes.append(gene.split("_")[0])
		# now check if the same files are in align
		for directory in [h + "." + hashes[mode]["per"][h] for h in hashes[mode]["per"].keys()]: #we need to avoid the yaml files.
			if os.path.isdir("results/alignments/full/"+directory):
				finished_alignments = [os.path.basename(al).split("_")[0] for al in glob.glob("results/alignments/full/"+directory+"/*.fas")]
				missing = [al for al in genes if al not in finished_alignments]
				if missing:
					print("\tMissing", directory, "alignments for", len(missing), "of", len(genes), "genes.")
					if verbose:
						print("\t These alignments are missing:")
						print("\t", ",".join(missing))
				else:
					print("\t"+str(len(genes)),"alignments for", directory, "done.")
			else:
					print("\tMissing", directory, "alignments for", len(genes), "of", len(genes), "genes.")
	if mode == "filter-align":
		# work on trimmed alignments
		if not os.path.isdir("results/alignments/trimmed"):
			print("Alignment trimming has not yet started")
			return
		for trimmer in hashes[mode]["per"].keys():
			for aligner in hashes[mode]["per"][trimmer].keys():
				alitrim = aligner + "-" + trimmer + "." + hashes[mode]["per-trimming"][trimmer][aligner]		
				genes = []
				if os.path.isdir("results/alignments/full/" + aligner + "."+previous_hashes[previous_mode]["per"][aligner]):
					genes = [al.split("/")[-1].split(".")[0].split("_")[0] for al in glob.glob("results/alignments/full/" + aligner + "." + previous_hashes[previous_mode]["per"][aligner] + "/*.fas")]
				if os.path.isdir("results/alignments/trimmed/"+alitrim):
					finished_alignments = [os.path.basename(al).split("_")[0] for al in glob.glob("results/alignments/trimmed/" + alitrim + "/*.fas")]
					missing = [al for al in genes if al not in finished_alignments]
					if missing:
						print("\tMissing", alitrim, "trimmed alignments for", len(missing), "of", len(genes), "genes.")
						if verbose:
							print("\t These trimmed alignments are missing:")
							print("\t", ",".join(missing))
					else:
						print("\t"+str(len(genes)),"trimmed alignments for", alitrim, "done.")
				else:
					print("\tMissing", alitrim, "trimmed alignments for", len(genes), "of", len(genes), "genes.")
		print()
		# work on filtered alignments
		if not os.path.isdir("results/alignments/filtered"):
			print("\tAlignment filtering has not yet started")
			return
		for trimmer in hashes[mode]["per"].keys():
			for aligner in hashes[mode]["per"][trimmer].keys():
				alitrim = aligner + "-" + trimmer + "." + hashes[mode]["per"][trimmer][aligner]		
				#print("results/alignments/full/" + aligner + "."+previous_hashes[previous_mode]["per"][aligner])			
				genes = []
				if os.path.isdir("results/alignments/full/" + aligner + "."+previous_hashes[previous_mode]["per"][aligner]):
					genes = [al.split("/")[-1].split(".")[0].split("_")[0] for al in glob.glob("results/alignments/full/" + aligner + "." + previous_hashes[previous_mode]["per"][aligner] + "/*.fas")]
				if os.path.isdir("results/alignments/filtered/"+alitrim):
					finished_alignments = [os.path.basename(al).split("_")[0] for al in glob.glob("results/alignments/filtered/" + alitrim + "/*.fas")]
					missing = [al for al in genes if al not in finished_alignments]
					if missing:
						print("\tFound", len(genes) - len(missing) ,"filtered alignments for", alitrim)
						if verbose:
							print("\t These filtered alignments are missing or did not survive filtering:")
							print("\t", ",".join(missing))
					else:
						print("\t"+str(len(genes)),"filtered alignments for", alitrim, "done.")
				else:
					print("\tFound", len(genes) - len(missing), "filtered alignments for", alitrim)
	if mode == "modeltest":
		if not os.path.isdir("results/checkpoints/modeltest"):
			print("Modeltesting has not yet started")
			return
		for trimmer in hashes[mode]["per"]["iqtree"].keys():
			for aligner in hashes[mode]["per"]["iqtree"][trimmer].keys():
				alitrim = aligner + "-" + trimmer + "." + previous_hashes[previous_mode]["per"][trimmer][aligner]
				mt = aligner + "-" + trimmer + "." + hashes[mode]["per"]["iqtree"][trimmer][aligner]
				genes_to_check = [gene.split("/")[-1].split("_")[0] for gene in glob.glob("results/alignments/filtered/"+alitrim+"/*.fas")]
				if os.path.isdir("results/modeltest/" + mt):
					finished_genetrees = [tree.split("/")[-1].split(".")[0] for tree in glob.glob("results/modeltest/"+mt+"/*/*.treefile")]
					missing = [tr for tr in genes_to_check if tr not in finished_genetrees]
					if missing:
						print("\tMissing", mt, "genetrees for", len(missing), "of", len(genes_to_check), "genes.")
						if verbose:
							print("\t These genetrees are missing:")
							print("\t", ",".join(missing))
					else:
						print("\t"+str(len(genes_to_check)),"genetrees for", mt, "done.")
				else:
					print("\tMissing", mt, "genetrees for", len(genes_to_check), "of", len(genes_to_check), "genes.")
	if mode == "njtree":
		#print(hashes[mode])
		#print("  ")
		#print(previous_hashes[previous_mode])

		for bs in hashes[mode]["per"].keys():
			print(" Bootstrap-cutoff:", bs)
			if not os.path.isdir("results/phylogeny/bootstrap-cutoff-"+bs+"/quicktree"):
				print("\tAnalysis for boostrap-cutoff "+bs+" not yet started.")
				continue
			missing = []
			for trimmer in hashes[mode]["per"][bs]["quicktree"]["iqtree"].keys(): #not sure why iqtree is in this
				for aligner in hashes[mode]["per"][bs]["quicktree"]["iqtree"][trimmer].keys(): #not sure why iqtree is in this
					alitrim = aligner + "-" + trimmer + "." + previous_hashes[previous_mode]["per"]["iqtree"][trimmer][aligner]
					genes_to_check = [gene.split("/")[-1].split("_")[0] for gene in glob.glob("results/alignments/filtered/"+alitrim+"/*.fas")]
					njtree = aligner + "-" + trimmer + "." + hashes[mode]["per"][bs]["quicktree"]["iqtree"][trimmer][aligner]	
					if os.path.isfile("results/phylogeny/quicktree/bootstrap-cutoff-"+bs+"/"+njtree+"/njtree.tre"):
						if verbose:
							print("\tnjtree for", njtree, "done.")
					else:
						missing.append(njtree)
						if verbose:
							print("\tMissing njtree for", njtree)
			if len(missing) > 0:
				print("\tTotal number of missing QUICKTREE trees:", len(missing))
	if mode == "mltree":
		for bs in hashes[mode]["per"].keys():
			print(" Bootstrap-cutoff:", bs)
			missing = []
			if "iqtree" in hashes[mode]["per"][bs].keys(): # this can be made more efficient!
				if not os.path.isdir("results/phylogeny/iqtree/bootstrap-cutoff-"+bs):
					print("\tIQ-Tree analysis for boostrap-cutoff "+bs+" not yet started.")
				else:
					for trimmer in hashes[mode]["per"][bs]["iqtree"]["iqtree"].keys(): #not sure why iqtree is in this
						for aligner in hashes[mode]["per"][bs]["iqtree"]["iqtree"][trimmer].keys(): #not sure why iqtree is in this
							alitrim = aligner + "-" + trimmer + "." + previous_hashes[previous_mode]["per"]["iqtree"][trimmer][aligner]
							genes_to_check = [gene.split("/")[-1].split("_")[0] for gene in glob.glob("results/modeltest/"+alitrim+"/*.treefile")]
							mltree = aligner + "-" + trimmer + "." + hashes[mode]["per"][bs]["iqtree"]["iqtree"][trimmer][aligner]	
							if os.path.isfile("results/phylogeny/iqtree/bootstrap-cutoff-"+bs+"/"+mltree+"/concat.contree"):
								if verbose:
									print("\tIQ-Tree for", mltree, "done.")
							else:
								missing.append(mltree)
								if verbose:
									print("\tMissing IQ-Tree for", mltree)
					if len(missing) > 0:
						print("\tTotal number of missing IQ-Tree trees:", len(missing))
					else:
						print("\tAll IQ-Tree trees finished.")
			print()
			missing = []
			if "raxml" in hashes[mode]["per"][bs].keys(): # this can be made more efficient!
				if not os.path.isdir("results/phylogeny/raxml/bootstrap-cutoff-"+bs):
					print("\tRAxML analysis for boostrap-cutoff-"+bs+" not yet started.")
				else:
					for trimmer in hashes[mode]["per"][bs]["raxml"]["iqtree"].keys(): #not sure why iqtree is in this
						for aligner in hashes[mode]["per"][bs]["raxml"]["iqtree"][trimmer].keys(): #not sure why iqtree is in this
							alitrim = aligner + "-" + trimmer + "." + previous_hashes[previous_mode]["per"]["iqtree"][trimmer][aligner]
							genes_to_check = [gene.split("/")[-1].split("_")[0] for gene in glob.glob("results/modeltest/"+alitrim+"/*.treefile")]
							mltree = aligner + "-" + trimmer + "." + hashes[mode]["per"][bs]["raxml"]["iqtree"][trimmer][aligner]	
							if os.path.isfile("results/phylogeny/raxml/bootstrap-cutoff-"+bs+"/"+mltree+"/raxmlng.raxml.bestTree"):
								if verbose:
									print("\tRAxML-NG for", mltree, "done.")
							else:
								missing.append(mltree)
								if verbose:
									print("\tMissing RAxML-NG for", mltree)
					if len(missing) > 0:
						print("\tTotal number of missing RAxML-NG trees:", len(missing))
					else:
						print("\tAll RAxML-NG trees finished.")
	if mode == "bitree":
		for bs in hashes[mode]["per"].keys():
			print(" Bootstrap-cutoff:", bs)
			missing = []
			if "phylobayes" in hashes[mode]["per"][bs].keys(): # this can be made more efficient!
				if not os.path.isdir("results/phylogeny/phylobayes/bootstrap-cutoff-"+bs):
					print("\tPhylobayes analysis for boostrap-cutoff "+bs+" not yet started.")
				else:
					for trimmer in hashes[mode]["per"][bs]["phylobayes"]["iqtree"].keys(): #not sure why iqtree is in this
						for aligner in hashes[mode]["per"][bs]["phylobayes"]["iqtree"][trimmer].keys(): #not sure why iqtree is in this
							alitrim = aligner + "-" + trimmer + "." + previous_hashes[previous_mode]["per"]["iqtree"][trimmer][aligner]
							genes_to_check = [gene.split("/")[-1].split("_")[0] for gene in glob.glob("results/modeltest/"+alitrim+"/*.treefile")]
							mltree = aligner + "-" + trimmer + "." + hashes[mode]["per"][bs]["phylobayes"]["iqtree"][trimmer][aligner]	
							if os.path.isfile("results/phylogeny/phylobayes/bootstrap-cutoff-"+bs+"/"+mltree+"/concat.contree"):
								if verbose:
									print("\tPhylobayes tree for", bitree, "done.")
							else:
								missing.append(bitree)
								if verbose:
									print("\tMissing Phylobayes tree for", bitree)
					if len(missing) > 0:
						print("\tTotal number of missing Phylobayes trees:", len(missing))
					else:
						print("\tAll Phylobayes trees finished.")
			print()
	if mode == "speciestree":
		for bs in hashes[mode]["per"].keys():
			print(" Bootstrap-cutoff:", bs)
			if not os.path.isdir("results/phylogeny/bootstrap-cutoff-"+bs):
				print("\tAnalysis for boostrap-cutoff "+bs+" not yet started.")
				continue
			missing = []
			for trimmer in hashes[mode]["per"][bs]["astral"]["iqtree"].keys(): #not sure why iqtree is in this
				for aligner in hashes[mode]["per"][bs]["astral"]["iqtree"][trimmer].keys(): #not sure why iqtree is in this
					alitrim = aligner + "-" + trimmer + "." + previous_hashes[previous_mode]["per"]["iqtree"][trimmer][aligner]
					genes_to_check = [gene.split("/")[-1].split("_")[0] for gene in glob.glob("results/modeltest/"+alitrim+"/*.treefile")]
					sptree = aligner + "-" + trimmer + "." + hashes[mode]["per"][bs]["astral"]["iqtree"][trimmer][aligner]	
					if os.path.isfile("results/phylogeny/astral/bootstrap-cutoff-"+bs+"/"+sptree+"/species_tree.tre"):
						if verbose:
							print("\tSpeciestree for", sptree, "done.")
					else:
						missing.append(sptree)
						if verbose:
							print("\tMissing speciestree for", sptree)
			if len(missing) > 0:
				print("\tTotal number of missing ASTRAL trees:", len(missing))
