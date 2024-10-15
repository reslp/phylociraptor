import os
import sys
import glob
sys.path.insert(0, os.getcwd()+"/bin") #needed so that phylociraptor module can be imported.
import yaml
from libphylociraptor.hashing import *

# get list of containers to use:
with open("data/containers.yaml", "r") as yaml_stream:
    containers = yaml.safe_load(yaml_stream)

#configfi=str(sys.argv[sys.argv.index("--configfile")+1])
configfi=os.environ["CONFIG"]
print(now(), "CONFIGFILE:", configfi)

def get_modeltest_checkpoint(wildcards):
	if os.path.exists("results/modeltest/parameters.modeltest."+wildcards.aligner+"-"+wildcards.alitrim+"."+modeltest_hashes["iqtree"][wildcards.alitrim][wildcards.aligner]+".yaml"):
		return "results/checkpoints/modeltest/aggregate_best_models_"+wildcards.aligner+"_"+wildcards.alitrim+"."+modeltest_hashes["iqtree"][wildcards.alitrim][wildcards.aligner]+".done"
	else:
		return []

def get_input_genes(wildcards):
	bs_cutoff = int(wildcards.bootstrap)
	list_of_genes = []
	if os.path.exists("results/modeltest/genetree_filter_" + wildcards.aligner + "_" + wildcards.alitrim + "."+modeltest_hashes["iqtree"][wildcards.alitrim][wildcards.aligner]+".txt"):
		# path when modeltest was run before
		with open("results/modeltest/genetree_filter_" + wildcards.aligner + "_" + wildcards.alitrim + "."+modeltest_hashes["iqtree"][wildcards.alitrim][wildcards.aligner]+".txt") as file:
			for line in file:
				gene = line.split("\t")[0]
				bs_value = int(line.strip().split("\t")[-1])
				if bs_value >= bs_cutoff:
					list_of_genes.append(gene)
	else:
		# path when modeltest was not run
		if os.path.exists("results/alignments/filtered/" + wildcards.aligner + "-" + wildcards.alitrim + "." + trimmer_hashes[wildcards.alitrim][wildcards.aligner] + "/"):
			alignments = glob.glob("results/alignments/filtered/" + wildcards.aligner + "-" + wildcards.alitrim + "." + trimmer_hashes[wildcards.alitrim][wildcards.aligner] + "/*.fas")
			list_of_genes = [alignment.split("_")[0].split("/")[-1] for alignment in alignments]
	return list_of_genes		

def get_aligners():
	aligners = config["alignment"]["method"]
	if isinstance(aligners, str):
		if ", " in aligners:
			return aligners.split(", ")
		elif " ," in aligners:
			return aligners.split(" ,")
		elif "," in aligners:
			return aligners.split(",")
		elif " " in aligners:
			return aligners.split(" ")
		else:
			return aligners
	elif isinstance(aligners, list):
		return aligners
	else:
		print("phylociraptor filter-align: The alignment methods need to be either specified as string or as python list. When this is provided as string use a space( ) or comma (,) as separator. ")
		sys.exit(1)

def get_trimmers():
	trimmers = config["trimming"]["method"]
	if isinstance(trimmers, str):
		if ", " in trimmers:
			return trimmers.split(", ")
		elif " ," in trimmers:
			return trimmers.split(" ,")
		elif "," in trimmers:
			return trimmers.split(",")
		elif " " in trimmers:
			return trimmers.split(" ")
		else:
			return trimmers
	elif isinstance(trimmers, list):
		return trimmers
	else:
		print("phylociraptor filter-align: The trimming methods need to be either specified as string or as python list. When this is provided as string use a space( ) or comma (,) as separator. ")
		sys.exit(1)

def get_treemethods():
	trees = config["mltree"]["method"]
	if isinstance(trees, str):
		if ", " in trees:
			return trees.split(", ")
		elif " ," in trees:
			return trees.split(" ,")
		elif "," in trees:
			return trees.split(",")
		elif " " in trees:
			return trees.split(" ")
		else:
			return trees
	elif isinstance(trees, list):
		return trees
	else:
		print("phylociraptor mltree: The tree methods need to be either specified as string or as python list. When this is provided as string use a space( ) or comma (,) as separator. ")
		sys.exit(1)

def get_bootstrap_cutoffs():
	bscut = config["genetree_filtering"]["bootstrap_cutoff"]
	if isinstance(bscut, str):
		if ", " in bscut:
			return bscut.split(", ")
		elif " ," in bscut:
			return bscut.split(" ,")
		elif "," in bscut:
			return bscut.split(",")
		elif " " in bscut:
			return bscut.split(" ")
		else:
			return bscut
	elif isinstance(bscut, list):
		return bscut
	else:
		print("phylociraptor: The bootstrap cutoffs need to be either specified as string or as python list. When this is provided as string use a space( ) or comma (,) as separator. ")
		sys.exit(1)

def get_genetree_methods():
	genetrees = config["filtering"]["bootstrap_cutoff"]
	if isinstance(genetrees, str):
		if ", " in genetrees:
			return genetrees.split(", ")
		elif " ," in genetrees:
			return genetrees.split(" ,")
		elif "," in genetrees:
			return genetrees.split(",")
		elif " " in genetrees:
			return genetrees.split(" ")
		else:
			return genetrees
	elif isinstance(genetrees, list):
		return genetrees
	else:
		print("phylociraptor modeltest: The genetree method needs to be either specified as string or as python list. When this is provided as string use a space( ) or comma (,) as separator. ")
		sys.exit(1)

def get_bichains():
	return [i for i in range(1, int(config["bitree"]["chains"]["phylobayes"]))]
