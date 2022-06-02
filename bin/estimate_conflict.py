#!/usr/bin/env python3

#import dendropy as dp
import sys
#from ete3 import Tree
import time
from itertools import combinations, permutations, product
import random
import multiprocessing
import pandas as pd
import datetime
import argparse
import glob
import os
import math
from rpy2.robjects.packages import importr
from rpy2 import robjects as ro

def reduce_trees(treelist, quat):
	sorted_trees = []
	i = 1
	for tree in treelist:
		tiplabels = list(tree.rx2("tip.label"))
		tf_list = []
		for tip in tiplabels:
			if tip in quat:
				continue
			else:
				tf_list.append(tip)
		tipss =ro.StrVector(tf_list)
		redtree = ape.drop_tip(tree, tipss)
		treestr2 = ape.write_tree(ape.drop_tip(tree, tipss))
		treestr = treestr2[0].strip(";")
		if len(list(redtree.rx2("tip.label"))) != 4:
			missing = [tip for tip in quat if tip not in list(redtree.rx2("tip.label"))] 
			print(quat, "is not properly represented in tree", i, "Missing tips:", missing)
		tr = [el.split(")") for el in treestr.split("(") if el]
		tt = [tip.lstrip(",").strip(",") for sublist in tr for tip in sublist if tip]
		tt = [el.strip(",").strip(":") for el in tt]
		sorted_tree = []
		tips = ""
		second_element = 0
		for el in tt:
			if "," in el:
				sorted_tree.append(",".join(sorted(el.split(","))))
			else:
				tips += el
				if second_element == 1:
					sorted_tree.append(",".join(sorted(tips.split(","))))
				if second_element == 0:
					tips += ","
					second_element = 1
		sorted_trees.append("-".join(sorted(sorted_tree)))
		i += 1
	return sorted_trees

def reformat_quat(quat):
	return quat[0] + "," + quat[1] + "-" + quat[2] + "," + quat[3]

def compare_trees(treelist, combinations):
	comparison = {}
	for comb in combinations:
		if treelist[comb[0]-1] == treelist[comb[1]-1]:
			comparison["T" + str(comb[0]) + "-T" + str(comb[1])] = 1
			#print("MATCH:", comb, treelist[comb[0]], treelist[comb[1]])	
		else:
			comparison["T" + str(comb[0]) + "-T" + str(comb[1])] = 0
			#print("DIFFERENCE:", comb, treelist[comb[0]-1], treelist[comb[1]-1])	
	return comparison


def quats_in_trees(tree, quat):
	tstring = Tree(tree.extract_tree_with_taxa_labels(quat).as_string(schema="newick"))
	tr = treestr.split("(")
	tt = [tip.lstrip(",").strip(",") for tip in sublist if tip]
	sorted_tree = []
	tips = ""
	second_element = 0
	for el in tt:
		if "," in el:
			sorted_tree.append(",".join(sorted(el.split(","))))
		else:
			tips += el
			if second_element == 1:
				sorted_tree.append(",".join(sorted(tips.split(","))))
			if second_element == 0:
				tips += ","
				second_element = 1
			
						
	sorted_trees.append("-".join(sorted(sorted_tree)))
	quat_formated = quat[0] + "," +quat[1]+"-"+quat[2]+","+quat[3]
	print(quat_formated)
	print(sorted_trees)


def random_combination(iterable, r):
	i = 0
	alltips = tuple(iterable)
	ntips = len(alltips)
	tipsrange = range(ntips)
	while i < r:
		i += 1
		yield [alltips[j] for j in random.sample(tipsrange, r)]

nn = 1
def compute_single_quat(tl,i, quat, combinations):
	#print("Computing quat:", i, end="\r")
	tl2 = reduce_trees(tl, quat)
	if tl2 == 0:
		return 0
	result = compare_trees(tl2, combinations)
	return (tl2[0], result)
	

if __name__ == '__main__':
	parser = argparse.ArgumentParser(prog="compare_trees.py", description = """This script compares trees based on quartets""", epilog = """written by Philipp Resl""")
	parser.add_argument("-i","--intrees", action="store", default="all")
	parser.add_argument("-o", "--outfile", action="store", default="output") 
	parser.add_argument("-s", "--seed", action="store")
	parser.add_argument("-t", "--threads", action="store")
	parser.add_argument("-q", "--nquartets", action="store")
	parser.add_argument("-l", "--lineagefile", action="store")
	parser.add_argument("-b", "--stopby", action="store")
	parser.add_argument("--selecttaxa", action="store")
	args = parser.parse_args() 
	start = time.time()
	# read in all parameters
	if args.seed != "random":
		seed = int(args.seed)
#		seed = 111123 # here for testing
		random.seed(seed)
		print("Random seed:", seed)
	outfile = args.outfile
	if not args.threads:
		ncpus = multiprocessing.cpu_count()
	else:
		ncpus = int(args.threads)
	if args.intrees == "all":
		iqtree_trees = glob.glob("results/phylogeny-*/iqtree/*/concat.contree")
		raxml_trees =  glob.glob("results/phylogeny-*/raxml/*/raxmlng.raxml.support")
		astral_trees = glob.glob("results/phylogeny-*/astral/*/species_tree.tre")
		#all_trees = ",".join(iqtree_trees + raxml_trees + astral_trees)
		treenames = iqtree_trees +raxml_trees + astral_trees
	
	else:
		#treenames = ["species_tree.tre", "concat.contree", "species_tree.tre2", "concat.contree2", "concat.contree3", "concat.contree3"]
		treenames = args.intrees.split(",")
	i = 1
	tl = []
	ape = importr('ape')
	for name in treenames:
		print("Tree", i, "->", name)
		tree = ape.read_tree(name)
		tree.rx[2] = [ro.r("NULL")] * len(tree.rx2("node.label"))
		tree.rx[4] = [ro.r("NULL")] * len(tree.rx2("node.label"))
		tl.append(tree)
		#tl.append(dp.Tree.get(path=name, schema="newick"))
		i += 1
	# create selection of tips to be analyzed:
	all_tips = [tip for tree in tl for tip in list(tree.rx2("tip.label"))]
	print("Total number of taxa in all trees:", len(set(all_tips)))
	if args.selecttaxa:
		if not args.lineagefile:
			print("Taxon selection requires a lineage file")
			sys.exit(0)
		else:
			if os.path.exists(args.lineagefile):
				lineage_df = pd.read_csv(args.lineagefile, sep="\t")
				rank, taxa = args.selecttaxa.split("=")
				selected_tips = lineage_df.loc[lineage_df[rank].isin(taxa.split(",")), 'name'].to_list()
				taxon_list = [tip for tip in selected_tips if all_tips.count(tip)]
			else:
				print("Lineage file does not exist")
				sys.exit(0)
	else:
		taxon_list = [taxon for taxon in set(all_tips)]
	
	random.shuffle(taxon_list)

	print("No. of trees:", len(tl))
	combinations = combinations(range(1, len(tl)+1), 2)
	combinations = list(combinations)
	results = 0.0
	i = 0
	print("No. of Taxa for calculating quartets:", len(taxon_list))
	print("Using", ncpus, "CPUs.")	
	quartet_list = []
	print("Calculating quartets, this may take some time...")
	# decide how many quartets should be calculated
	if args.nquartets: # 1. by specifying a max number
		input_list=[]
		nquartets = int(args.nquartets)
		total_quartets = math.factorial(len(taxon_list)) / (math.factorial(len(taxon_list) - 4) * math.factorial(4))
		if total_quartets < nquartets:
			print("WARNING: The total number of possible quartets (",total_quartets,") is smaller than the specified nquartets.")
			#nquartets = int(total_quartets)
			
		print("No. of quartets to be calculated:", nquartets)
		for i in range(nquartets):
			quartet = next(random_combination(taxon_list, 4))
			quartet_list += quartet
			input_list.append([tl, i, quartet, combinations])
		pool = multiprocessing.Pool(ncpus)
		results = pool.starmap_async(compute_single_quat, input_list)
		final_results = results.get()
		pool.close()
	elif args.stopby: # 2. by stopping criterion (can still be extended)
		which, num = args.stopby.split("=")
		num = int(num)
		if which == "conflicts":
			print("Stopping criterion is number of conflicts:", num)
			nconflicts = 0
			ncalculations = 0
			step = 100
			final_results = []
			while nconflicts < num:
				input_list = []
				for i in range(step):
					quartet = next(random_combination(taxon_list, 4))
					quartet_list += quartet
					input_list.append([tl, i, quartet, combinations])
				pool = multiprocessing.Pool(ncpus)
				results = pool.starmap_async(compute_single_quat, input_list)
				intermed_results = results.get()
				pool.close()
				for result in intermed_results:
					if any(v == 0 for v in list(result[1].values())):
						nconflicts += 1
				final_results += intermed_results
				ncalculations += step	
				print("Calculated quartets:", ncalculations, "- Identified conflicts:",nconflicts,"of",num)
			nquartets = ncalculations
		elif which == "tipcoverage":
			print("Stopping criterion is tip coverage:", num)
			samplingdepth = 0
			ncalculations = 0
			step = 100
			final_results = []
			while True:
				input_list=[]
				for i in range(step):
					quartet = next(random_combination(taxon_list, 4))
					quartet_list += quartet
					input_list.append([tl, i, quartet, combinations])
				pool = multiprocessing.Pool(ncpus)
				results = pool.starmap_async(compute_single_quat, input_list)
				intermed_results = results.get()
				pool.close()
				if not all(quartet_list.count(v) >= num for v in taxon_list):
					depthlist = [quartet_list.count(tip) for tip in taxon_list]
					samplingdepth = sum(depthlist) / len(depthlist)
					final_results += intermed_results
					ncalculations += step	
					print("Calculated quartets:", ncalculations, "- Times each tip was sampled (mean):",round(samplingdepth, 2))
				else:
					break
			nquartets = ncalculations
			
		else:
			print("Stopping criterion not recognized.")
			sys.exit(1)
	else:
		print("Either use a stopping criterion or a maximum number of quartets")
		sys.exit(0)			
	print("Calculating quartets finished.\n")	
	print("Number of results:", len(final_results))
	results_dict = {}
	i = 1
	for result in final_results:
		if result != 0:
			results_dict[result[0]] = result[1]
			i += 1
	df = pd.DataFrame.from_dict(results_dict)	
	print("Calculating sampling coverage of tips:")
	with open(outfile +".tip_sampling_coverage.tsv", "w") as f:
		for tip in taxon_list:
			print(tip, quartet_list.count(tip), file=f)

	print("Saving results. Prefix:", outfile)
	df.transpose().to_csv(outfile + ".quartets.csv", index_label="quartet", quoting=False) 
	
	print("Calculating pairwise comparison of trees (% similarity based on quartetts)...")
	dif = len(df.columns.to_list()) # get total number of actually computed quartets. This can be smaller as nquartets eg. when there are few tips!
	df_comp = df.sum(axis=1).sort_index()/dif
	df_comp = df_comp.sort_index()	
	combinations = df_comp.index.to_list()
	index_names = ["T"+ str(i+1) for i in range(0,len(tl))]
	values = [1 for item in index_names]
	df_dict = {}
	for combination in index_names:
		df_dict[combination] = values
	df = pd.DataFrame(df_dict)
	df.index = index_names 
	for tr1 in index_names:
		for tr2 in index_names:
			if "-".join([tr1, tr2]) in df_comp.index.to_list():
				df.loc[tr1, tr2] = df_comp.transpose()[tr1+"-"+tr2]
				df.loc[tr2, tr1] = df_comp.transpose()[tr1+"-"+tr2]
	df.to_csv(outfile + ".similarity_matrix.csv")
	print("Writing list of trees...")
	with open(outfile + ".treelist.tsv", "w") as f:
		i = 1
		for tree in treenames:
			print("T"+str(i), tree, file=f, sep="\t")
			i += 1
	end = time.time()
	print("Tree comparison is done.")
	print("\nThe following output has been created:")
	print(outfile+".quartets.csv contains the raw quartet occurences.")
	print(outfile+".similarity_matrix.csv contains the % similarity between pairs of trees.")
	print(outfile+".treelist.tsv contains all trees used and the corresponding tree number.")
	print(outfile+".tip_sampling_coverage.tsv contains counts how often each tip was sampled.")
	print()
	print("Time to execute:", datetime.timedelta(seconds=round(end-start)))


