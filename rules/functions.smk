import sys
import hashlib
import yaml

#configfi=str(sys.argv[sys.argv.index("--configfile")+1])
configfi="data/config.yaml"

###all kinds of functions related to the hash functionality

def get_hash(previous, string_of_dict_paths, yamlfile):
	import collections
	print("### GET HASH: "+string_of_dict_paths,yamlfile+" ###")
	dict_to_hash = {}
	with open(yamlfile) as f:
		my_dict = yaml.safe_load(f)
	dic_list = string_of_dict_paths.split(" ")
	for d in dic_list:
		print(d)
		l = d.split(",")
		if len(l) == 2:
			print("level 2")
			if not l[0] in dict_to_hash:
				dict_to_hash[l[0]] = {l[1]: str(my_dict[l[0]][l[1]])}
			elif not l[1] in dict_to_hash[l[0]]:
				dict_to_hash[l[0]][l[1]] = str(my_dict[l[0]][l[1]])
		if len(l) == 3:
			print("level 3")
			if not l[0] in dict_to_hash:
				dict_to_hash[l[0]] = {l[1]: {l[2]: str(my_dict[l[0]][l[1]][l[2]])}}
			elif not l[1] in dict_to_hash[l[0]]:
				dict_to_hash[l[0]][l[1]] = {l[2]: str(my_dict[l[0]][l[1]][l[2]])}
			elif not l[2] in dict_to_hash[l[0]][l[1]]:
				dict_to_hash[l[0]][l[1]][l[2]] = str(my_dict[l[0]][l[1]][l[2]])
			
	print("FINAL: "+str(dict_to_hash))
	ordered = collections.OrderedDict(dict_to_hash)
	print(str(ordered))
	combined = str(previous+str(ordered))
	print(combined)
	hash = hashlib.shake_256(combined.encode()).hexdigest(5)
	print(hash)
	print("### DONE HASH ###")
	return hash



def collect_hashes(mode):
	hashes = {}

	#orthology
	hashes['orthology'] = get_hash("", "orthology,method orthology,busco_options,set orthology,busco_options,version orthology,busco_options,mode orthology,busco_options,augustus_species orthology,busco_options,additional_parameters", configfi)
	if mode == "orthology":
		print("Gathered hashes until 'orthology'")
		print(hashes['orthology'])
		return hashes
	
	#filter-orthology
	hashes['filter-orthology'] = get_hash(hashes['orthology'], "filtering,dupseq filtering,cutoff filtering,minsp filtering,seq_type filtering,exclude_orthology", configfi)

	if mode == "filter-orthology":
		print("Gathered hashes until 'filter-orthology'")
		print(hashes['filter-orthology'])
		return hashes

	#align
	hashes['align'] = {}
	for a in config["alignment"]["method"]:
		hashes['align'][a] = get_hash(hashes['filter-orthology'], "alignment,options,"+a, configfi)

	if mode == "align":
		print("Gathered hashes until 'align'")
		print(hashes['align'])
		return hashes

	#filter-alignment
	hashes['filter-align'] = {}
	for t in config["trimming"]["method"]:
		hashes['filter-align'][t] = {}
		for a in hashes['align'].keys():
			if os.path.isfile("results/alignments/full/"+a+"."+hashes['align'][a]+"/parameters.yaml"):
				hashes['filter-align'][t][a] = get_hash(hashes['align'][a], "trimming,options,"+t+" trimming,min_parsimony_sites", configfi)
			else:
				print("Please doublecheck if the stage 'align' was run with the parameters currently specified in "+configfi)
				sys.exit()

	if mode == "filter-align":
		print("Gathered hashes until 'filter-align'")
		print(hashes['filter-align'])
		return hashes

	#modeltest
	hashes['modeltest'] = {}
	for m in config["modeltest"]["method"]:
		hashes['modeltest'][m] = {}
		for t in config["trimming"]["method"]:
			hashes['modeltest'][m][t] = {}
			for a in hashes['align'].keys():
				if os.path.isfile("results/alignments/parameters."+a+"-"+t+"."+hashes['filter-align'][t][a]+".yaml"):
					hashes['modeltest'][m][t][a] = get_hash(hashes['filter-align'][t][a], "modeltest,options,"+m+" modeltest,bootstrap", configfi)
				else:
					print("Please doublecheck if the stage 'filter-align' was run with the parameters currently specified in "+configfi)
					sys.exit()

	if mode == "modeltest":
		print("Gathered hashes until 'modeltest'")
		print(hashes['modeltest'])
		return hashes

	hashes['tree_inference'] = {}
	#speciestree
	if mode == "speciestree":
		for i in config["speciestree"]["method"]:
			hashes['tree_inference'][i] = {}
			for m in hashes['modeltest'].keys():
				hashes['tree_inference'][i][m] = {}
				for t in hashes['filter-align'].keys():
					hashes['tree_inference'][i][m][t] = {}
					for a in hashes['align'].keys():
						if os.path.isfile("results/modeltest/parameters."+a+"-"+t+"."+hashes['modeltest'][m][t][a]+".yaml"):
							hashes['tree_inference'][i][m][t][a] = get_hash(hashes['modeltest'][m][t][a], "speciestree,options,"+i+" speciestree,include", configfi)
						else:
							print("Please doublecheck if the stage 'modeltest' was run with the parameters currently specified in "+configfi)
							sys.exit()

		print("Gathered hashes until 'speciestree'")
		print(hashes['tree_inference'])
		return hashes

	###################################
	#mltree	
	if mode == "mltree":
		for i in config["mltree"]["method"]:
			hashes['tree_inference'][i] = {}
			for m in hashes['modeltest'].keys():
				hashes['tree_inference'][i][m] = {}
				for t in hashes['filter-align'].keys():
					hashes['tree_inference'][i][m][t] = {}
					for a in hashes['align'].keys():
						if os.path.isfile("results/modeltest/parameters."+a+"-"+t+"."+hashes['modeltest'][m][t][a]+".yaml"):
							hashes['tree_inference'][i][m][t][a] = get_hash(hashes['modeltest'][m][t][a], "speciestree,options,"+i+" speciestree,include", configfi)
						else:
							print("Please doublecheck if the stage 'modeltest' was run with the parameters currently specified in "+configfi)
							sys.exit()

		print("Gathered hashes until 'mltree'")
		print(hashes['mltree'])
		return hashes


def trigger(current_yaml, config_yaml):
	#the function triggers the start of a submodule if parameter files aren't there or when there are differences to the config file
	if os.path.isfile(current_yaml):
#		print("Reading in file")
		with open(current_yaml) as f:
			my_dict = yaml.safe_load(f)
		print(str(my_dict))
		if find_top(my_dict, config):
			print("ALL GOOD")
			return "None"
		else:
			return config_yaml
	else:
		print("File "+current_yaml+" is not there")
		return config_yaml


def get_all_keys(d):
	#finds all keys in a nested dictionary
	outlist=[]
	for key, value in d.items():
		yield key
		if isinstance(value, dict):
			yield from get_all_keys(value)

def search(d, k, path=None):
	#finds the path to all keys in a nested dictionary
	if path is None:
		path = []
    
	# Reached bottom of dict - no good
	if not isinstance(d, dict):
		return False
    
	# Found it!
	if k in d.keys():
		path.append(k)
		return path
    
	else:
		check = list(d.keys())
		# Look in each key of dictionary
		while check:
			first = check[0]
			# Note which we just looked in
			path.append(first)
			if search(d[first], k, path) is not False:
				break
			else:
				# Move on
				check.pop(0)
				path.pop(-1)
		else:
			return False
		return path

def find_top(t, against):
	#finds all terminal values and the paths in a nested dictionary and compares to a second dictionary
	#returns True if no difference
	#return False if difference
	for x in get_all_keys(t):
		path = search(t, x)
		if len(path) == 1:
			if not isinstance(t[path[0]], dict):
				print("1 - key: "+x+" - "+str(search(t, x))+" - "+t[path[0]])
		if len(path) == 2:
			if not isinstance(t[path[0]][path[1]], dict):
				print("2 - key: "+x+" - "+str(search(t, x))+" - "+t[path[0]][path[1]])
				print("config: "+str(against[path[0]][path[1]]))
				from_file = str(t[path[0]][path[1]])
				from_config = str(against[path[0]][path[1]])
				if from_file == from_config:
					print("EQUAL")
				else:
					print("NOT EQUAL")
					return False
		if len(path) == 3:
			if not isinstance(t[path[0]][path[1]][path[2]], dict):
				print("3 - key: "+x+" - "+str(search(t, x))+" - "+t[path[0]][path[1]][path[2]])
				print("config: "+str(against[path[0]][path[1]][path[2]]))
				from_file = str(t[path[0]][path[1]][path[2]])
				from_config = str(against[path[0]][path[1]][path[2]])
				if from_file == from_config:
					print("EQUAL")
				else:
					print("NOT EQUAL")
					return False
	return True


####
def get_input_genes(wildcards):
	bs_cutoff = int(wildcards.bootstrap)
	list_of_genes = []
	with open("results/modeltest/genetree_filter_" + wildcards.aligner + "_" + wildcards.alitrim + "."+modeltest_hashes["iqtree"][wildcards.alitrim][wildcards.aligner]+".txt") as file:
		for line in file:
			gene = line.split("\t")[0]
			bs_value = int(line.strip().split("\t")[-1])
			if bs_value >= bs_cutoff:
				list_of_genes.append(gene)
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
	bscut = config["filtering"]["bootstrap_cutoff"]
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
