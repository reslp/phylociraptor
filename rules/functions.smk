import os
import sys
import hashlib
import yaml

#configfi=str(sys.argv[sys.argv.index("--configfile")+1])
configfi=os.environ["CONFIG"]
print("CONFIGFILE:", configfi)

def get_modeltest_checkpoint(wildcards):
	return "results/checkpoints/modeltest/aggregate_best_models_"+wildcards.aligner+"_"+wildcards.alitrim+"."+modeltest_hashes["iqtree"][wildcards.alitrim][wildcards.aligner]+".done"


###all kinds of functions related to the hash functionality

def read_file_from_yaml(read, file):
	import pandas as pd
	if not os.path.exists(file) or file == None:
		print("File not found or not specified")
		return []
	if file:
		print("reading in file:", file, read )
		if len(read) > 0:
			df = pd.read_csv(file)
			lis = read.split(",")
			for i in reversed(range(len(lis))):
				if not lis[i] in df:
					del(lis[i])
			df = df[lis]
		else:
			df = pd.read_csv(file, header=None)
		return sorted(df.values.tolist())

def get_hash(previous, string_of_dict_paths=None, yamlfile=None, returndict=False):
	import collections
	import pandas as pd
	if not string_of_dict_paths and not yamlfile:
		hash = hashlib.shake_256(str(collections.OrderedDict(previous)).encode()).hexdigest(5)
		print("### DIRECT HASH ###:", hash)
		return hash
	
	print("### GET HASH: "+string_of_dict_paths,yamlfile+" ###")
	dict_to_hash = {}
	with open(yamlfile) as f:
		my_dict = yaml.safe_load(f)
	dic_list = string_of_dict_paths.split(" ")
	for d in dic_list:
		print(d)
		read = ""
		if d.startswith("<"):
			(read,d) = d.split(">")
		l = d.split(",")
		print(l)
		lev = len(l)
		if len(l) == 1:
			print("level 1")
			string = str(my_dict[l[0]])
			dict_to_hash[l[0]] = string
			if read:
				read = read.replace("<","")
				dict_to_hash[l[0]] = {string: read_file_from_yaml(read, string)}

		if len(l) == 2:
			print("level 2")
			string = str(my_dict[l[0]][l[1]])
			if not l[0] in dict_to_hash:
				dict_to_hash[l[0]] = {l[1]: string}
			elif not l[1] in dict_to_hash[l[0]]:
				dict_to_hash[l[0]][l[1]] = string
			if read:
				read = read.replace("<","")
				dict_to_hash[l[0]][l[1]] = {string: read_file_from_yaml(read, string)}

		if len(l) == 3:
			print("level 3")
			if isinstance(my_dict[l[0]][l[1]], list):
				string = l[2]
				if not l[0] in dict_to_hash:
					dict_to_hash[l[0]] = {l[1]: string}
				elif not l[1] in dict_to_hash[l[0]]:
					dict_to_hash[l[0]][l[1]] = string
			else:
				string = str(my_dict[l[0]][l[1]][l[2]])
				if not l[0] in dict_to_hash:
					dict_to_hash[l[0]] = {l[1]: {l[2]: string}}
				elif not l[1] in dict_to_hash[l[0]]:
					dict_to_hash[l[0]][l[1]] = {l[2]: string}
				elif not l[2] in dict_to_hash[l[0]][l[1]]:
					dict_to_hash[l[0]][l[1]][l[2]] = string
			if read:
				read = read.replace("<","")
				dict_to_hash[l[0]][l[1]][l[2]] = {string: read_file_from_yaml(read, string)}
	if returndict:
		return dict_to_hash

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
	if config["orthology"]["method"] == "busco":
		hashes['orthology'] = get_hash("", "<species,web_local,mode>species orthology,method <>orthology,exclude orthology,busco_options,set orthology,busco_options,version orthology,busco_options,mode orthology,busco_options,augustus_species orthology,busco_options,additional_parameters", configfi)

	if mode == "orthology":
		print("Gathered hashes until 'orthology':", hashes['orthology'],"\n")
		return hashes
	
	#filter-orthology
	if not os.path.isfile("results/orthology/busco/params.orthology."+hashes['orthology']+".yaml"):
		print("Please doublecheck if the stage 'orthology' was run with the parameters currently specified in "+configfi)
		print("I am looking for: results/orthology/busco/params.orthology."+hashes['orthology']+".yaml")
		sys.exit()
	else:
		hashes['filter-orthology'] = get_hash(hashes['orthology'], "filtering,dupseq filtering,cutoff filtering,minsp filtering,seq_type filtering,exclude_orthology", configfi)

	if mode == "filter-orthology":
		print("Gathered hashes until 'filter-orthology'")
		print(hashes['filter-orthology'])
		return hashes

	#align
	hashes['align'] = {"global": "", "per": {}}
	if not os.path.isfile("results/orthology/busco/params.filter-orthology."+hashes['filter-orthology']+".yaml"):
		print("Please doublecheck if the stage 'filter-orthology' was run with the parameters currently specified in "+configfi)
		print("I am looking for: results/orthology/busco/params.filter-orthology."+hashes['filter-orthology']+".yaml")
		sys.exit()
	else:
		hashes['align']["global"] = get_hash(hashes['filter-orthology'], "alignment,options alignment,method", configfi)
		for a in config["alignment"]["method"]:
			hashes['align']["per"][a] = get_hash(hashes['filter-orthology'], "alignment,options,"+a, configfi)

	if mode == "align":
		print("Gathered hashes until 'align'")
		print(hashes['align'])
		return hashes

	#filter-alignment
	hashes['filter-align'] = {"global": "", "per": {}}
	if not os.path.isfile("results/alignments/full/parameters.align."+hashes['align']["global"]+".yaml"):
		print("Please doublecheck if the stage 'align' was run with the parameters currently specified in "+configfi)
		print("I am looking for: results/alignments/full/parameters.align."+hashes['align']["global"]+".yaml")
		sys.exit()
	else:
		hashes['filter-align']["global"] = get_hash(hashes['align']["global"], "trimming,method trimming,options trimming,min_parsimony_sites", configfi)
		for t in config["trimming"]["method"]:
			hashes['filter-align']["per"][t] = {}
			for a in hashes['align']["per"].keys():
				hashes['filter-align']["per"][t][a] = get_hash(hashes['align']["per"][a], "trimming,options,"+t, configfi)

	if mode == "filter-align":
		print("Gathered hashes until 'filter-align'")
		print(hashes['filter-align'])
		return hashes

	#modeltest
	hashes['modeltest'] = {"global": "", "per": {}}
	if not os.path.isfile("results/alignments/trimmed/parameters.filter-align."+hashes['filter-align']["global"]+".yaml"):
		print("Please doublecheck if the stage 'filter-align' was run with the parameters currently specified in "+configfi)
		sys.exit()
	else:
		hashes['modeltest']["global"] = get_hash(hashes['filter-align']["global"], "seed modeltest,method modeltest,options modeltest,bootstrap", configfi)
		for m in config["modeltest"]["method"]:
			hashes['modeltest']["per"][m] = {}
			for t in config["trimming"]["method"]:
				hashes['modeltest']["per"][m][t] = {}
				for a in hashes['align']["per"].keys():
					hashes['modeltest']["per"][m][t][a] = get_hash(hashes['filter-align']["per"][t][a], "seed modeltest,options,"+m+" modeltest,bootstrap", configfi)

	if mode == "modeltest":
		print("Gathered hashes until 'modeltest'")
		print(hashes['modeltest'])
		return hashes

	hashes['tree_inference'] = {"global": "", "per": {}}
	#speciestree
	if mode == "speciestree":
		if not os.path.isfile("results/modeltest/parameters.modeltest."+hashes['modeltest']["global"]+".yaml"):
			print("Please doublecheck if the stage 'modeltest' was run with the parameters currently specified in "+configfi)
			sys.exit()
		else:
			hashes['tree_inference']["global"] = get_hash(hashes['modeltest']["global"], "seed genetree_filtering,bootstrap_cutoff speciestree,method speciestree,options speciestree,include", configfi)
			for c in config["genetree_filtering"]["bootstrap_cutoff"]:
				c = str(c)
				hashes['tree_inference']["per"][c] = {}
				for i in config["speciestree"]["method"]:
					hashes['tree_inference']["per"][c][i] = {}
					for m in hashes['modeltest']["per"].keys():
						hashes['tree_inference']["per"][c][i][m] = {}
						for t in hashes['filter-align']["per"].keys():
							hashes['tree_inference']["per"][c][i][m][t] = {}
							for a in hashes['align']["per"].keys():
								hashes['tree_inference']["per"][c][i][m][t][a] = get_hash(hashes['modeltest']["per"][m][t][a], "seed speciestree,options,"+i+" speciestree,include", configfi)

		print("Gathered hashes until 'speciestree'")
		print(hashes['tree_inference'])
		return hashes

	###################################
	#mltree	
	if mode == "mltree":
		if not os.path.isfile("results/modeltest/parameters.modeltest."+hashes['modeltest']["global"]+".yaml"):
			print("Please doublecheck if the stage 'modeltest' was run with the parameters currently specified in "+configfi)
			sys.exit()
		else:
			hashes['tree_inference']["global"] = get_hash(hashes['modeltest']["global"], "seed genetree_filtering,bootstrap_cutoff mltree,method mltree,bootstrap mltree,options", configfi)
			for c in config["genetree_filtering"]["bootstrap_cutoff"]:
				c = str(c)
				hashes['tree_inference']["per"][c] = {}
				for i in config["mltree"]["method"]:
					hashes['tree_inference']["per"][c][i] = {}
					for m in hashes['modeltest']["per"].keys():
						hashes['tree_inference']["per"][c][i][m] = {}
						for t in hashes['filter-align']["per"].keys():
							hashes['tree_inference']["per"][c][i][m][t] = {}
							for a in hashes['align']["per"].keys():
								hashes['tree_inference']["per"][c][i][m][t][a] = get_hash(hashes['modeltest']["per"][m][t][a], "seed genetree_filtering,bootstrap_cutoff,"+c+" mltree,bootstrap,"+i+" mltree,options,"+i, configfi)

		print("Gathered hashes until 'mltree'")
		print(hashes['tree_inference'])
		return hashes

	#njtree	
	if mode == "njtree":
		if not os.path.isfile("results/modeltest/parameters.modeltest."+hashes['modeltest']["global"]+".yaml"):
			print("Please doublecheck if the stage 'modeltest' was run with the parameters currently specified in "+configfi)
			sys.exit()
		else:
			hashes['tree_inference']["global"] = get_hash(hashes['modeltest']["global"], "genetree_filtering,bootstrap_cutoff njtree,method njtree,options", configfi)
			for c in config["genetree_filtering"]["bootstrap_cutoff"]:
				c = str(c)
				hashes['tree_inference']["per"][c] = {}
				for i in config["njtree"]["method"]:
					hashes['tree_inference']["per"][c][i] = {}
					for m in hashes['modeltest']["per"].keys():
						hashes['tree_inference']["per"][c][i][m] = {}
						for t in hashes['filter-align']["per"].keys():
							hashes['tree_inference']["per"][c][i][m][t] = {}
							for a in hashes['align']["per"].keys():
								hashes['tree_inference']["per"][c][i][m][t][a] = get_hash(hashes['modeltest']["per"][m][t][a], "genetree_filtering,bootstrap_cutoff,"+c+" njtree,options,"+i, configfi)

		print("Gathered hashes until 'njtree'")
		print(hashes['tree_inference'])
		return hashes

def trigger(current_yaml, config_yaml, optional=[]):
	#the function triggers the start of a submodule if parameter files aren't there or when there are differences to the config file
	if not os.path.isfile(current_yaml):
		# that's the simplest case - output file is not yet there
		print("File "+current_yaml+" is not there")
		if optional:
			lis = [config_yaml] + optional
			return lis
		else:
			return config_yaml
	else:
		# that's a bit more tricky
		# could be the relevant yaml file has already been written, but there was a change made in the config file
		# which results in a newer time stamp so the rule would be rerun, but
		# could be that the file has changed but not in the parts that are relevant for the current step.
		# in this case we don't want it to be triggered so we return the current yaml file, which does not have a newer time stamp, so nothing happens
		print("File "+current_yaml+" is there")
		return current_yaml

		#this part actually reads in the two yaml files and compares the relevant parts
		# dont' think this is needed, since if the relevant parts have changed also the hash would change which would trigger a new run
#		print("Reading in file")
#		with open(current_yaml) as f:
#			my_dict = yaml.safe_load(f)
#		print(str(my_dict))
#		if find_top(my_dict, config):
#			print("ALL GOOD")
#			#in this case return the current yaml file. this was written before so it should not trigger the read_params rule
#			return current_yaml
#		else:
#			return config_yaml

def compare(yaml_one, yaml_two, optional=[]):
        print("Comparing:", yaml_one, "vs. ", yaml_two)
        return [trigger(yaml_one, yaml_two, optional=optional)]

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
