import sys

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
	trees = config["tree"]["method"]
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
