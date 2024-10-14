import os
import yaml
import time

def now(): #defined here again so that we don't have to import the whole hashing module
        return time.strftime("%Y-%m-%d %H:%M") + " -"

def parse_config_file(cf, debug=False):
	if debug:
		print(now(), "Will try to load yaml file:", cf)
	if not os.path.isfile(cf):
		print(now(), "Specified config file not found:", cf)
		sys.exit(1)
	with open(cf) as f:
		data = yaml.load(f, Loader=yaml.FullLoader)
		for key in ["align", "trimming", "modeltest"]: # check if all entries in different categories also have an options field in the yaml file.
			if key in data.keys():
				for k in data[key]["method"]:
					if k not in data[key]["options"].keys():
						data[key]["options"][k] = ""
	return data
