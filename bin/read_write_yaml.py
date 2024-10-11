#!/usr/bin/env python

import sys
import yaml

with open(sys.argv[1], "r") as yaml_stream:
	config = yaml.safe_load(yaml_stream)
	for key in ["align", "trimming", "modeltest"]: # check if all entries in different categories also have an options field in the yaml file.
		if key in config.keys():
			for k in config[key]["method"]:
				if k not in config[key]["options"].keys():
					config[key]["options"][k] = ""
#	for key in config:
#		print(key,config[key])

out_dict = {}

if len(sys.argv) > 2:
	for target in sys.argv[3:]:
		t = target.split(",")
#		print(str(t))
		if len(t) == 1:
			print("level 1")
			string = str(config[t[0]])
			out_dict[t[0]] = string
		if len(t) == 2:
			print("level 2")
#			print(config[t[0]][t[1]])
			if not t[0] in out_dict:
				out_dict[t[0]] = {t[1]: str(config[t[0]][t[1]])}
			elif not t[1] in out_dict[t[0]]:
				out_dict[t[0]][t[1]] = str(config[t[0]][t[1]])
		if len(t) == 3:
			print("level 3")
			if isinstance(config[t[0]][t[1]], list):
				string = t[2]
				if not t[0] in out_dict:
					out_dict[t[0]] = {t[1]: string}
				elif not t[1] in out_dict[t[0]]:
					out_dict[t[0]][t[1]] = string
			else:
				if not t[0] in out_dict:
					out_dict[t[0]] = {t[1]: {t[2]: str(config[t[0]][t[1]][t[2]])}}
				elif not t[1] in out_dict[t[0]]:
					out_dict[t[0]][t[1]] = {t[2]: str(config[t[0]][t[1]][t[2]])}
				elif not t[2] in out_dict[t[0]][t[1]]:
					out_dict[t[0]][t[1]][t[2]] = str(config[t[0]][t[1]][t[2]])

with open(sys.argv[2], 'w') as file:
	yaml.dump(out_dict, file)
