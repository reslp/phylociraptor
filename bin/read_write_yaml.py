#!/usr/bin/env python

import sys
import yaml
from libphylociraptor.filehandling import *

config = parse_config_file(sys.argv[1])

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
