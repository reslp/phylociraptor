#!/usr/bin/env python

import sys, os, io

if sys.version_info[0] < 3:
	raise Exception("Must be using Python 3")

import argparse
import subprocess


#some basic variables needed:
njobs = "10001"
latency_wait = "10" #in seconds
singularity_bindpoints = "-B $(pwd)/.usr_tmp:/usertmp"

default_help = """
Welcome to phylociraptor, the rapid phylogenomic tree calculator

Usage: phylociraptor <command> <arguments>

Commands:
	setup			Setup pipeline
	orthology		Infer orthologs in a set of genomes
	filter-orthology	Filter orthology results
	align			Create alignments for orthologous genes
	filter-align		Filter alignments
	model			Perform modeltesting
	tree			Calculate phylogenomic trees
	speciestree		Calculate species tree
	njtree			Calculate NJ tree

-v, --version Print version

"""

setup_help = """
phylociraptor setup - will prepare your analysis

Usage: phylociraptor setup <arguments>

Argumemts:
	-t, --cluster		Specify cluster type. Options: slurm, sge, torque. Default: None
	-c, --cluster-config	Specify Cluster config file path. Default: None
	-f, --force		Force run this part in case it has already bin run.
	
	--dry			Make a dry run.
	--verbose		Display more output.
	-h, --help		Display help.
	
"""

orthology_help = """
phylociraptor orthology - Will infer orthologous genes in a set of genomes.

Usage: phylociraptor orthology <arguments>

Argumemts:
        -t, --cluster           Specify cluster type. Options: slurm, sge, torque, serial. Default: serial.
        -c, --cluster-config    Specify Cluster config file path. Default: None
        -f, --force             Force run this part in case it has already bin run.
        
        --dry                   Make a dry run.
        --verbose               Display more output.
        -h, --help              Display help.
        
"""

forthology_help = """
phylociraptor filter-orthology - Will filter orthology results produced by phylociraptor orthology.

Usage: phylociraptor filter-orthology <arguments>

Argumemts:
        -t, --cluster           Specify cluster type. Options: slurm, sge, torque, serial. Default: serial
        -c, --cluster-config    Specify Cluster config file path. Default: None
        -f, --force             Force run this part in case it has already bin run.
        
        --dry                   Make a dry run.
        --verbose               Display more output.
        -h, --help              Display help.
        
"""
def help_message(mes):
	return mes

def determine_submission_mode(flag):
	cmd = []
	if "serial" in flag:
		return []
	else:
		return ["--cluster", 'bin/immediate_submit.py {dependencies} %s' % flag, "--immediate-submit", "--jobs", njobs, "--notemp"]

def get_flags(flags):
	mapdict ={
	#"t": '--cluster "bin/immediate_submit.py {dependencies} ', "cluster": '--cluster "bin/immediate_submit.py {dependencies} ',
	"c": "--cluster-config", "cluster_config": "--cluster-config",
	"force": "-f",
	"dry": "-n"
	}
	cmd = []
	for flag in flags.keys():
		#print(flag)
		if flag in mapdict.keys() and flags[flag] != None:
			if flag == "t" or flag == "cluster": #handle cluster specification
				arg = mapdict[flag]
				arg = arg + " "+flags[flag]+'"'
				cmd.append(arg) 
			if flag == "c" or flag == "cluster_config": #handle cluster config file
				#print("here")
				#arg = mapdict[flag]
				#arg = arg + " " + flags[flag]
				cmd.append(mapdict[flag])
				cmd.append(flags[flag])
			else:
				if flags[flag]:
					cmd.append(mapdict[flag])
	
	return cmd

def check_required_files(runmode):
	outfile_dict = {
	"setup": ["data/config.yaml"],
	"orthology": ["data/config.yaml", ".phylogenomics_setup.done"],
	"filter-orthology": [".phylogenomics_setup.done", "checkpoints/orthology.done"]
	}
	
	for f in outfile_dict[runmode]:
		if not os.path.isfile(f):
			return f
	return

def execute_command(cmd, verbose):
	popen = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
	line = ""
	jobcounts = False
	njobs = 0
	for line in io.TextIOWrapper(popen.stdout, encoding="utf-8"):
		if verbose:
			yield line[:-1]
			#if char != "\n":
			#	line += char
			#else:
			#	result = line
			#	line = ""
			#	yield result
				#sys.stdout.flush()
		else:
			if char != "\n":
				line += char
			else: #now we have a complete line of output
				result = line
				#print(result)
				if result.startswith("Job counts"):
					jobcounts = True
					#yield result
				elif jobcounts and result != "":
					if len(result.split("\t")) == 2: # the line with the total number of jobs has two elements
						yield "Total number of tasks to run: "+ result.split("\t")[1] 
				elif jobcounts and result == "":
					jobcounts = False
				line=""
	
pars = argparse.ArgumentParser(prog="phylociraptor", description = """Description: phylociraptor, the rapid phylogenomic tree calculator""", epilog = """written by Philipp Resl and Christoph Hahn""", usage=help_message(default_help))
pars.add_argument('command', action='store')
pars.add_argument('arguments', action='store', nargs=argparse.REMAINDER)
args = pars.parse_args()

if args.command == "setup":
	setup_parser = argparse.ArgumentParser(prog="phylociraptor setup", description = """Description: phylociraptor setup - Will help to setup your analysis.""", epilog = """written by Philipp Resl and Christoph Hahn""", usage=help_message(setup_help), add_help=False)	
	setup_parser.add_argument("-f", "--force", action="store_true", help="Force setup in case this has been run before")
	setup_parser.add_argument("-t", "--cluster",  action="store", help="Submission Possible options: slurm, sge, torque, serial (no job submission)", required=True)
	setup_parser.add_argument("-c", "--cluster-config", action="store", help="Cluster config file.")
	setup_parser.add_argument("--dry", action="store_true", help="Only a dry run")
	setup_parser.add_argument("-h", "--help", action="store_true")
	setup_parser.add_argument("--verbose", action="store_true", default=True)
	setup_args = setup_parser.parse_args(args.arguments)
	
	if setup_args.help:
		print(help_message(setup_help))
		sys.exit(0)
	
	print("Preparing to run phylociraptor setup...")
	cmd = ["snakemake", "--use-singularity", "-r", "setup", "--latency-wait", latency_wait]
	cmd += get_flags(vars(setup_args))
	cmd += determine_submission_mode(setup_args.cluster)
	print(cmd)

	
	#execute_command(cmd)
	for line in execute_command(cmd, setup_args.verbose):
		print(line)
	#print(p.stdout)
	#print(p.stderr)
if args.command=="orthology":
	orthology_parser = argparse.ArgumentParser(prog="phylociraptor orthology", description = """Description: phylociraptor orthology - Will infer orthologous genes for a set of genomes""", epilog = """written by Philipp Resl and Christoph Hahn""", add_help=False)	
	orthology_parser.add_argument("-f", "--force", action="store_true", help="Force setup in case this has been run before")
	orthology_parser.add_argument("-t", "--cluster",  action="store", help="Submission system. Possible options: slurm, sge, torque, serial (no submission; default).", default="serial")
	orthology_parser.add_argument("-c", "--cluster-config", action="store", help="Cluster config file.")
	orthology_parser.add_argument("--dry", action="store_true", help="Only a dry run")
	orthology_parser.add_argument("-h", "--help", action="store_true")
	orthology_parser.add_argument("--verbose", action="store_true", default=True)	
	orthology_args = orthology_parser.parse_args(args.arguments)
	#if len(sys.argv) < 3: # too few argument
	#	print(help_message(orthology_help))
	#	sys.exit(0)
	if orthology_args.help: # help is specified
		print(help_message(orthology_help))
		sys.exit(0)
	if check_required_files(args.command):
		print("Some files are missing preventing this part to run.\nDid your run phylociraptor setup already? Missing file:", check_required_files(args.command))
	#print(orthology_args)
	cmd = ["snakemake", "--use-singularity", "--singularity-args", "%s"% singularity_bindpoints,"-r", "orthology", "--latency-wait", latency_wait]
	cmd += determine_submission_mode(orthology_args.cluster)
	cmd += get_flags(vars(orthology_args))
	print(cmd)
	
	print("Preparing to run phylociraptor orthology...")
	for line in execute_command(cmd, orthology_args.verbose):
                print(line)
if args.command=="filter-orthology":
	forthology_parser = argparse.ArgumentParser(prog="phylociraptor filter-orthology", description = """Description: phylociraptor filter-orthology - Will filter results from phylociraptor orthology""", epilog = """written by Philipp Resl and Christoph Hahn""", add_help=False)
	forthology_parser.add_argument("-f", "--force", action="store_true", help="Force setup in case this has been run before")
	forthology_parser.add_argument("-t", "--cluster",  action="store", help="Submission system. Possible options: slurm, sge, torque, serial (no submission; default).", default="serial")
	forthology_parser.add_argument("-c", "--cluster-config", action="store", help="Cluster config file.")
	forthology_parser.add_argument("--dry", action="store_true", help="Only a dry run")
	forthology_parser.add_argument("-h", "--help", action="store_true")
	forthology_parser.add_argument("--verbose", action="store_true", default=True)
	forthology_args = forthology_parser.parse_args(args.arguments)	
	print(sys.argv)
	#if len(sys.argv) < 3: # too few argument
	#	print(help_message(forthology_help))
	#	sys.exit(0)
	if forthology_args.help: # help is specified
		print(help_message(forthology_help))
		sys.exit(0)
	if check_required_files(args.command):
		print("Some files are missing preventing this part to run.\nDid your run phylociraptor setup and orthology already? Missing file:", check_required_files(args.command))
		
