#!/usr/bin/env python3
import sys
import argparse
import os, io
import subprocess

if sys.version_info[0] < 3:
	raise Exception("Must be using Python 3")
	exit(1)

pars = argparse.ArgumentParser(prog="parse-log.py", description = """This script will parse phylociraptor logfiles.""", epilog = """written by Philipp Resl""")
pars.add_argument('-f', '--logfile', dest="logfile", help="Phylociraptor logile path.")
pars.add_argument('-l', '--list', dest="llist", action='store_true', default=False, help="List job IDs") 
pars.add_argument('-c', '--cancel', dest="cancel", action='store_true', default=False, help="Cancel Jobs based on IDs.") 
pars.add_argument('-v', '--verbose', action='store_true', default=False, help="More output")
args=pars.parse_args()

logfile = args.logfile 

print("Extracting job IDs from logfile:", logfile)
with open(logfile, "r") as fh:
	jobids = []
	subsys = ""
	for line in fh:
		if line.startswith("PHYLOCIRAPTOR"):
			runmode = line.split(" ")[3]
			print("According to log file the phylociraptor runmode was:", runmode)
		if line.startswith("EXECUTED"):
			if "--cluster" in line:
				print("Cluster jobs submission was used.")
			else:
				print("Local job execution was used.")
		if line.startswith("Submission system") and not subsys:
			subsys = line.strip().split(":")[-1].split(" ")[-1]
		if line.startswith("Submitted job"):
			jobid = [id for id in line.split("'")[1::2][0].split(" ") if id.isnumeric()][0]
			jobids.append(jobid)
if len(jobids) == 0:
	print("Found no job IDs in logfile. Will stop")
	sys.exit(0)
print("Found", len(jobids), "job IDs.")
print("Submission system:", subsys)
if subsys == "sge" or subsys == "torque":
	if args.llist:
		print("Job IDs:\n", " ".join(jobids))
	if args.cancel:
		print("Will cancel jobs using qdel:")
		cmd = ["qdel"] + jobids
		print(" ".join(cmd))
		popen = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
		if args.verbose:
			for line in io.TextIOWrapper(popen.stdout, encoding="utf-8").readlines():
				print(line.rstrip())
	else:
		print("Will not cancel any jobs. If you really like to cancel the jobs add -c or --cancel to your phylociraptor command.")
		
elif subsys == "slurm":
	if args.llist:
		print("Job IDs:\n", " ".join(jobids))
	if args.cancel:
		print("Will cancel jobs using scancel:")
		cmd = ["scancel"] + jobids
		print(" ".join(cmd))
		popen = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
		if args.verbose:
			print("INFO: scancel typically does not give any output.")
			for line in io.TextIOWrapper(popen.stdout, encoding="utf-8").readlines():
				print(line.rstrip())
	else:
		print("Will not cancel any jobs. If you really like to cancel the jobs add -c or --cancel to your phylociraptor command.")
