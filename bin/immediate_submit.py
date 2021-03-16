#!/usr/bin/env python3

import os
import sys
import shutil

from snakemake.utils import read_job_properties

# last command-line argument is the job script
jobscript = sys.argv[-1]

# last but one argument is the submission system
subs = sys.argv[-2]
#print("---------------------------------------------------------------------", file=sys.stderr)
print("Command line arguments are: ", file=sys.stderr)
print(sys.argv, file=sys.stderr)
# check command-line arguments for dependencies
dependencies = ""
if subs == "slurm":
	dependencies = set(sys.argv[1:-2])
elif subs == "sge":
	dependencies = dependencies.join(sys.argv[1:-2])
elif subs == "torque":
	dependencies = dependencies.join(sys.argv[1:-2])
else:
	print("Cannot get dependencies for submission system")
	sys.exit(1)
# parse the job script for the job properties that are encoded by snakemake within
# this also includes the information contained in the cluster-config file as job_properties["cluster"]
job_properties = read_job_properties(jobscript)

print("Job submision for job:", job_properties["rule"], file=sys.stderr)
print("Submission system:", subs,file=sys.stderr)

#some output for debugging:
#print(job_properties, file=sys.stderr)
#print(job_properties["wildcards"], file=sys.stderr)
#print("Dependencies:", file=sys.stderr)
#print(dependencies, file=sys.stderr)
#print(sys.argv, file=sys.stderr)

#list of items recognized by the -l (resource) flag of qsub:
sge_resources=["nodes", "pmem", "vmem", "mem", "walltime", "h_vmem", "ppn"] #my need to be extended for different clusters


cmdline=[]
# create list with command line arguments
if subs == "slurm":
	cmdline = ["sbatch"]
	if "species" in job_properties["wildcards"]:
                job_properties["cluster"]["J"] = job_properties["cluster"]["J"]+"-"+job_properties["wildcards"]["species"]
                prefix = job_properties["wildcards"]["species"] + "-" + job_properties["rule"] + "-slurm"
                job_properties["cluster"]["output"] = job_properties["cluster"]["output"].replace("slurm", prefix)
                job_properties["cluster"]["error"] = job_properties["cluster"]["error"].replace("slurm", prefix)
	elif "busco" in job_properties["wildcards"]:
                job_properties["cluster"]["J"] = job_properties["cluster"]["J"]+"-"+job_properties["wildcards"]["busco"]
                prefix = job_properties["wildcards"]["busco"] + "-" + job_properties["rule"] + "-slurm"
                job_properties["cluster"]["output"] = job_properties["cluster"]["output"].replace("slurm", prefix)
                job_properties["cluster"]["error"] = job_properties["cluster"]["error"].replace("slurm", prefix)
	else:		
		#properly assign job names for rules which don't have a sample wildcard
		job_properties["cluster"]["output"] = job_properties["cluster"]["output"].replace("slurm", job_properties["rule"])
		job_properties["cluster"]["error"] = job_properties["cluster"]["error"].replace("slurm", job_properties["rule"])
	
	#determine threads from the Snakemake profile, i.e. as determined in the Snakefile and the main config file respectively
	job_properties["cluster"]["ntasks"] = job_properties["threads"]
	job_properties["cluster"]["ntasks-per-node"] = job_properties["threads"]
	
	
	# create string for slurm submit options for rule
	# this accounts for -n (shared jobs) or -N (whole node jobs). -N overrules -n.
	if "N" in job_properties["cluster"]:
		slurm_args = "--partition={partition} --time={time} --qos={qos} --ntasks={ntasks} --ntasks-per-node={ntasks-per-node} --hint={hint} --output={output} --error={error} -N {N} -J {J} --mem={mem}".format(**job_properties["cluster"])
	else:
		slurm_args = "--partition={partition} --time={time} --qos={qos} --ntasks={ntasks} --ntasks-per-node={ntasks-per-node} --hint={hint} --output={output} --error={error} -n {n} -J {J} --mem={mem}".format(**job_properties["cluster"])
	cmdline.append(slurm_args)
	
	# now work on dependencies
	print("\nDependencies for this job:", file=sys.stderr)
	if dependencies:
		cmdline.append("--dependency")
		# only keep numbers (which are the jobids) in dependencies list. this is necessary because slurm returns more than the jobid. For other schedulers this could be different!
		dependencies = [x for x in dependencies if x.isdigit()]
		cmdline.append("afterok:" + ",".join(dependencies))
		print(dependencies, file=sys.stderr)
	else:
		print("none", file=sys.stderr)	
elif subs == "sge":
	#print("Job properties:", file=sys.stderr)
	#print(job_properties["cluster"], file=sys.stderr)
	# set informative job name:
	if "species" in job_properties["wildcards"]:
		name =  job_properties["cluster"]["N"] + "-" + job_properties["wildcards"]["species"]
	elif "busco" in job_properties["wildcards"]:
		name =  job_properties["cluster"]["N"] + "-" + job_properties["wildcards"]["busco"]
	else:
		name = job_properties["cluster"]["N"]
	job_properties["cluster"]["N"] = name
	# set name for cluster log files:
	prefix = "comparative-" + job_properties["rule"] + "-sge"
	job_properties["cluster"]["o"] = job_properties["cluster"]["o"].replace("sge", prefix).replace("%j",name)
	job_properties["cluster"]["e"] = job_properties["cluster"]["e"].replace("sge", prefix).replace("%j",name)
	cmdline = ["qsub"]
	
	# extract info an requested resources
	sge_args = "-cwd -V"
	
	if "pe" in job_properties["cluster"] and "nodes" in job_properties["cluster"]:
		print("WARNING: You have set pe and nodes in your cluster config file. Is this really what you want?", file=sys.stderr)

	for element in job_properties["cluster"]:
		#print(element, file=sys.stderr)
		
		if element in sge_resources:
			if "nodes" == element: #special case when nodes are specified, in this case threads are ppns
				sge_args += " -l %s=%s:ppn=%s" % ("nodes", job_properties["cluster"]["nodes"], job_properties["threads"])
			else:			
				sge_args += " -l %s=%s" % (element, job_properties["cluster"][element])
		else:
			if element == "pe": # set parallel environment and correct number of threads according to config.yaml
				sge_args += " -pe %s %s" % (job_properties["cluster"]["pe"], job_properties["threads"])
			else:
				sge_args += " -%s %s" % (element, job_properties["cluster"][element])
	#job_properties["sge_resources"] = resources
	# TODO: add correct thread handling for SGE clusters
	#sge_args = "-cwd -V -q {queue} -l h_vmem={mem} -pe {pe} {ntasks} -o {output} -e {error} -N {N}".format(**job_properties["cluster"])	
	cmdline.append(sge_args)

	#now work on dependencies
	print("\nDependencies for this job:", file=sys.stderr)
	if dependencies:
		cmdline.append("-hold_jid")
		# only keep numbers (which are the jobids) in dependencies list. this is necessary because slurm returns more than the jobid. For other schedulers this could be different!
		dependencies = [x for x in dependencies.split(" ") if x.isdigit()]
		cmdline.append(",".join(dependencies))
		print(dependencies, file=sys.stderr)
	else:
		print("none", file=sys.stderr)	
	# add @:
	#cmdline.append("-@")
elif subs == "torque":
	#print("Job properties:", file=sys.stderr)
	#print(job_properties["cluster"], file=sys.stderr)
	# set informative job name:
	if "species" in job_properties["wildcards"]:
		name =  job_properties["cluster"]["N"] + "-" + job_properties["wildcards"]["species"]
	elif "busco" in job_properties["wildcards"]:
		name =  job_properties["cluster"]["N"] + "-" + job_properties["wildcards"]["busco"]
	else:
		name = job_properties["cluster"]["N"]
	job_properties["cluster"]["N"] = name
	# set name for cluster log files:
	prefix =  job_properties["rule"] + "-torque"
	job_properties["cluster"]["o"] = job_properties["cluster"]["o"].replace("torque", prefix).replace("%j",name)
	job_properties["cluster"]["e"] = job_properties["cluster"]["e"].replace("torque", prefix).replace("%j",name)
	cmdline = ["qsub"]
	
	# extract info an requested resources
	sge_args = "-w $(pwd) -V"
	
	if "pe" in job_properties["cluster"] and "nodes" in job_properties["cluster"]:
		print("WARNING: You have set pe and nodes in your cluster config file. Is this really what you want?", file=sys.stderr)

	for element in job_properties["cluster"]:
		#print(element, file=sys.stderr)
		
		if element in sge_resources:
			if "nodes" == element: #special case when nodes are specified, in this case threads are ppns
				sge_args += " -l %s=%s:ppn=%s" % ("nodes", job_properties["cluster"]["nodes"], job_properties["threads"])
			else:			
				sge_args += " -l %s=%s" % (element, job_properties["cluster"][element])
		else:
			if element == "pe": # set parallel environment and correct number of threads according to config.yaml
				sge_args += " -pe %s %s" % (job_properties["cluster"]["pe"], job_properties["threads"])
			else:
				sge_args += " -%s %s" % (element, job_properties["cluster"][element])
	#job_properties["sge_resources"] = resources
	cmdline.append(sge_args)

	#now work on dependencies
	print("\nDependencies for this job:", file=sys.stderr)
	print(dependencies, file=sys.stderr)
	if dependencies:
		deps = "-W depend=afterok:"
		#cmdline.append("-W depend=afterok")
		# only keep numbers (which are the jobids) in dependencies list. this is necessary because slurm returns more than the jobid. For other schedulers this could be different!
		dependencies = [x for x in dependencies.split(" ") if x.isdigit()]
		deps += ":".join(dependencies)
		cmdline.append(deps)
		print(dependencies, file=sys.stderr)
		print(deps, file=sys.stderr)
	else:
		print("none", file=sys.stderr)	
	# add @:
	#cmdline.append("-@")
else:
	#print("Immediate submit error: Unkown submission system!")
	sys.exit(1)


cmdline.append(jobscript)


#now write final commandback to the system
print("\nSubmission command:", file=sys.stderr)
print(" ".join(cmdline), file=sys.stderr)
print
os.system(" ".join(cmdline))

print("---------------------------------------------------------------------", file=sys.stderr)
