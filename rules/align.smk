include: "functions.smk"
import subprocess
import glob
import yaml
import sys
import hashlib

# get list of containers to use:
with open("data/containers.yaml", "r") as yaml_stream:
    containers = yaml.safe_load(yaml_stream)

#create new hashes for current stage (alignment) by combining the previuos hash with a newly generated
hashes = collect_hashes("align")
current_hash = hashes["align"]["global"]
mafft_hash = hashes["align"]["per"]["mafft"]
muscle_hash = hashes["align"]["per"]["muscle"]
clustalo_hash = hashes["align"]["per"]["clustalo"]
#previous hash
previous_hash = hashes['filter-orthology']


BUSCOS, = glob_wildcards("results/orthology/busco/busco_sequences_deduplicated."+hashes["filter-orthology"]+"/{busco}_all.fas")

def determine_concurrency_limit():
	fname = "results/orthology/busco/busco_sequences_deduplicated."+previous_hash
	if os.path.isdir(fname):
		ngenes = glob.glob("results/orthology/busco/busco_sequences_deduplicated"+"/*.fas")
		ngenes = len(ngenes)
		if ngenes < config["concurrency"]:
			return range(1, ngenes + 1)
		else:
			return range(1, config["concurrency"] + 1)
	else:
		return

batches = determine_concurrency_limit()

aligners = get_aligners()
def compare_align(wildcards):
	return [trigger("results/alignments/full/{aligner}.{hash}/parameters.align.{aligner}.{hash}.yaml".format(aligner=wildcards.aligner, hash=wildcards.hash), configfi)]	

rule read_params_per:
	input:
		trigger = compare_align,
		previous = "results/orthology/busco/params.filter-orthology."+hashes["filter-orthology"]+".yaml"
	output:
		"results/alignments/full/{aligner}.{hash}/parameters.align.{aligner}.{hash}.yaml"
	shell:
		"""
		bin/read_write_yaml.py {input.trigger} {output} alignment,options,{wildcards.aligner}
		cat {input.previous} >> {output}
		"""

rule read_params_global:
	input:
		trigger = compare("results/alignments/full/parameters.align."+current_hash+".yaml", configfi),
		previous = "results/orthology/busco/params.filter-orthology."+hashes["filter-orthology"]+".yaml"
	output:
		"results/alignments/full/parameters.align."+current_hash+".yaml"
	shell:
		"""
		bin/read_write_yaml.py {input.trigger} {output} alignment,options alignment,method
		cat {input.previous} >> {output}
		"""

rule clustalo:
		input:
			"results/alignments/full/clustalo."+clustalo_hash+"/parameters.align.clustalo."+clustalo_hash+".yaml",
			sequence_file = "results/orthology/busco/busco_sequences_deduplicated."+hashes["filter-orthology"]+"/{busco}_all.fas",
		output:
			alignment = "results/alignments/full/clustalo."+clustalo_hash+"/{busco}_aligned.fas",
		benchmark:
			"results/statistics/benchmarks/align/clustalo_align_{busco}."+clustalo_hash+".txt"
		log:
			"log/align/clustalo/clustalo_align_{busco}."+clustalo_hash+".log.txt"
		singularity:
			containers["clustalo"]
		threads:
			int(config["alignment"]["threads"])
		params:
			config["alignment"]["options"]["clustalo"]
		shell:
			"""
			clustalo -i {input.sequence_file} --threads={threads} $(if [[ "{params}" != "None" ]]; then echo {params}; fi) 1> {output.alignment} 2> {log}
			"""

rule mafft:
		input:
			"results/alignments/full/mafft."+mafft_hash+"/parameters.align.mafft."+mafft_hash+".yaml",
			sequence_file = "results/orthology/busco/busco_sequences_deduplicated."+hashes["filter-orthology"]+"/{busco}_all.fas",
		output:
			alignment = "results/alignments/full/mafft."+mafft_hash+"/{busco}_aligned.fas",
		benchmark:
			"results/statistics/benchmarks/align/mafft_align_{busco}."+mafft_hash+".txt"
		log:
			"log/align/mafft/mafft_align_{busco}.log.txt"
		singularity:
			containers["mafft"]
		threads:
			int(config["alignment"]["threads"])
		params:
			config["alignment"]["options"]["mafft"]
		shell:
			"""
			mafft --thread {threads} $(if [[ "{params}" != "None" ]]; then echo {params}; fi) {input.sequence_file} 1> {output.alignment} 2> {log}
			"""
rule muscle:
		input:
			"results/alignments/full/muscle."+muscle_hash+"/parameters.align.muscle."+muscle_hash+".yaml",
			sequence_file = "results/orthology/busco/busco_sequences_deduplicated."+hashes["filter-orthology"]+"/{busco}_all.fas",
		output:
			alignment = "results/alignments/full/muscle."+muscle_hash+"/{busco}_aligned.fas",
		benchmark:
			"results/statistics/benchmarks/align/muscle_align_{busco}."+muscle_hash+".txt"
		log:
			"log/align/muscle/muscle_align_{busco}."+muscle_hash+".log.txt"
		singularity:
			containers["muscle"]
		threads:
			int(config["alignment"]["threads"])
		params:
			config["alignment"]["options"]["muscle"]
		shell:
			"""
			muscle -super5 {input.sequence_file} -threads {threads} $(if [[ "{params}" != "None" ]]; then echo {params}; fi) -output {output.alignment} 2> {log}
			"""

def per_aligner(wildcards):
	lis = []
#	print(wildcards)
	dict = {"muscle": muscle_hash, "mafft": mafft_hash, "clustalo": clustalo_hash}
	for b in BUSCOS:
		lis.append("results/alignments/full/"+wildcards.aligner+"."+dict[wildcards.aligner]+"/"+str(b)+"_aligned.fas")
	return lis

rule aggregate_alignments:
	input:
		per_aligner
	output:
		checkpoint = "results/checkpoints/{aligner}_aggregate_align.{hash}.done"
	shell:
		"""
		touch {output.checkpoint}
		"""

rule get_alignment_statistics:
	input:
		rules.aggregate_alignments.output.checkpoint
	output:
		statistics_alignment = "results/statistics/align-{aligner}.{hash}/{aligner}_statistics_alignments-{batch}-"+str(config["concurrency"])+".txt",
		overview_statistics = "results/statistics/align-{aligner}.{hash}/{aligner}_overview-{batch}-"+str(config["concurrency"])+".txt"
	params:
		wd = os.getcwd(),
		ids = config["species"],
		datatype = config["filtering"]["seq_type"],
		alignment_method = "{aligner}",
		mafft_alignment_params = config["alignment"]["options"]["mafft"],
		clustalo_alignment_params = config["alignment"]["options"]["clustalo"],
		muscle_alignment_params = config["alignment"]["options"]["muscle"],
		pars_sites = config["trimming"]["min_parsimony_sites"],
		nbatches = config["concurrency"],
		set = config["orthology"]["busco_options"]["set"],
		orthology_hash = hashes['orthology']
	log:	"log/align/{aligner}_{batch}_get_aligment_statistics.{hash}.txt"
	singularity: containers["concat"] 
	shadow: "minimal"
	shell:
		"""
		# here the ids for the alignments need to be filtered as well first. maybe this can be changed in the concat.py script, so that an id file is not needed anymore.
		concat.py -i $(ls -1 {params.wd}/results/alignments/full/{wildcards.aligner}.{wildcards.hash}/* | sed -n '{wildcards.batch}~{params.nbatches}p' | tr '\\n' ' ') -t <(for name in $(ls -1 {params.wd}/results/orthology/busco/busco_runs.{params.set}.{params.orthology_hash}); do echo "${{name%.*}}"; done) --runmode concat -o results/statistics/ --biopython --statistics --seqtype {params.datatype} --noseq 2>&1 | tee {log}
		mv results/statistics/statistics.txt {output.statistics_alignment}
		# make this output tab delimited so it is easier to parse
		ovstats="{params.alignment_method}"
		if [[ "{wildcards.aligner}" == "mafft" ]]; then 
			ovstats="${{ovstats}}\t{params.mafft_alignment_params}"
		fi
		if [[ "{wildcards.aligner}" == "clustalo" ]]; then 
			ovstats="${{ovstats}}\t{params.clustalo_alignment_params}"
		fi
		if [[ "{wildcards.aligner}" == "muscle" ]]; then 
			ovstats="${{ovstats}}\t{params.muscle_alignment_params}"
		fi
		ovstats="${{ovstats}}\t{params.pars_sites}"
		echo -e $ovstats > {output.overview_statistics}
		"""

def pull(wildcards):
	lis = []
	dict = {"muscle": muscle_hash, "mafft": mafft_hash, "clustalo": clustalo_hash}
	for ali in config["alignment"]["method"]:
		for i in range(1, config["concurrency"] + 1):
			lis.append("results/statistics/align-"+str(ali)+"."+str(dict[ali])+"/"+str(ali)+"_statistics_alignments-"+str(i)+"-"+str(config["concurrency"])+".txt")
		lis.append("results/alignments/full/"+str(ali)+"."+str(dict[ali])+"/parameters.align."+str(ali)+"."+str(dict[ali])+".yaml")
	return lis	

rule align:
	input:
		pull,
		rules.read_params_global.output,
		checkpoint = "results/checkpoints/modes/filter_orthology."+previous_hash+".done",
	output:
		"results/checkpoints/modes/align."+current_hash+".done"
	shell:
		"""
		touch {output}
		echo "$(date) - phylociraptor align done." >> results/statistics/runlog.txt
		"""
