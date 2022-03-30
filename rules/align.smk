include: "functions.smk"
import subprocess
import glob
import yaml
import sys

# get list of containers to use:
with open("data/containers.yaml", "r") as yaml_stream:
    containers = yaml.safe_load(yaml_stream)

BUSCOS, = glob_wildcards("results/orthology/busco/busco_sequences_deduplicated/{busco}_all.fas")

def determine_concurrency_limit():
	fname = "results/orthology/busco/busco_sequences_deduplicated"
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

rule clustalo:
		input:
			checkpoint = "results/checkpoints/modes/filter_orthology.done",
			sequence_file = "results/orthology/busco/busco_sequences_deduplicated/{busco}_all.fas"
		output:
			alignment = "results/alignments/full/clustalo/{busco}_aligned.fas",
		benchmark:
			"results/statistics/benchmarks/align/clustalo_align_{busco}.txt"
		log:
			"results/statistics/benchmarks/align/clustalo_align_{busco}.log.txt"
		singularity:
			containers["clustalo"]
		threads:
			int(config["alignment"]["threads"])
		params:
			config["alignment"]["clustalo_parameters"]
		shell:
			"""
			clustalo -i {input.sequence_file} --threads={threads} $(if [[ "{params}" != "None" ]]; then echo {params}; fi) 1> {output.alignment} 2> {log}
			"""

rule mafft:
		input:
			checkpoint = "results/checkpoints/modes/filter_orthology.done",
			sequence_file = "results/orthology/busco/busco_sequences_deduplicated/{busco}_all.fas"
		output:
			alignment = "results/alignments/full/mafft/{busco}_aligned.fas",
		benchmark:
			"results/statistics/benchmarks/align/mafft_align_{busco}.txt"
		log:
			"results/statistics/benchmarks/align/mafft_align_{busco}.log.txt"
		singularity:
			containers["mafft"]
		threads:
			int(config["alignment"]["threads"])
		params:
			config["alignment"]["mafft_parameters"]
		shell:
			"""
			mafft --thread {threads} $(if [[ "{params}" != "None" ]]; then echo {params}; fi) {input.sequence_file} 1> {output.alignment} 2> {log}
			"""
rule muscle:
		input:
			checkpoint = "results/checkpoints/modes/filter_orthology.done",
			sequence_file = "results/orthology/busco/busco_sequences_deduplicated/{busco}_all.fas"
		output:
			alignment = "results/alignments/full/muscle/{busco}_aligned.fas",
		benchmark:
			"results/statistics/benchmarks/align/muscle_align_{busco}.txt"
		log:
			"results/statistics/benchmarks/align/muscle_align_{busco}.log.txt"
		singularity:
			containers["muscle"]
		threads:
			int(config["alignment"]["threads"])
		params:
			config["alignment"]["muscle_parameters"]
		shell:
			"""
			muscle -super5 {input.sequence_file} -threads {threads} $(if [[ "{params}" != "None" ]]; then echo {params}; fi) -output {output.alignment} 2> {log}
			"""

rule aggregate_alignments:
	input:
		expand("results/alignments/full/{{aligner}}/{bus}_aligned.fas", aligner=config["alignment"]["method"], bus=BUSCOS)
	output:
		checkpoint = "results/checkpoints/{aligner}_aggregate_align.done"
	shell:
		"""
		touch {output.checkpoint}
		"""
rule get_alignment_statistics:
	input:
		rules.aggregate_alignments.output.checkpoint
	output:
		statistics_alignment = "results/statistics/align-{aligner}/{aligner}_statistics_alignments-{batch}-"+str(config["concurrency"])+".txt",
		overview_statistics = "results/statistics/align-{aligner}/{aligner}_overview-{batch}-"+str(config["concurrency"])+".txt"
	params:
		wd = os.getcwd(),
		ids = config["species"],
		datatype = config["filtering"]["seq_type"],
		alignment_method = "{aligner}",
		mafft_alignment_params = config["alignment"]["mafft_parameters"],
		clustalo_alignment_params = config["alignment"]["clustalo_parameters"],
		pars_sites = config["filtering"]["min_parsimony_sites"],
		nbatches = config["concurrency"],
	singularity: containers["concat"] 
	shadow: "minimal"
	shell:
		"""
		# here the ids for the alignments need to be filtered as well first. maybe this can be changed in the concat.py script, so that an id file is not needed anymore.
		concat.py -i $(ls -1 {params.wd}/results/alignments/full/{wildcards.aligner}/* | sed -n '{wildcards.batch}~{params.nbatches}p' | tr '\\n' ' ') -t <(ls -1 {params.wd}/results/orthology/busco/busco_runs/) --runmode concat -o results/statistics/ --biopython --statistics --seqtype {params.datatype} --noseq
		mv results/statistics/statistics.txt {output.statistics_alignment}
		# make this output tab delimited so it is easier to parse
		ovstats="{params.alignment_method}"
		if [[ "{wildcards.aligner}" == "mafft" ]]; then 
			ovstats="${{ovstats}}\t{params.mafft_alignment_params}"
		fi
		if [[ "{wildcards.aligner}" == "clustalo" ]]; then 
			ovstats="${{ovstats}}\t{params.clustalo_alignment_params}"
		fi
		ovstats="${{ovstats}}\t{params.pars_sites}"
		echo -e $ovstats > {output.overview_statistics}
		"""
rule align:
	input:
#		expand("results/alignments/full/{aligner}/{busco}_aligned.fas", aligner=config["alignment"]["method"], busco=BUSCOS),
#		expand("results/checkpoints/{aligner}_aggregate_align.done", aligner=config["alignment"]["method"]),
		expand("results/statistics/align-{aligner}/{aligner}_statistics_alignments-{batch}-"+str(config["concurrency"])+".txt", aligner=aligners, batch=batches)
	output:
		"results/checkpoints/modes/align.done"
	shell:
		"""
		touch {output}
		echo "$(date) - Pipeline part 2 (align) done." >> results/statistics/runlog.txt
		"""
