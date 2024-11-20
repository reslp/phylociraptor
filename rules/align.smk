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
hashes = collect_hashes("align", config, configfi, wd=os.getcwd())
current_hash = hashes["align"]["global"]
aligner_hash = hashes["align"]["per"]
#previous hash
previous_hash = hashes['filter-orthology']["global"]


BUSCOS, = glob_wildcards("results/orthology/single-copy-orthologs_deduplicated." + hashes["filter-orthology"]["global"] + "/{busco}_all.fas")

def determine_concurrency_limit():
	fname = "results/orthology/orthology/single-copy-orthologs_deduplicated."+previous_hash
	if os.path.isdir(fname):
		ngenes = glob.glob("results/orthology/single-copy-orthologs_deduplicated." + previous_hash + "/*.fas")
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
		previous = "results/orthology/busco/parameters.filter-orthology."+hashes["filter-orthology"]["global"]+".yaml"
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
		previous = "results/orthology/busco/parameters.filter-orthology."+hashes["filter-orthology"]["global"]+".yaml"
	output:
		"results/alignments/full/parameters.align."+current_hash+".yaml"
	shell:
		"""
		bin/read_write_yaml.py {input.trigger} {output} alignment,options alignment,method
		cat {input.previous} >> {output}
		"""


def per_aligner(wildcards):
	lis = []
	for b in BUSCOS:
		lis.append("results/alignments/full/"+wildcards.aligner+"."+aligner_hash[wildcards.aligner]+"/"+str(b)+"_aligned.fas")
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
def get_aligner_params(wildcards):
	if wildcards.aligner in config["alignment"]["options"].keys():
		return config["alignment"]["options"][wildcards.aligner]
	else:
		return ""


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
		aligner_params = get_aligner_params,
		mafft_alignment_params = config["alignment"]["options"]["mafft"],
		clustalo_alignment_params = config["alignment"]["options"]["clustalo"],
		muscle_alignment_params = config["alignment"]["options"]["muscle"],
		pars_sites = config["trimming"]["min_parsimony_sites"],
		nbatches = config["concurrency"],
		set = config["orthology"]["busco_options"]["set"],
		orthology_hash = hashes['filter-orthology']["global"],
		mode = config["orthology"]["method"]
	log:	"log/align/{aligner}_{batch}_get_aligment_statistics.{hash}.txt"
	singularity: containers["concat"] 
	shadow: "minimal"
	shell:
		"""
		# here the ids for the alignments need to be filtered as well first. maybe this can be changed in the concat.py script, so that an id file is not needed anymore.
		if [[ "{params.mode}" == "orthofinder" ]]; then
			concat.py -i $(ls -1 {params.wd}/results/alignments/full/{wildcards.aligner}.{wildcards.hash}/* | sed -n '{wildcards.batch}~{params.nbatches}p' | tr '\\n' ' ') --runmode concat -o results/statistics/ --biopython --statistics --seqtype {params.datatype} --noseq 2>&1 | tee {log}
		
		else
			concat.py -i $(ls -1 {params.wd}/results/alignments/full/{wildcards.aligner}.{wildcards.hash}/* | sed -n '{wildcards.batch}~{params.nbatches}p' | tr '\\n' ' ') -t <(for name in $(ls -1 {params.wd}/results/orthology/busco/busco_runs.{params.set}.{params.orthology_hash}); do echo "${{name%.*}}"; done) --runmode concat -o results/statistics/ --biopython --statistics --seqtype {params.datatype} --noseq 2>&1 | tee {log}

		fi
		mv results/statistics/statistics.txt {output.statistics_alignment}
		# make this output tab delimited so it is easier to parse
		ovstats="{params.alignment_method}\t{params.aligner_params}"
		ovstats="${{ovstats}}\t{params.pars_sites}"
		echo -e $ovstats > {output.overview_statistics}
		"""

def pull(wildcards):
	lis = []
	for ali in config["alignment"]["method"]:
		for i in range(1, config["concurrency"] + 1):
			lis.append("results/statistics/align-" + str(ali) + "." + str(aligner_hash[ali]) + "/" + str(ali) + "_statistics_alignments-" + str(i) + "-" + str(config["concurrency"]) + ".txt")
		lis.append("results/alignments/full/" + str(ali) + "." + str(aligner_hash[ali]) + "/parameters.align." + str(ali) + "." + str(aligner_hash[ali]) + ".yaml")
	#print(lis, file=sys.stderr)
	return lis	

rule align:
	input:
		pull,
		rules.read_params_global.output,
		checkpoint = "results/checkpoints/modes/filter_orthology."+previous_hash+".done"
	output:
		"results/checkpoints/modes/align."+current_hash+".done"
	shell:
		"""
		touch {output}
		echo "$(date) - phylociraptor align done." >> results/statistics/runlog.txt
		"""

for aligner in config["alignment"]["method"]:
	include: "aligners/" + aligner + ".smk"
