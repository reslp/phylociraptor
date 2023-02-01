include: "functions.smk"
import yaml
import sys
import hashlib


# get list of containers to use:
with open("data/containers.yaml", "r") as yaml_stream:
    containers = yaml.safe_load(yaml_stream)


#get hash for current
hashes = collect_hashes("filter-orthology", config, configfi, wd=os.getcwd())
current_hash = hashes["filter-orthology"]["global"]
previous_hash = hashes["orthology"]["global"]

def determine_concurrency_limit():
	fname = "results/orthology/orthology_table."+previous_hash+".txt"
	if os.path.isfile(fname):
		handle = open(fname, "r")
		ngenes = 0
		ngenes = len(handle.readline().split("\t")) - 2 #minus 2 because there are two extra columns
		if ngenes < config["concurrency"]:
			return range(1, ngenes + 1)
		else:
			return range(1, config["concurrency"] + 1)
	else:
		return

batches = determine_concurrency_limit()

rule read_params_global:
	input:
		trigger = compare("results/orthology/busco/parameters.filter-orthology."+current_hash+".yaml", configfi),
		previous = "results/orthology/busco/parameters.orthology."+previous_hash+".yaml"
	output:
		"results/orthology/busco/parameters.filter-orthology."+current_hash+".yaml"
	shell:
		"""
		bin/read_write_yaml.py {input.trigger} {output} filtering,dupseq filtering,cutoff filtering,minsp filtering,seq_type filtering,exclude_orthology
		cat {input.previous} >> {output}
		"""

rule create_sequence_files:
	input:
		table = "results/orthology/orthology_table."+previous_hash+".txt",
		params = rules.read_params_global.output
	output:
		sequence_dir=directory("results/orthology/busco/busco_sequences."+current_hash+"/{batch}-"+str(config["concurrency"])),
		checkpoint = "results/checkpoints/create_sequence_files_{batch}-"+str(config["concurrency"])+"."+current_hash+".done",
		genome_statistics = "results/statistics/orthology_filtering/orthology_filtering_genomes_statistics_{batch}-"+str(config["concurrency"])+"."+current_hash+".txt",
		gene_statistics = "results/statistics/orthology_filtering/orthology_filtering_gene_statistics_{batch}-"+str(config["concurrency"])+"."+current_hash+".txt"
	benchmark:
		"results/statistics/benchmarks/create_seq_files/create_sequence_files_{batch}-"+str(config["concurrency"])+"."+current_hash+".txt"
	params:
		cutoff=config["filtering"]["cutoff"],
		minsp=config["filtering"]["minsp"],
		busco_dir = "results/orthology/busco/busco_runs."+str(config["orthology"]["busco_options"]["set"])+"."+previous_hash,
		seqtype = config["filtering"]["seq_type"],
		exclude = config["filtering"]["exclude_orthology"],
		nbatches = config["concurrency"],
	singularity:
		containers["biopython"]
	log: "log/filter-orthology/exSeqfiles_{batch}-"+str(config["concurrency"])+"."+current_hash+".log"
	shell:
		""" 
		if [[ ! -d {output.sequence_dir} ]]; then mkdir -p {output.sequence_dir}; fi
		exclude=""
		if [[ "{params.exclude}" != "None" ]]
		then
			exclude="--exclude {params.exclude}"
		fi
		python bin/create_sequence_files.py --type {params.seqtype} --busco_table {input.table} --busco_results {params.busco_dir} --cutoff {params.cutoff} --outdir {output.sequence_dir} --minsp {params.minsp} --genome_statistics {output.genome_statistics} --gene_statistics {output.gene_statistics} --batchID {wildcards.batch} --nbatches {params.nbatches} --fix_aa_codes $exclude &> {log}
		touch {output.checkpoint}   
		"""

rule remove_duplicated_sequence_files:
	input:
		checkpoint = expand("results/checkpoints/create_sequence_files_{batch}-"+str(config["concurrency"])+"."+current_hash+".done", batch=batches)
	output:
		checkpoint = "results/checkpoints/remove_duplicated_sequence_files."+current_hash+".done",
		statistics = "results/statistics/duplicated_sequences_handling_information."+current_hash+".txt"
	benchmark:
		"results/statistics/benchmarks/remdupseq/remove_duplicated_sequence_files."+current_hash+".txt"
	singularity: containers["biopython"] 
	params:
		wd = os.getcwd(),
		minsp=config["filtering"]["minsp"],
		dupseq=config["filtering"]["dupseq"],
		current_hash = current_hash
	shell:
		"""
		if [[ -d results/orthology/busco/busco_sequences_deduplicated.{params.current_hash} ]]; then
			rm -rf results/orthology/busco/busco_sequences_deduplicated.{params.current_hash}
		fi
		mkdir results/orthology/busco/busco_sequences_deduplicated.{params.current_hash}
		
		rm -f {output.statistics}
		
		if [[ {params.dupseq} == "persample" ]];
		then
			echo "$(date) - BUSCO files will be filtered on a per-sample basis. This could lower the number of species in the final tree." >> results/statistics/runlog.txt
		else
			echo "$(date) - BUSCO files will be filtered on a per BUSCO gene basis. This could lower the number of genes used to calculate the final tree." >> results/statistics/runlog.txt
		fi
		
		for file in results/orthology/busco/busco_sequences.{params.current_hash}/*/*.fas;
			do
			if [[ {params.dupseq} == "persample" ]];
			then
				# per sequence filtering
				python bin/filter_alignments.py --alignments {params.wd}/$file --outdir "{params.wd}/results/orthology/busco/busco_sequences_deduplicated.{params.current_hash}" --per_sample --minsp {params.minsp} >> {output.statistics}
			else
				# whole alignment filtering
				python bin/filter_alignments.py --alignments {params.wd}/$file --outdir "{params.wd}/results/orthology/busco/busco_sequences_deduplicated.{params.current_hash}" >> {output.statistics}
			fi
		done
		#gather runtime statistics currently not activated. This needs some more work.
		#for file in $(ls results/statistics/benchmarks/busco/run_busco_*); do printf $file"\t"; sed '2q;d' results/statistics/benchmarks/busco/$file; done > results/statistics/benchmark_all_busco_runs.bench	

		# get statistics files for report. This is still a bit hacky and could be solved better, especially for the genomes file:
		cat results/statistics/orthology_filtering/orthology_filtering_gene_* > results/statistics/orthology_filtering_gene_statistics.{params.current_hash}.txt
		files=$(ls results/statistics/orthology_filtering/orthology_filtering_genomes_*)
		cat $(echo $files | awk '{{print $1}}') > results/statistics/orthology_filtering_genomes_statistics.{params.current_hash}.txt		
		
		echo "$(date) - Number of BUSCO sequence files: $(ls results/orthology/busco/busco_sequences.{params.current_hash}/*/*.fas | wc -l)" >> results/statistics/runlog.txt
		echo "$(date) - Number of deduplicated BUSCO sequence files: $(ls results/orthology/busco/busco_sequences_deduplicated.{params.current_hash}/*.fas | wc -l)" >> results/statistics/runlog.txt
		touch {output.checkpoint}
		"""

rule filter_orthology:
	input:
		"results/checkpoints/modes/orthology."+previous_hash+".done",
		"results/checkpoints/remove_duplicated_sequence_files."+current_hash+".done"
	output:
		"results/checkpoints/modes/filter_orthology."+current_hash+".done"
	shell:
		"""
		echo "$(date) - Phylociraptor filter-orthology done." >> results/statistics/runlog.txt
		touch {output}
		"""
