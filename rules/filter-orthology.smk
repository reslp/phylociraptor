import yaml

# get list of containers to use:
with open("data/containers.yaml", "r") as yaml_stream:
    containers = yaml.safe_load(yaml_stream)

def determine_concurrency_limit():
	fname = "results/orthology/busco/busco_table.txt"
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

rule create_sequence_files:
	input:
		table = "results/orthology/busco/busco_table.txt"
	output:
		sequence_dir=directory("results/orthology/busco/busco_sequences/{batch}-"+str(config["concurrency"])),
		checkpoint = "results/checkpoints/create_sequence_files_{batch}-"+str(config["concurrency"])+".done",
		genome_statistics = "results/statistics/orthology_filtering/orthology_filtering_genomes_statistics_{batch}-"+str(config["concurrency"])+".txt",
		gene_statistics = "results/statistics/orthology_filtering/orthology_filtering_gene_statistics_{batch}-"+str(config["concurrency"])+".txt"
	benchmark:
		"results/statistics/benchmarks/create_seq_files/create_sequence_files_{batch}-"+str(config["concurrency"])+".txt"
	params:
		cutoff=config["filtering"]["cutoff"],
		minsp=config["filtering"]["minsp"],
		busco_dir = "results/orthology/busco/busco_runs",
		seqtype = config["filtering"]["seq_type"],
		nbatches = config["concurrency"],
	singularity:
		containers["biopython"]
	log: "log/exSeqfiles_{batch}-"+str(config["concurrency"])+".log"
	shell:
		""" 
		if [[ ! -d {output.sequence_dir} ]]; then mkdir -p {output.sequence_dir}; fi
		# remove files in case there are already some:
		rm -f {output.sequence_dir}/*
		python bin/create_sequence_files.py --type {params.seqtype} --busco_table {input.table} --busco_results {params.busco_dir} --cutoff {params.cutoff} --outdir {output.sequence_dir} --minsp {params.minsp} --genome_statistics {output.genome_statistics} --gene_statistics {output.gene_statistics} --batchID {wildcards.batch} --nbatches {params.nbatches} &> {log}
		touch {output.checkpoint}   
		"""

rule remove_duplicated_sequence_files:
	input:
		checkpoint = expand(rules.create_sequence_files.output.checkpoint, batch=batches)
	output:
		checkpoint = "results/checkpoints/remove_duplicated_sequence_files.done",
		statistics = "results/statistics/duplicated_sequences_handling_information.txt"
	benchmark:
		"results/statistics/benchmarks/remdupseq/remove_duplicated_sequence_files.txt"
	singularity: containers["biopython"] 
	params:
		wd = os.getcwd(),
		minsp=config["filtering"]["minsp"],
		dupseq=config["filtering"]["dupseq"]
	shell:
		"""
		if [[ -d results/orthology/busco/busco_sequences_deduplicated ]]; then
			rm -rf results/orthology/busco/busco_sequences_deduplicated
		fi
		mkdir results/orthology/busco/busco_sequences_deduplicated
		
		rm -f {output.statistics}
		
		if [[ {params.dupseq} == "persample" ]];
		then
			echo "$(date) - BUSCO files will be filtered on a per-sample basis. This could lower the number of species in the final tree." >> results/statistics/runlog.txt
		else
			echo "$(date) - BUSCO files will be filtered on a per BUSCO gene basis. This could lower the number of genes used to calculate the final tree." >> results/statistics/runlog.txt
		fi
		
		for file in results/orthology/busco/busco_sequences/*/*.fas;
			do
			if [[ {params.dupseq} == "persample" ]];
			then
				# per sequence filtering
				python bin/filter_alignments.py --alignments {params.wd}/$file --outdir "{params.wd}/results/orthology/busco/busco_sequences_deduplicated" --per_sample --minsp {params.minsp} >> {output.statistics}
			else
				# whole alignment filtering
				python bin/filter_alignments.py --alignments {params.wd}/$file --outdir "{params.wd}/results/orthology/busco/busco_sequences_deduplicated" >> {output.statistics}
			fi
		done
		#gather runtime statistics currently not activated. This needs some more work.
		#for file in $(ls results/statistics/benchmarks/busco/run_busco_*); do printf $file"\t"; sed '2q;d' results/statistics/benchmarks/busco/$file; done > results/statistics/benchmark_all_busco_runs.bench	

		# get statistics files for report. This is still a bit hacky and could be solved better, especially for the genomes file:
		cat results/statistics/orthology_filtering/orthology_filtering_gene_* > results/statistics/orthology_filtering_gene_statistics.txt
		files=$(ls results/statistics/orthology_filtering/orthology_filtering_genomes_*)
		cat $(echo $files | awk '{{print $1}}') > results/statistics/orthology_filtering_genomes_statistics.txt		
		
		echo "$(date) - Number of BUSCO sequence files: $(ls results/orthology/busco/busco_sequences/*/*.fas | wc -l)" >> results/statistics/runlog.txt
		echo "$(date) - Number of deduplicated BUSCO sequence files: $(ls results/orthology/busco/busco_sequences_deduplicated/*.fas | wc -l)" >> results/statistics/runlog.txt
		touch {output.checkpoint}
		"""

rule filter_orthology:
	input:
		"results/checkpoints/remove_duplicated_sequence_files.done"
	output:
		"results/checkpoints/modes/filter_orthology.done"
	shell:
		"""
		echo "$(date) - Pipeline part filter-orthology done." >> results/statistics/runlog.txt
		touch {output}
		"""
