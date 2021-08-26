configfile:  "data/config.yaml"
import subprocess


BUSCOS, = glob_wildcards("results/orthology/busco/busco_sequences_deduplicated/{busco}_all.fas")

rule align_clustalo:
		input:
			checkpoint = "checkpoints/filter_orthology.done",
			sequence_file = "results/orthology/busco/busco_sequences_deduplicated/{busco}_all.fas"
		output:
			alignment = "results/alignments/full/clustalo/{busco}_aligned.fas",
		benchmark:
			"results/statistics/benchmarks/align/clustalo_align_{busco}.txt"
		singularity:
			"docker://reslp/clustalo:1.2.4"
		threads:
			int(config["alignment"]["threads"])
		params:
			config["alignment"]["clustalo_parameters"]
		shell:
			"""
			clustalo -i {input.sequence_file} --threads={threads} $(if [[ "{params}" != "None" ]]; then echo {params}; fi) > {output.alignment} 
			"""

rule align_mafft:
		input:
			checkpoint = "checkpoints/filter_orthology.done",
			sequence_file = "results/orthology/busco/busco_sequences_deduplicated/{busco}_all.fas"
		output:
			alignment = "results/alignments/full/mafft/{busco}_aligned.fas",
		benchmark:
			"results/statistics/benchmarks/align/mafft_align_{busco}.txt"
		singularity:
			"docker://reslp/mafft:7.464"
		threads:
			int(config["alignment"]["threads"])
		params:
			config["alignment"]["mafft_parameters"]
		shell:
			"""
			mafft --thread {threads} $(if [[ "{params}" != "None" ]]; then echo {params}; fi) {input.sequence_file} > {output.alignment}
			"""

rule aggregate_align:
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
		rules.aggregate_align.output.checkpoint
	output:
		statistics_alignment = "results/statistics/{aligner}_statistics_alignments-{batch}-"+str(config["filtering"]["concurrency"])+".txt",
		overview_statistics = "results/statistics/{aligner}_align_trim_overview_statistics-{batch}-"+str(config["filtering"]["concurrency"])+".txt"
	params:
		wd = os.getcwd(),
		ids = config["species"],
		datatype = config["filtering"]["seq_type"],
		alignment_method = "{aligner}",
		mafft_alignment_params = config["alignment"]["mafft_parameters"],
		clustalo_alignment_params = config["alignment"]["clustalo_parameters"],
		pars_sites = config["filtering"]["min_parsimony_sites"],
		nbatches = config["filtering"]["concurrency"],
	singularity: "docker://reslp/concat:0.21"
	shadow: "minimal"
	shell:
		"""
		# here the ids for the alignments need to be filtered as well first. maybe this can be changed in the concat.py script, so that an id file is not needed anymore.
		concat.py -i $(ls -1 {params.wd}/results/alignments/full/{wildcards.aligner}/* | sed -n '{wildcards.batch}~{params.nbatches}p' | tr '\\n' ' ') -t <(ls -1 {params.wd}/results/orthology/busco/busco_runs/) --runmode concat -o results/statistics/ --biopython --statistics --seqtype {params.datatype} --noseq
		mv results/statistics/statistics.txt {output.statistics_alignment}
		echo "Alignment method:\t{params.alignment_method}" > {output.overview_statistics}
		if [[ "{wildcards.aligner}" == "mafft" ]]; then echo "Alignment parameters:\t{params.mafft_alignment_params}" >> {output.overview_statistics}; fi
		if [[ "{wildcards.aligner}" == "clustalo" ]]; then echo "Alignment parameters:\t{params.clustalo_alignment_params}" >> {output.overview_statistics}; fi
		echo "Min Parssites:\t{params.pars_sites}" >> {output.overview_statistics}
		"""
rule all_align:
	input:
#		expand("results/alignments/full/{aligner}/{busco}_aligned.fas", aligner=config["alignment"]["method"], busco=BUSCOS),
#		expand("results/checkpoints/{aligner}_aggregate_align.done", aligner=config["alignment"]["method"]),
		expand("results/statistics/{aligner}_statistics_alignments-{batch}-"+str(config["filtering"]["concurrency"])+".txt", aligner=config["alignment"]["method"], batch=range(1 , config["filtering"]["concurrency"]+1)),
	output:
		"checkpoints/align.done"
	shell:
		"""
		touch {output}
		echo "$(date) - Pipeline part 2 (align) done." >> results/statistics/runlog.txt
		"""
