configfile:  "data/config.yaml"
import subprocess


BUSCOS, = glob_wildcards("results/busco_sequences_deduplicated/{busco}_all.fas")
if config["alignment"]["method"] == "clustalo":
	rule align:
		input:
			checkpoint = "results/checkpoints/remove_duplicated_sequence_files.done",
			sequence_file = "results/busco_sequences_deduplicated/{busco}_all.fas"
		output:
			alignment = "results/alignments/{busco}_aligned.fas",
		benchmark:
			"results/statistics/benchmarks/align/align_{busco}.txt"
		singularity:
			"docker://reslp/clustalo:1.2.4"
		threads:
			int(config["alignment"]["threads"])
		params:
			config["alignment"]["parameters"]
		shell:
			"""
			clustalo -i {input.sequence_file} --threads={threads} {params} > {output.alignment} 
			"""

else:
	rule align:
		input:
			checkpoint = "results/checkpoints/remove_duplicated_sequence_files.done",
			sequence_file = "results/busco_sequences_deduplicated/{busco}_all.fas"
		output:
			alignment = "results/alignments/{busco}_aligned.fas",
		benchmark:
			"results/statistics/benchmarks/align/align_{busco}.txt"
		singularity:
			"docker://reslp/mafft:7.464"
		threads:
			int(config["alignment"]["threads"])
		params:
			config["alignment"]["parameters"]
		shell:
			"""
			mafft {params} {input.sequence_file} > {output.alignment}
			"""


rule aggregate_align:
	input:
		expand("results/alignments/{bus}_aligned.fas", bus=BUSCOS)
	output:
		checkpoint = "results/checkpoints/aggregate_align.done"
	shell:
		"""
		touch {output.checkpoint}
		"""
rule get_alignment_statistics:
	input:
		rules.aggregate_align.output.checkpoint
	output:
		statistics_alignment = "results/statistics/statistics_alignments.txt",
		overview_statistics = "results/statistics/align_trim_overview_statistics.txt"
	params:
		ids = config["species"],
		datatype = config["filtering"]["seq_type"],
		alignment_method = config["alignment"]["method"],
		alignment_params = config["alignment"]["parameters"],
		trimming_method = config["trimming"]["method"],
		trimming_params = config["trimming"]["parameters"],
		pars_sites = config["filtering"]["min_parsimony_sites"]
	singularity: "docker://reslp/concat:0.21"
	shell:
		"""
		tail -n +2 {params.ids} | awk -F "," '{{print $1;}}' | sed 's/^"//g' | sed 's/"$//g' | sed 's/ /_/g' > results/statistics/ids_alignments.txt
		# here the ids for the alignments need to be filtered as well first. maybe this can be changed in the concat.py script, so that an id file is not needed anymore.
		concat.py -d results/alignments/ -t results/statistics/ids_alignments.txt --runmode concat -o results/statistics/ --biopython --statistics --seqtype aa --noseq
		mv results/statistics/statistics.txt {output.statistics_alignment}
		echo "Alignment method:\t{params.alignment_method}" > {output.overview_statistics}
		echo "Alignment parameters:\t{params.alignment_params}" >> {output.overview_statistics}
		echo "Trimming method:\t{params.trimming_method}" >> {output.overview_statistics}
		echo "Trimming parameters:\t{params.trimming_params}" >> {output.overview_statistics}
		echo "Min Parssites:\t{params.pars_sites}" >> {output.overview_statistics}
		"""
rule all_align:
	input:
		"results/checkpoints/aggregate_align.done",
		"results/statistics/statistics_alignments.txt"
	output:
		"checkpoints/align_trim.done"
	shell:
		"""
		touch {output}
		echo "$(date) - Pipeline part 2 (align) done." >> results/statistics/runlog.txt
		"""
