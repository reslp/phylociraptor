import subprocess

BUSCOS, = glob_wildcards("results/busco_sequences_deduplicated/{busco}_all.fas")

rule filter_alignments:
	input:
		alignments = expand("results/trimmed_alignments/{bus}_aligned_trimmed.fas", bus=BUSCOS),
		stats = rules.get_alignment_statistics.output.statistics_trimmed
	output:
		checkpoint = "results/checkpoints/filter_alignments.done",
		filter_info = "results/statistics/alignment_filter_information.txt"
	benchmark:
		"results/statistics/benchmarks/align/filter_alignments.txt"
	singularity:
		"docker://reslp/biopython_plus:1.77"
	params:
		wd = os.getcwd(),
		trimming_method = config["trimming"]["method"],
		min_pars_sites = config["filtering"]["min_parsimony_sites"]
	shell:
		"""
		if [[ -d results/filtered_alignments ]]; then
			rm -rf results/filtered_alignments
		fi
		mkdir results/filtered_alignments
		
		for file in results/trimmed_alignments/*.fas;
		do
			python bin/filter_alignments.py --alignments {params.wd}/$file --outdir "{params.wd}/results/filtered_alignments" --statistics-file {input.stats} --min-parsimony {params.min_pars_sites} >> {output.filter_info}
		done
		echo "$(date) - Number of alignments: $(ls results/alignments/*.fas | wc -l)" >> results/statistics/runlog.txt
		echo "$(date) - Number of trimmed alignments: $(ls results/filtered_alignments/*.fas | wc -l)" >> results/statistics/runlog.txt
		echo "$(date) - Number of alignments after filtering: $(ls results/filtered_alignments/*.fas | wc -l)" >> results/statistics/runlog.txt
		touch {output.checkpoint}
		"""

rule get_filter_statistics:
	input:
		rules.filter_alignments.output.checkpoint
	output:
		statistics_filtered = "results/statistics/statistics_filtered.txt"	
	params:
		ids = config["species"],
		datatype = config["filtering"]["seq_type"],
	singularity: "docker://reslp/concat:0.21"
	shell:
		"""
		tail -n +2 {params.ids} | awk -F "," '{{print $1;}}' | sed 's/^"//g' | sed 's/"$//g' | sed 's/ /_/g' > results/statistics/ids_alignments.txt
		# here the ids for the alignments need to be filtered as well first. maybe this can be changed in the concat.py script, so that an id file is not needed anymore.
		concat.py -d results/filtered_alignments/ -t results/statistics/ids_alignments.txt --runmode concat -o results/statistics/ --biopython --statistics --seqtype {params.datatype} --noseq
		mv results/statistics/statistics.txt {output.statistics_filtered}
		"""
