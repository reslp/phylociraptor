rule create_sequence_files:
	input:
		busco_table = rules.extract_busco_table.output.busco_table,
		checkpoint = rules.orthology.output
	output:
		sequence_dir = directory("results/busco_sequences"),
		checkpoint = "results/checkpoints/create_sequence_files.done",
		genome_statistics = "results/statistics/orthology_filtering_genomes_statistics.txt",
		gene_statistics = "results/statistics/orthology_filtering_gene_statistics.txt"
	benchmark:
		"results/statistics/benchmarks/busco/create_sequence_files.txt"
	params:
		cutoff=config["filtering"]["cutoff"],
		minsp=config["filtering"]["minsp"],
		busco_dir = "results/busco",
		seqtype = config["filtering"]["seq_type"]
	singularity:
		"docker://reslp/biopython_plus:1.77"
	shell:
		""" 
		mkdir -p {output.sequence_dir}
		# remove files in case there are already some:
		rm -f {output.sequence_dir}/*
		python bin/create_sequence_files.py --type {params.seqtype} --busco_table {input.busco_table} --busco_results {params.busco_dir} --cutoff {params.cutoff} --outdir {output.sequence_dir} --minsp {params.minsp} --genome_statistics {output.genome_statistics} --gene_statistics {output.gene_statistics}
		touch {output.checkpoint}   
		"""

rule remove_duplicated_sequence_files:
	input:
		checkpoint = rules.create_sequence_files.output.checkpoint
	output:
		checkpoint = "results/checkpoints/remove_duplicated_sequence_files.done",
		statistics = "results/statistics/duplicated_sequences_handling_information.txt"
	benchmark:
		"results/statistics/benchmarks/busco/remove_duplicated_sequence_files.txt"
	singularity: "docker://reslp/biopython_plus:1.77"
	params:
		wd = os.getcwd(),
		dupseq=config["filtering"]["dupseq"]
	shell:
		"""
		if [[ -d results/busco_sequences_deduplicated ]]; then
			rm -rf results/busco_sequences_deduplicated
		fi
		mkdir results/busco_sequences_deduplicated
		
		rm -f {output.statistics}
		
		if [[ {params.dupseq} == "persample" ]];
		then
			echo "$(date) - BUSCO files will be filtered on a per-sample basis. This could lower the number of species in the final tree." >> results/statistics/runlog.txt
		else
			echo "$(date) - BUSCO files will be filtered on a per BUSCO gene basis. This could lower the number of genes used to calculate the final tree." >> results/statistics/runlog.txt
		fi
		
		for file in results/busco_sequences/*.fas;
			do
			if [[ {params.dupseq} == "persample" ]];
			then
				# per sequence filtering
				python bin/filter_alignments.py --alignments {params.wd}/$file --outdir "{params.wd}/results/busco_sequences_deduplicated" --per_sample >> {output.statistics}
			else
				# whole alignment filtering
				python bin/filter_alignments.py --alignments {params.wd}/$file --outdir "{params.wd}/results/busco_sequences_deduplicated" >> {output.statistics}
			fi
		done
		#gather runtime statistics
		for file in $(ls results/statistics/benchmarks/busco/); do printf $file"\t"; sed '2q;d' results/statistics/benchmarks/busco/$file; done > results/statistics/benchmark_all_busco_runs.bench	
		echo "$(date) - Number of BUSCO sequence files: $(ls results/busco_sequences/*.fas | wc -l)" >> results/statistics/runlog.txt
		echo "$(date) - Number of deduplicated BUSCO sequence files: $(ls results/busco_sequences_deduplicated/*.fas | wc -l)" >> results/statistics/runlog.txt
		touch {output.checkpoint}
		"""
