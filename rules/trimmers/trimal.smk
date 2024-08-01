rule trimal:
		input:
			checkpoint = return_aligner_checkpoint,
			params = return_trimal_params,
			alignment = per_aligner
		output:
			trimmed_alignment = "results/alignments/trimmed/{aligner}-trimal.{hash}/{busco}_aligned_trimmed.fas",
		benchmark:
			"results/statistics/benchmarks/align/{aligner}_trimal_{busco}.{hash}.txt"
		params:
			trimmer = config["trimming"]["options"]["trimal"]
		log:
			"log/filter-align/trimal/trimal_{aligner}_{busco}.{hash}.txt"
		singularity:
			containers["trimal"]	
		shell:
			"""
			trimal {params.trimmer} -in {input.alignment} -out {output.trimmed_alignment} 2>&1 | tee {log}
			#it can happen that the alignment is empty after trimming. In this case trimal doesn't produce an error, but does not write the file so we do that manually
			if [[ ! -f "{output.trimmed_alignment}" ]]; then echo "Making dummy file" 2>&1 | tee {log}; touch {output.trimmed_alignment}; fi
			"""
