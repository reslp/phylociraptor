rule bmge:
		input:
			checkpoint = return_aligner_checkpoint,
			params = return_bmge_params,
			alignment = per_aligner
		output:
			trimmed_alignment = "results/alignments/trimmed/{aligner}-bmge.{hash}/{busco}_aligned_trimmed.fas",
		benchmark:
			"results/statistics/benchmarks/align/{aligner}_bmge_{busco}.{hash}.txt"
		params:
			trimmer = config["trimming"]["options"]["bmge"],
			type = config["filtering"]["seq_type"] 
		log:
			"log/filter-align/bmge/bmge_{aligner}_{busco}.{hash}.txt"
		singularity:
			containers["bmge"]	
		shell:
			"""
			if [[ "{params.type}" == "aa" ]]
			then
				echo "Setting sequence type to AA" 2>&1 | tee {log}
				seqtype="AA"
			elif [[ "{params.type}" == "nu" ]]
			then
				seqtype="DNA"
				echo "Setting sequence type to DNA" 2>&1 | tee {log}
			fi
			bmge -i {input.alignment} -t $seqtype -of {output.trimmed_alignment} 2>&1 | tee -a {log}
			"""
