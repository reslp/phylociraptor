rule muscle:
		input:
			"results/alignments/full/muscle.{hash}/parameters.align.muscle.{hash}.yaml",
			sequence_file = "results/orthology/single-copy-orthologs." + hashes["filter-orthology"]["global"] + "/{busco}_all.fas",
		output:
			alignment = "results/alignments/full/muscle.{hash}/{busco}_aligned.fas",
		benchmark:
			"results/statistics/benchmarks/align/muscle_align_{busco}.{hash}.txt"
		log:
			"log/align/muscle/muscle_align_{busco}.{hash}.log.txt"
		singularity:
			containers["muscle"]
		threads:
			int(config["alignment"]["threads"])
		params:
			config["alignment"]["options"]["muscle"]
		shell:
			"""
			muscle -super5 {input.sequence_file} -threads {threads} $(if [[ "{params}" != "None" ]]; then echo {params}; fi) -output {output.alignment} 2> {log}
			"""
