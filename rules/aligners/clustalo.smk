rule clustalo:
		input:
			"results/alignments/full/clustalo.{hash}/parameters.align.clustalo.{hash}.yaml",
			sequence_file = "results/orthology/single-copy-orthologs." + hashes["filter-orthology"]["global"] + "/{busco}_all.fas",
		output:
			alignment = "results/alignments/full/clustalo.{hash}/{busco}_aligned.fas",
		benchmark:
			"results/statistics/benchmarks/align/clustalo_align_{busco}.{hash}.txt"
		log:
			"log/align/clustalo/clustalo_align_{busco}.{hash}.log.txt"
		singularity:
			containers["clustalo"]
		threads:
			int(config["alignment"]["threads"])
		params:
			config["alignment"]["options"]["clustalo"]
		shell:
			"""
			clustalo -i {input.sequence_file} --threads={threads} $(if [[ "{params}" != "None" ]]; then echo {params}; fi) 1> {output.alignment} 2> {log}
			"""
