rule mafft:
		input:
			"results/alignments/full/mafft.{hash}/parameters.align.mafft.{hash}.yaml",
			sequence_file = "results/orthology/single-copy-orthologs_deduplicated." + hashes["filter-orthology"]["global"] + "/{busco}_all.fas",
		output:
			alignment = "results/alignments/full/mafft.{hash}/{busco}_aligned.fas",
		benchmark:
			"results/statistics/benchmarks/align/mafft_align_{busco}.{hash}.txt"
		log:
			"log/align/mafft/mafft_align_{busco}.{hash}.log.txt"
		singularity:
			containers["mafft"]
		threads:
			int(config["alignment"]["threads"])
		params:
			config["alignment"]["options"]["mafft"]
		shell:
			"""
			mafft --thread {threads} $(if [[ "{params}" != "None" ]]; then echo {params}; fi) {input.sequence_file} 1> {output.alignment} 2> {log}
			"""
