rule MYALIGNER:
		input:
			"results/alignments/full/MYALIGNER.{hash}/parameters.align.MYALIGNER.{hash}.yaml",
			sequence_file = "results/orthology/single-copy-orthologs_deduplicated." + hashes["filter-orthology"]["global"] + "/{busco}_all.fas",
		output:
			alignment = "results/alignments/full/MYALIGNER.{hash}/{busco}_aligned.fas",
		benchmark:
			"results/statistics/benchmarks/align/MYALIGNER_align_{busco}.{hash}.txt"
		log:
			"log/align/MYALIGNER/MYALIGNER_align_{busco}.{hash}.log.txt"
		singularity:
			"docker://repository/MYALIGNER:TAG"
		threads:
			int(config["alignment"]["threads"])
		params:
			config["alignment"]["options"]["MYALIGNER"]
		shell:
			"""
			MYALIGNER -i {input.sequence_file} -t {threads} $(if [[ "{params}" != "None" ]]; then echo {params}; fi) -o {output.alignment}
			"""
