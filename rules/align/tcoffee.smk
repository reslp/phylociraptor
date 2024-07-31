rule tcoffee:
		input:
			"results/alignments/full/tcoffee.{hash}/parameters.align.tcoffee.{hash}.yaml",
			sequence_file = "results/orthology/single-copy-orthologs." + hashes["filter-orthology"]["global"] + "/{busco}_all.fas",
		output:
			alignment = "results/alignments/full/tcoffee.{hash}/{busco}_aligned.fas",
		benchmark:
			"results/statistics/benchmarks/align/tcoffee_align_{busco}.{hash}.txt"
		log:
			"log/align/tcoffee/tcoffee_align_{busco}.{hash}.log.txt"
		singularity:
			containers["tcoffee"]
		threads:
			int(config["alignment"]["threads"])
		params:
			config["alignment"]["options"]["tcoffee"]
		shadow:
			"shallow"
		shell:
			"""
			t_coffee {input.sequence_file} -quiet -thread {threads} $(if [[ "{params}" != "None" ]]; then echo {params}; fi) -output=fasta_aln -outfile=stdout -out_lib=stderr 1> {output.alignment} 2> /dev/null

			"""
