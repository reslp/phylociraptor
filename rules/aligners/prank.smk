rule prank:
		input:
			"results/alignments/full/prank.{hash}/parameters.align.prank.{hash}.yaml",
			sequence_file = "results/orthology/single-copy-orthologs." + hashes["filter-orthology"]["global"] + "/{busco}_all.fas",
		output:
			alignment = "results/alignments/full/prank.{hash}/{busco}_aligned.fas",
		benchmark:
			"results/statistics/benchmarks/align/prank_align_{busco}.{hash}.txt"
		log:
			"log/align/prank/prank_align_{busco}.{hash}.log.txt"
		singularity:
			containers["prank"]
		threads:
			int(config["alignment"]["threads"])
		params:
			config["alignment"]["options"]["prank"]
		shadow:
			"shallow"
		shell:
			"""
			prank -d={input.sequence_file} -o={output.alignment} $(if [[ "{params}" != "None" ]]; then echo {params}; fi) 

			# This next part is necessary because prank allways appends .best.fas to the specified output file.
			echo "Moving file: {output.alignment}.best.fas => {output.alignment}"
			mv {output.alignment}.best.fas {output.alignment}
			echo "done"

			"""
