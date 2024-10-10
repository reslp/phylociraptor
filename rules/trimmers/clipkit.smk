def return_clipkit_params(wildcards):
	return "results/alignments/trimmed/"+wildcards.aligner+"-clipkit."+trimmed_hashes["clipkit"][wildcards.aligner]+"/parameters.filter-align."+wildcards.aligner+"-clipkit."+trimmed_hashes["clipkit"][wildcards.aligner]+".yaml"


rule clipkit:
		input:
			checkpoint = return_aligner_checkpoint,
			params = return_clipkit_params,
			alignment = per_aligner
		output:
			trimmed_alignment = "results/alignments/trimmed/{aligner}-clipkit.{hash}/{busco}_aligned_trimmed.fas",
		benchmark:
			"results/statistics/benchmarks/align/{aligner}_clipkit_{busco}.{hash}.txt"
		params:
			trimmer = config["trimming"]["options"]["clipkit"]
		log:
			"log/filter-align/clipkit/clipkit_{aligner}_{busco}.{hash}.txt"
		singularity:
			containers["clipkit"]	
		shell:
			"""
			clipkit {input.alignment} {params.trimmer} -o {output.trimmed_alignment} 1> {log}
			"""
