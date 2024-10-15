def return_untrimmed_params(wildcards):
	return "results/alignments/trimmed/"+wildcards.aligner+"-untrimmed."+trimmed_hashes["untrimmed"][wildcards.aligner]+"/parameters.filter-align."+wildcards.aligner+"-untrimmed."+trimmed_hashes["untrimmed"][wildcards.aligner]+".yaml"

rule untrimmed:
		input:
			checkpoint = return_aligner_checkpoint,
			params = return_untrimmed_params,
			alignment = per_aligner
		output:
			trimmed_alignment = "results/alignments/trimmed/{aligner}-untrimmed.{hash}/{busco}_aligned_trimmed.fas",
		benchmark:
			"results/statistics/benchmarks/align/{aligner}_untrimmed_{busco}.{hash}.txt"
		log:
			"log/filter-align/untrimmed/untrimmed_{aligner}_{busco}.{hash}.txt"
		shell:
			"""
			echo "Aligned remains untrimmed: {input.alignment}" > {log}
			cp {input.alignment} {output.trimmed_alignment}
			"""
