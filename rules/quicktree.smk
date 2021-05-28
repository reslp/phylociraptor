configfile: "data/config.yaml"

include: "concatenate.smk"

rule njtree:
	input:
		rules.concatenate.output.stockholm_alignment
	output:
		"results/phylogeny/njtree/njtree.tre"
	singularity: "docker://reslp/quicktree:2.5"
	shell:
		"""
		quicktree -in a {input} > {output}
		echo "$(date) - Pipeline part ntree (njtree) done." >> results/statistics/runlog.txt
		"""
		
