configfile: "data/config.yaml"

include: "concatenate.smk"

rule njtree:
	input:
		rules.concatenate.output.stockholm_alignment
	output:
		"results/phylogeny-{bootstrap}/njtree/{aligner}-{alitrim}/njtree.tre"
	singularity: "docker://reslp/quicktree:2.5"
	shell:
		"""
		quicktree -in a {input} > {output}
		"""
		
rule all_njtree:
	input:
		expand("results/phylogeny-{bootstrap}/njtree/{aligner}-{alitrim}/njtree.tre", aligner=config["alignment"]["method"], alitrim=config["trimming"]["method"], bootstrap=config["filtering"]["bootstrap_cutoff"])
	output:
		"results/checkpoints/modes/njtree.done"
	shell:
		"""
		echo "$(date) - Pipeline part ntree (njtree) done." >> results/statistics/runlog.txt
		touch {output}
		"""
