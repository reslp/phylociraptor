import yaml

# get list of containers to use:
with open("data/containers.yaml", "r") as yaml_stream:
    containers = yaml.safe_load(yaml_stream)

include: "concatenate.smk"

rule njtree:
	input:
		rules.concatenate.output.stockholm_alignment
	output:
		"results/phylogeny-{bootstrap}/njtree/{aligner}-{alitrim}/njtree.tre"
	singularity: containers["quicktree"] 
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
