include: "functions.smk"
import yaml

# get list of containers to use:
with open("data/containers.yaml", "r") as yaml_stream:
    containers = yaml.safe_load(yaml_stream)

include: "concatenate.smk"

aligners = get_aligners()		
trimmers = get_trimmers()		
bscuts = get_bootstrap_cutoffs()

rule quicktree:
	input:
		rules.concatenate.output.stockholm_alignment
	output:
		"results/phylogeny-{bootstrap}/njtree/{aligner}-{alitrim}/njtree.tre"
	singularity: containers["quicktree"] 
	shell:
		"""
		quicktree -in a {input} > {output}
		"""
		
rule njtree:
	input:
		expand("results/phylogeny-{bootstrap}/njtree/{aligner}-{alitrim}/njtree.tre", aligner=aligners, alitrim=trimmers, bootstrap=bscuts)
	output:
		"results/checkpoints/modes/njtree.done"
	shell:
		"""
		echo "$(date) - Pipeline part ntree (njtree) done." >> results/statistics/runlog.txt
		touch {output}
		"""
