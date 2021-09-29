configfile:  "data/config.yaml"

import os.path
import glob

BUSCOS, = glob_wildcards("results/orthology/busco/busco_sequences_deduplicated/{busco}_all.fas")

def return_trees(wildcards):
	lis = []
        for busco in BUSCOS:
		if os.path.isfile("results/alignments/filtered/"+wildcards.aligner+"-"+wildcards.alitrim+"/"+str(busco)+"_aligned_trimmed.fas"):
			lis.append("results/phylogeny/gene_trees/"+wildcards.aligner+"-"+wildcards.alitrim+"/"+busco+"/"+busco+"_gt.treefile")
	return lis

rule iqtree_gene_trees:
	input:
		"results/statistics/filter-{aligner}-{alitrim}/alignment_filter_information_{alitrim}_{aligner}.txt"
	output:
		checkpoint = "results/checkpoints/gene_trees/{aligner}-{alitrim}/{busco}_genetree.done",
		trees = "results/phylogeny/gene_trees/{aligner}-{alitrim}/{busco}/{busco}_gt.treefile"
	benchmark:
		"results/statistics/benchmarks/speciestree/iqtree_gene_tree_{busco}_{aligner}_{alitrim}.txt"
	params:
		wd = os.getcwd(),
		models = "results/modeltest/best_models_{aligner}_{alitrim}.txt",
		maxmem = config["iqtree"]["maxmem"],
		busco = "{busco}",
		bb = config["genetree"]["boostrap"],
		additional_params = config["iqtree"]["additional_params"]
	threads:
		config["genetree"]["threads"]
	singularity:
		"docker://reslp/iqtree:2.0.7"
	shell:
		"""
		cd results/phylogeny/gene_trees/{wildcards.aligner}-{wildcards.alitrim}/{params.busco}
		cp {params.wd}/results/alignments/filtered/{wildcards.aligner}-{wildcards.alitrim}/{params.busco}_aligned_trimmed.fas {params.busco}_aligned_trimmed.fas
		
		# here we decide how iqtree should be run. In case modeltesting was run before, this will not be repeated here.
		if [[ -f {params.wd}/{params.models} && {params.wd}/checkpoints/modeltest.done ]]; then
			model=$(cat {params.wd}/{params.models} | grep -P "^{params.busco}\\t" | awk '{{print $2}}')
			echo "$(date) - phylociraptor was run with modeltesting before. Will run iqtree gene tree for {params.busco} with best model: $model" >> {params.wd}/results/statistics/runlog.txt
			if [[ -z "{params.maxmem}" ]]; then
				iqtree -s {params.busco}_aligned_trimmed.fas -m $model --prefix {params.busco}_gt -nt {threads} -redo $(if [[ "{params.bb}" != "None" ]]; then echo "-bb {params.bb}"; fi) {params.additional_params}
			else
				iqtree -s {params.busco}_aligned_trimmed.fas -m $model --prefix {params.busco}_gt -nt {threads} -redo -mem {params.maxmem} $(if [[ "{params.bb}" != "None" ]]; then echo "-bb {params.bb}"; fi) {params.additional_params}
			fi
		else
			echo "$(date) - phylociraptor will run iqtree gene tree for {params.busco}  now, with automated model testing." >> {params.wd}/results/statistics/runlog.txt
			if [[ -z "{params.maxmem}" ]]; then
				iqtree -s {params.busco}_aligned_trimmed.fas --prefix {params.busco}_gt -nt {threads} -m MFP -redo $(if [[ "{params.bb}" != "None" ]]; then echo "-bb {params.bb}"; fi) {params.additional_params}
			else
				iqtree {params.busco}_aligned_trimmed.fas --prefix {params.busco}_gt -nt {threads} -m MFP -redo -mem {params.maxmem} $(if [[ "{params.bb}" != "None" ]]; then echo "-bb {params.bb}"; fi) {params.additional_params}
			fi
		fi
		rm {params.busco}_aligned_trimmed.fas
		touch {params.wd}/{output.checkpoint}
		"""

rule aggregate_gene_trees:
	input:
		treefiles = return_trees
	output:
		trees = "results/phylogeny/astral/{aligner}-{alitrim}/trees_{aligner}_{alitrim}.tre",
		checkpoint = "results/checkpoints/aggregate_gene_trees_{aligner}_{alitrim}.done"
	params:
		include = config["speciestree"]["include"]
	shell:
		"""
		if [[ "{params.include}" == "None" ]]
		then
			cat $(ls -1 results/phylogeny/gene_trees/{wildcards.aligner}-{wildcards.alitrim}/*/*_gt.treefile) > {output.trees}
		else
			cat $(cat {params.include}) > {output.trees}
		fi
		touch {output.checkpoint}
		"""
#		for treefile in $treefiles
#		do
#			cat $treefile
#		done > {output.trees}
#		cat {input.treefiles} > {output.trees}


rule astral_species_tree:
	input:
		trees = "results/phylogeny/astral/{aligner}-{alitrim}/trees_{aligner}_{alitrim}.tre" 
#		trees = rules.aggregate_gene_trees.output.trees,
#		checkpoint = rules.aggregate_gene_trees.output.checkpoint
	output:
		species_tree = "results/phylogeny/astral/{aligner}-{alitrim}/species_tree.tre",
		checkpoint = "results/checkpoints/astral_species_tree_{aligner}_{alitrim}.done"
	benchmark:
		"results/statistics/benchmarks/speciestree/astral_species_tree_{aligner}_{alitrim}.txt"
	params:
		wd = os.getcwd()
	singularity:
		"docker://reslp/astral:5.7.1"
	shell:
		"""
		java -jar /ASTRAL-5.7.1/Astral/astral.5.7.1.jar -i {input.trees} -o {output.species_tree}
		touch {output.checkpoint}
		"""
rule all_speciestree:
	input:
		expand("results/checkpoints/astral_species_tree_{aligner}_{alitrim}.done", aligner=config["alignment"]["method"], alitrim=config["trimming"]["method"])
	output:
		"results/checkpoints/modes/speciestree.done"
	shell:
		"""
		echo "$(date) - Speciestree reconstruction done." >> results/statistics/runlog.txt
		touch {output}
		"""

