rule iqtree_gene_trees:
	input:
		rules.part2.output,
		busco = "results/filtered_alignments/{busco}_aligned_trimmed.fas"
	output:
		checkpoint = "results/checkpoints/gene_trees/{busco}_genetree.done",
		trees = "results/phylogeny/gene_trees/{busco}/{busco}.treefile"
	params:
		wd = os.getcwd(),
		maxmem = config["iqtree"]["maxmem"],
		busco = "{busco}"
	threads:
		config["iqtree"]["threads"]
	singularity:
		"docker://reslp/iqtree:2.0.7"
	shell:
		"""
		cd results/phylogeny/gene_trees/{params.busco}
		cp {params.wd}/results/filtered_alignments/{params.busco}_aligned_trimmed.fas {params.busco}_aligned_trimmed.fas
		
		# here we decide how iqtree should be run. In case modeltesting was run before, this will not be repeated here.
		if [[ -f {params.wd}/results/modeltest/best_models.txt && {params.wd}/checkpoints/part_model.done ]]; then
			model=$(cat {params.wd}/results/modeltest/best_models.txt | grep {params.busco} | awk '{{print $2}}')
			echo "$(date) - phylociraptor was run with modeltesting before. Will run iqtree gene tree for {params.busco} with best model: $model" >> {params.wd}/results/report.txt
			if [[ -z "{params.maxmem}" ]]; then
				iqtree -s {params.busco}_aligned_trimmed.fas -m $model --prefix {params.busco}_gt -nt AUTO -ntmax {threads} -redo
			else
				iqtree -s {params.busco}_aligned_trimmed.fas -m $model --prefix {params.busco}_gt -nt AUTO -ntmax {threads} -redo -mem {params.maxmem}
			fi
		else
			echo "$(date) - phylociraptor will run iqtree gene tree for {params.busco}  now, with automated model testing." >> {params.wd}/results/report.txt
			if [[ -z "{params.maxmem}" ]]; then
				iqtree -s {params.busco}_aligned_trimmed.fas --prefix {params.busco}_gt -nt AUTO -ntmax {threads} -m MFP -redo
			else
				iqtree {params.busco}_aligned_trimmed.fas --prefix {params.busco}_gt -nt AUTO -ntmax {threads} -m MFP -redo -mem {params.maxmem}
			fi
		fi
		touch {params.wd}/{output.checkpoint}
		"""
BUSCOS, = glob_wildcards("results/filtered_alignments/{busco}_aligned_trimmed.fas")

rule aggregate_gene_trees:
	input:
		treefiles = expand("results/phylogeny/gene_trees/{busco}/{busco}.treefile", busco=BUSCOS),
		checkpoint = expand("results/checkpoints/gene_trees/{busco}_genetree.done", busco=BUSCOS)	
	output:
		trees = "results/gene_trees/all_gene_trees.tre",
		checkpoint = "results/checkpoints/aggregate_gene_trees.done"
	shell:
		"""
		cat {input.treefiles} > {output.trees}
		touch {output.checkpoint}
		"""

rule astral_species_tree:
	input:
		trees = rules.aggregate_gene_trees.output.trees,
		checkpoint = rules.aggregate_gene_trees.output.checkpoint
	output:
		species_tree = "results/phylogeny/astral/species_tree.tre",
		checkpoint = "results/checkpoints/astral_species_tree.done"
	params:
		wd = os.getcwd()
	singularity:
		"docker://reslp/astral:5.7.1"
	shell:
		"""
		java -jar /ASTRAL-5.7.1/Astral/astral.5.7.1.jar -i {input.trees} -o {output.species_tree}
		touch {output.checkpoint}
		"""
