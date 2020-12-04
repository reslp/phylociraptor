rule iqtree_gene_trees:
	input:
		rules.part2.output
	output:
		checkpoint = "results/checkpoints/iqtree_gene_trees.done",
		trees = "results/phylogeny/gene_trees/loci.treefile"
	params:
		wd = os.getcwd(),
		maxmem = config["iqtree"]["maxmem"]
	threads:
		config["iqtree"]["threads"]
	singularity:
		"docker://reslp/iqtree:2.0.7"
	shell:
		"""
		rm -rf results/phylogeny/gene_trees/algn
		mkdir -p results/phylogeny/gene_trees/algn
		cd results/phylogeny/gene_trees
		cp {params.wd}/results/filtered_alignments/*.fas algn
		# here we decide how iqtree should be run. In case modetesting was run, this will not be repeated here.
		if [[ -f {params.wd}/results/modeltest/best_models.txt && {params.wd}/checkpoints/part_model.done ]]; then
			echo "$(date) - phylociraptor was run with -model before. Will run iqtree gene trees with best models." >> {params.wd}/results/report.txt
			echo "Will create NEXUS partition file for gene trees with model information now."
			echo "#nexus" > loci.nex
			echo "begin sets;" >> loci.nex
			cat {params.wd}/results/modeltest/best_models.txt | awk '{{print "charset part"NR" = algn/"$1"_aligned_trimmed.fas:*;"}}' >> loci.nex
			printf "charpartition mine = " >> loci.nex
			cat {params.wd}/results/modeltest/best_models.txt | awk '{{printf($2":part"NR", ")}}' | sed 's/\\(.*\\), /\\1;\\n/' >> loci.nex
			echo "end;" >> loci.nex
			echo "$(date) - nexus file for iqtree written." >> {params.wd}/results/report.txt
			if [[ -z "{params.maxmem}" ]]; then
				iqtree -S loci.nex --prefix loci -nt AUTO -ntmax {threads} -redo
			else
				iqtree -S loci.nex --prefix loci -nt AUTO -ntmax {threads} -redo -mem {params.maxmem}
			fi
		else
			echo "$(date) - phylociraptor will run iqtree gene trees now, with model testing as specified in the config.yaml file" >> {params.wd}/results/report.txt
			if [[ -z "{params.maxmem}" ]]; then
				iqtree -S algn/ --prefix loci -nt AUTO -ntmax {threads} -m MFP -redo
			else
				iqtree -S algn/ --prefix loci -nt AUTO -ntmax {threads} -m MFP -redo -mem {params.maxmem}
			fi
		fi
		rm -r algn
		cd {params.wd}
		touch {output.checkpoint}
		"""

rule astral_species_tree:
	input:
		trees = rules.iqtree_gene_trees.output.trees,
		checkpoint = rules.iqtree_gene_trees.output.checkpoint
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
