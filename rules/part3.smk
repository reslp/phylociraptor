if config["phylogeny"]["concat"] == "yes":
	rule iqtree:
		input:
			rules.part2.output
		output:
			checkpoint = "results/checkpoints/iqtree.done"
		singularity:
			"docker://reslp/iqtree:2.0rc2"
		params:
			wd = os.getcwd(),
			nt = "AUTO",
			bb = "1000",
			m = "WAG"
		threads:
			config["iqtree"]["threads"]
		shell:
			"""
			rm -rf results/phylogeny/concatenated/algn
			mkdir -p results/phylogeny/concatenated/
			cd results/phylogeny/concatenated/
			mkdir algn
			cp {params.wd}/results/filtered_alignments/*.fas algn
			iqtree -p algn/ --prefix concat -bb {params.bb} -nt {params.nt} -m {params.m} -redo -T {threads}
			rm -r algn
			cd {params.wd}
			touch {output.checkpoint}
			"""
	rule concatenate:
		input:
			rules.part2.output
		output:
			checkpoint = "results/checkpoints/concatenate.done",
			alignment = "results/phylogeny/raxmlng/concat.fas"
		params:
			wd = os.getcwd(),
			ids = config["species"]
		singularity:
			"docker://reslp/concat:0.2"
		shell:
			"""
			tail -n +2 {params.ids} | awk -F "," '{{print $1;}}' | sed 's/ /_/g' > results/phylogeny/raxmlng/ids.txt	
			concat.py -d results/filtered_alignments/ -t results/phylogeny/raxmlng/ids.txt --runmode concat -o results/phylogeny/raxmlng/
			touch {output.checkpoint}
			"""	
	rule raxmlng:
		input:
			rules.concatenate.output.alignment
		output:
			checkpoint = "results/checkpoints/raxmlng.done"
		singularity:
			"docker://reslp/raxml-ng:1.0.0"
		params:
			threads = config["raxmlng"]["threads"],
			bs = config["raxmlng"]["bootstrap"],
			wd = os.getcwd(),
			model = config["raxmlng"]["model"]
		shell:
			"""
			cd results/phylogeny/raxmlng
			raxml-ng --msa concat.fas --prefix raxmlng -threads {params.threads} --bs-trees {params.bs} --model {params.model} --all
			cd {params.wd}
			touch {output.checkpoint}
			"""

else:  #checkpoint files need to be created anyway
	rule iqtree:
		output:
			checkpoint = "results/checkpoints/iqtree.done"
		shell:
			"""
			touch {output.checkpoint}
			"""
	rule concatenate:
		output:
			checkpoint = "results/checkpoints/concatenate.done",
			alignment = "results/phylogeny/raxmlng/concat.fas"
		shell:
			"""
			touch {output.checkpoint}
			touch {output.alignment}
			"""
	rule raxmlng:
		output: 
			checkpoint = "results/checkpoints/raxmlng.done"
		shell:
			"""
			touch {output.checkpoint}
			"""

if config["phylogeny"]["species_tree"] == "yes":
	rule iqtree_gene_trees:
		input:
			rules.part2.output,
		output:
			checkpoint = "results/checkpoints/iqtree_gene_trees.done",
			trees = "results/phylogeny/gene_trees/loci.treefile"
		params:
			wd = os.getcwd(),
		threads:
			config["iqtree"]["threads"]
		singularity:
			"docker://reslp/iqtree:2.0rc2"
		shell:
			"""
			rm -rf results/phylogeny/gene_trees/algn
			mkdir -p results/phylogeny/gene_trees/algn
			cd results/phylogeny/gene_trees
			cp {params.wd}/results/filtered_alignments/*.fas algn
			iqtree -S algn/ --prefix loci -nt AUTO -m WAG -redo -T {threads}
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
			
else: #checkpoint files need to be created anyway
	rule iqtree_gene_trees:
		output:
			checkpoint = "results/checkpoints/iqtree_gene_trees.done"
		shell:
			"""
			touch {output.checkpoint}
			"""
	rule astral_species_tree:
		output:
			checkpoint = "results/checkpoints/astral_species_tree.done"
		shell:
			"""
			touch {output.checkpoint}
			"""

