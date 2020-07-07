import pandas as pd
configfile: "data/config.yaml"

wd = os.getcwd()

sample_data = pd.read_csv(config["species"]).set_index("species", drop=False)
samples = [sample.replace(" ", "_") for sample in sample_data["species"].tolist()]
#print(samples)
	
def get_species_names(wildcards):
	names = [name.replace(" ", "_") for name in sample_data["species"].to_list()]
	names= ",".join(names)
	print(names)
	return names

rule all:
	input:
		#expand("results/checkpoints/download_genome_{sample}.done", sample=samples)
		"results/checkpoints/download_genomes.done",
		"results/checkpoints/download_busco_set.done",
		"results/checkpoints/prepare_augustus.done",
		#expand("results/checkpoints/busco_{sample}.done", sample=samples),
		"results/checkpoints/busco_table.done",
		"results/checkpoints/create_sequence_files.done",
		"results/checkpoints/get_all_trimmed_files.done",
		"results/checkpoints/iqtree.done",
		"results/checkpoints/iqtree_gene_trees.done",
		"results/checkpoints/astral_species_tree.done"	
rule setup:
	input:
		"results/checkpoints/download_genomes.done",
		"results/checkpoints/download_busco_set.done",
		"results/checkpoints/prepare_augustus.done"
	output:
		".phylogenomics_setup.done"
	shell:
		"""
		touch {output}
		"""
rule part1:
	input:
		"results/checkpoints/extract_busco_table.done",
		"results/checkpoints/create_sequence_files.done"
	output:
		"checkpoints/part1.done"
	shell:
		"""
		touch {output}
		"""
rule part2:
	input:
		"results/trimmed_alignments/busco_list.txt"
	output:
		"checkpoints/part2.done"
	shell:
		"""
		touch {output}	
		"""

rule part3:
	input:
		"results/checkpoints/iqtree.done",
		"results/checkpoints/iqtree_gene_trees.done",
		"results/checkpoints/astral_species_tree.done"
	output:
		"checkpoints/part3.done"
	shell:
		"""
		touch {output}
		"""

include: "rules/setup.smk"
include: "rules/part1.smk"
include: "rules/part2.smk"
include: "rules/part3.smk"

