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
		"results/checkpoints/create_busco_table.done",
		"results/checkpoints/create_sequence_files.done"
	output:
		"checkpoints/part1.done"
	shell:
		"""
		touch {output}
		"""


include: "rules/setup.smk"
include: "rules/part1.smk"

		
checkpoint create_sequence_files:
	input:
		busco_table = rules.extract_busco_table.output.busco_table,
		busco_dir = "results/busco/"
	output:
		sequence_dir = directory("results/busco_sequences"),
		checkpoint = "results/checkpoints/create_sequence_files.done"
	params:
		cutoff=config["extract_sequences"]["cutoff"],
		minsp=config["extract_sequences"]["minsp"]
	singularity:
		"docker://continuumio/miniconda3:4.7.10"
	conda:
		"envs/biopython.yml"
	shell:
		"""
		mkdir -p {output.sequence_dir}
		python bin/create_sequence_files.py --busco_table {input.busco_table} --busco_results {input.busco_dir} --cutoff {params.cutoff} --outdir {output.sequence_dir} --minsp {params.minsp}
		touch {output.checkpoint}
		"""

rule align:
	input:
		checkpoint = rules.create_sequence_files.output.checkpoint,
		sequence_file = "results/busco_sequences/{busco}_all.fas"
	output:
		alignment = "results/alignments/{busco}_aligned.fas",
		checkpoint = "results/checkpoints/mafft/{busco}_aligned.done"
	singularity:
		"docker://continuumio/miniconda3:4.7.10"
	conda:
		"envs/mafft.yml"
	threads:
		config["mafft"]["threads"]
	params:
		config["mafft"]["parameters"]
	shell:
		"""
		mafft {params} {input.sequence_file} > {output.alignment}
		touch {output.checkpoint}
		"""

checkpoint trim:
	input:
		rules.align.output.alignment
	output:
		trimmed_alignment = "results/trimmed_alignments/{busco}_aligned_trimmed.fas",
		checkpoint = "results/checkpoints/trimmed/{busco}_trimmed.done"
	singularity:
		"docker://continuumio/miniconda3:4.7.10"
	conda:
		"envs/trimal.yml"
	shell:
		"""
		trimal -gappyout -in {input} -out {output.trimmed_alignment}
		touch {output.checkpoint}
		"""

def aggregate_trimmed_alignments(wildcards):
    checkpoint_output = checkpoints.create_sequence_files.get(**wildcards).output[0]  
    file_names = expand("results/trimmed_alignments/{busco}_aligned_trimmed.fas", busco = glob_wildcards(os.path.join(checkpoint_output, "{busco}_all.fas")).busco)
    return file_names


rule get_all_trimmed_files:
	input:
		aggregate_trimmed_alignments
	output:
		checkpoint = "results/checkpoints/get_all_trimmed_files.done"
	shell:
		"""
		touch {output.checkpoint}
		"""
		
if config["phylogeny"]["concat"] == "yes":
	rule iqtree:
		input:
			rules.get_all_trimmed_files.output.checkpoint
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
			cp {params.wd}/results/trimmed_alignments/* algn
			iqtree -p algn/ --prefix concat -bb {params.bb} -nt {params.nt} -m {params.m} -redo -T {threads}
			rm -r algn
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

if config["phylogeny"]["species_tree"] == "yes":
	rule iqtree_gene_trees:
		input:
			rules.get_all_trimmed_files.output.checkpoint,
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
			cp {params.wd}/results/trimmed_alignments/* algn
			iqtree -S algn --prefix loci -nt AUTO -m WAG -redo -T {threads}
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
