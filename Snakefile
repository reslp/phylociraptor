import pandas as pd
configfile: "data/config.yaml"

wd = os.getcwd()

sample_data = pd.read_csv(config["species"]).set_index("species", drop=False)
samples = [sample.replace(" ", "_") for sample in sample_data["species"].tolist()]
print("Samples which will be analysed")
print(samples)
	
def get_species_names(wildcards):
	names = [name.replace(" ", "_") for name in sample_data.loc[sample_data["web_local"] == "web", "species"].to_list()]
	names= ",".join(names)
	return names

def get_species_names_rename(wildcards):
        names = [name.replace(" ", "_") for name in sample_data.loc[sample_data["web_local"] == "web", "species"].to_list()]
        names= " ".join(names)
        return names


def get_local_species_names_rename(wildcards):
	names = [name.replace(" ","_") for name in sample_data.loc[sample_data["web_local"] != "web", "species"].to_list()]
	assembly = [name for name in sample_data.loc[sample_data["web_local"] != "web", "web_local"].to_list()]
	assembly = [name for name in assembly if str(name) != "nan"]
	name_asse_pair = [name + "," + ass for name, ass in zip(names,assembly)]
	names= " ".join(name_asse_pair)
	#print(names)
	return names

def get_all_species_names(wildcards):
	print("All species names:")
	names = [name for name in sample_data.loc["species"].to_list()]
	#print(names)
	return names

rule all:
	input:
		#expand("results/checkpoints/download_genome_{sample}.done", sample=samples)
		#"results/checkpoints/download_genomes.done",
		#"results/checkpoints/download_busco_set.done",
		#"results/checkpoints/prepare_augustus.done",
		#expand("results/checkpoints/busco_{sample}.done", sample=samples),
		#"results/checkpoints/busco_table.done",
		#"results/checkpoints/create_sequence_files.done",
		#"results/checkpoints/get_all_trimmed_files.done",
		#"results/checkpoints/iqtree.done",
		#"results/checkpoints/iqtree_gene_trees.done",
		#"results/checkpoints/astral_species_tree.done"
		".phylogenomics_setup.done",
		"checkpoints/part1.done",
		"checkpoints/part2.done",
		"checkpoints/part3.done"
rule setup:
	input:
		"results/checkpoints/download_genomes.done",
		"results/checkpoints/download_busco_set.done",
		"results/checkpoints/prepare_augustus.done",
		"results/checkpoints/rename_assemblies.done",
		expand("results/assemblies/{species}.fna", species=samples)
	output:
		".phylogenomics_setup.done"

	shell:
		"""
		touch {output}
		touch "results/report.txt"
		echo "$(date) - Pipeline setup done." >> results/report.txt
		"""
rule part1:
	input:
		expand("results/checkpoints/busco_{species}.done", species=samples),
		"results/checkpoints/extract_busco_table.done",
		"results/checkpoints/create_sequence_files.done",
		"results/checkpoints/remove_duplicated_sequence_files.done"
	output:
		"checkpoints/part1.done"
	shell:
		"""
		touch {output}
		echo "$(date) - Pipeline part 1 (busco) done." >> results/report.txt
		"""
rule part2:
	input:
		"results/checkpoints/get_all_trimmed_files.done"
	output:
		"checkpoints/part2.done"
	shell:
		"""
		touch {output}	
		echo "$(date) - Pipeline part 2 (align) done." >> results/report.txt
		"""

rule part3:
	input:
		"results/checkpoints/iqtree.done",
		"results/checkpoints/iqtree_gene_trees.done",
		"results/checkpoints/astral_species_tree.done",
		"results/checkpoints/raxmlng.done"
	output:
		"checkpoints/part3.done"
	shell:
		"""
		touch {output}
		echo "$(date) - Pipeline part 3 (tree) done." >> results/report.txt
		"""

include: "rules/setup.smk"
include: "rules/busco.smk"
include: "rules/align_trim.smk"
include: "rules/tree.smk"
include: "rules/report.smk"
