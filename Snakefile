import pandas as pd
#from pathlib import Path


configfile: "data/config.yaml"

wd = os.getcwd()

sample_data = pd.read_csv(config["species"]).set_index("species", drop=False)
samples = [sample.replace(" ", "_") for sample in sample_data["species"].tolist()]
print(samples)
	
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
		expand("results/checkpoints/busco_{sp}.done", sp=samples),
		"results/checkpoints/busco_table.done"

# needs to be run first before other rules can be run.
rule download_genomes:
	output:
		"results/checkpoints/download_genomes.done"	
	singularity:
		"docker://reslp/biomartr:0.9.2"
	params:
		species = get_species_names,
		wd = os.getcwd()
	shell:
		"""
		Rscript bin/genome_download.R {params.species} {params.wd}
		touch {output}
		"""		
		
rule download_busco_set:
	output:
		busco_set = directory("results/busco_set"),
		checkpoint = "results/checkpoints/download_busco_set.done"
	params:
		set = config["busco"]["set"]
	shell:
		"""
		echo {params.set}
		wget http://busco.ezlab.org/v2/datasets/{params.set}.tar.gz
		tar xfz {params.set}.tar.gz
		rm {params.set}.tar.gz
		if [ -d {output.busco_set} ]; then rm -rf {output.busco_set}; fi
		mv {params.set} {output.busco_set}
		touch {output.checkpoint}
		"""

rule prepare_augustus:
	output:
		augustus_config_path = directory("results/augustus_config_path"),
		checkpoint = "results/checkpoints/prepare_augustus.done"
	singularity:
		"docker://reslp/busco:3.0.2"
	shell:
		"""
		cp -r /opt/conda/config {output.augustus_config_path}
		touch {output.checkpoint}
		"""

rule run_busco:
	input:
		assembly = "results/downloaded_genomes/{species}_genomic_genbank.fna",
		augustus_config_path = rules.prepare_augustus.output.augustus_config_path,
		busco_set = rules.download_busco_set.output.busco_set,
	output:
		busco_dir = directory("results/busco/{species}"),
		checkpoint = "results/checkpoints/busco_{species}.done"
	threads: config["busco"]["threads"]
	params:
		wd = os.getcwd()
	singularity:
		"docker://reslp/busco:3.0.2"
	shell:
		"""
		if [ ! -d {output.busco_dir} ]; then mkdir {output.busco_dir}; fi
		cd {output.busco_dir}
		export AUGUSTUS_CONFIG_PATH={params.wd}/{input.augustus_config_path}
		echo $AUGUSTUS_CONFIG_PATH
		run_busco -i ../../../{input.assembly} --out busco -c {threads} --lineage_path ../../../{input.busco_set} -m genome
		cd ../../../
		touch {output.checkpoint}		
		"""
		
rule extract_hmms:
	input:
		busco_set = rules.download_busco_set.output.busco_set,
		busco_dir = directory("results/busco/")
	output:
		busco_table = "results/busco_table/busco_table.txt",
		checkpoint = "results/checkpoints/busco_table.done"
	shell:
		"""
		python bin/get_busco_table.py --hmm {input.busco_set}/hmms --busco_results {input.busco_dir} > {output.busco_table}
		touch {output.checkpoint}
		"""
		