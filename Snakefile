import pandas as pd
#from pathlib import Path


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
		"results/checkpoints/busco_table.done",
		"results/checkpoints/create_sequence_files.done",
		"results/checkpoints/get_all_trimmed_files.done",
		"results/checkpoints/iqtree.done"

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
		
rule extract_busco_table:
	input:
		busco_set = rules.download_busco_set.output.busco_set,
		busco_dir = "results/busco/"
	output:
		busco_table = "results/busco_table/busco_table.txt",
		checkpoint = "results/checkpoints/busco_table.done"
	singularity:
		"docker://continuumio/miniconda3:4.7.10"
	shell:
		"""
		python bin/extract_busco_table.py --hmm {input.busco_set}/hmms --busco_results {input.busco_dir} > {output.busco_table}
		touch {output.checkpoint}
		"""
		
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
		8
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
