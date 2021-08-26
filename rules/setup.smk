import pandas as pd
configfile: "data/config.yaml"

sample_data = pd.read_csv(config["species"]).set_index("species", drop=False)
samples = [sample.replace(" ", "_") for sample in sample_data["species"].tolist()]

def get_species_names(wildcards):
	names = [name.replace(" ", "_") for name in sample_data.loc[sample_data["web_local"].str.startswith("web"), "species"].to_list()]
	acc = sample_data.loc[sample_data["web_local"].str.startswith("web"), "web_local"].to_list()
	full = []
	for i in range(len(names)):
		if len(acc[i].split("=")) > 1:
			full.append(names[i]+"="+acc[i].split("=")[1])
		else:
			full.append(names[i])
	full= ",".join(full)
	return full

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

# needs to be run first before other rules can be run.
#rule download_genomes:
#	output:
#		checkpoint = "results/checkpoints/download_genomes.done",
#		download_overview = "results/downloaded_genomes/download_overview.txt",
#		success = "results/downloaded_genomes/successfully_downloaded.txt",
#		runlog = "results/statistics/runlog.txt"
#	benchmark:
#		"results/statistics/benchmarks/setup/download_genomes.txt"
#	singularity:
#		"docker://reslp/biomartr:0.9.2_exp"
#	log:
#		"log/download_genomes.log"
#	params:
#		species = get_species_names,
#                wd = os.getcwd()
#	shell:
#		"""
#		if [[ ! -f results/statistics/runlog.txt ]]; then touch results/statistics/runlog.txt; fi
#		if [[ "{params.species}" != "" ]]; then
#			echo "(date) - Setup: Will download species now" >> results/statistics/runlog.txt
#			Rscript bin/genome_download.R {params.species} {params.wd} 2>&1 | tee {log}
#		else
#			echo "(date) - Setup: No species to download." >> results/statistics/runlog.txt
#		fi
#		touch {output}
#		"""

rule download_genomes:
	output:
		checkpoint = "results/checkpoints/download_genomes.done",
		download_overview = "results/downloaded_genomes/download_overview.txt",
		success = "results/downloaded_genomes/successfully_downloaded.txt",	
		runlog = "results/statistics/runlog.txt"
	benchmark:
		"results/statistics/benchmarks/setup/download_genomes.txt"
	singularity: "docker://reslp/biopython_plus:1.77"
	log:
		"log/download_genomes.log"
	params:
		species = get_species_names,
		wd = os.getcwd(),
		email = config["email"]
	shell:
		"""
		if [[ ! -f results/statistics/runlog.txt ]]; then touch results/statistics/runlog.txt; fi
		if [[ "{params.species}" != "" ]]; then
			echo "(date) - Setup: Will download species now" >> results/statistics/runlog.txt
			python bin/genome_download.py --entrez_email {params.email} --outdir {params.wd}/results/downloaded_genomes/ --genomes {params.species} 2>&1 | tee {log}
		else
			# need to touch these files, since they are usually produced by the python script.
			touch {output.success}
			touch {output.download_overview}
			echo "(date) - Setup: No species to download." >> results/statistics/runlog.txt
		fi
		touch {output.checkpoint}
		"""

rule rename_assemblies:
	input:
		success = rules.download_genomes.output.success,
		overview = rules.download_genomes.output.download_overview
	output:
		checkpoint = "results/checkpoints/rename_assemblies.done",
		statistics = "results/statistics/species_not_downloaded.txt",
		statistics_local = "results/statistics/local_species.txt"
	benchmark:
		"results/statistics/benchmarks/setup/rename_assemblies.txt"
	params:
		#downloaded_species = get_species_names_rename,
		local_species = get_local_species_names_rename,
		wd = os.getcwd()
	shell:
		"""
		mkdir -p results/assemblies
		for spe in $(grep ",success" {input.overview}); do
			sp=$(echo $spe | awk -F',' '{{print $1}}')
			if [[ -f {params.wd}/results/assemblies/"$sp".fna.gz ]]; then
				continue
			else
				link="{params.wd}/results/downloaded_genomes/"$sp"_genomic_genbank.fna.gz"
				echo $link
				if [[ ! -f "$link" ]]; then
					echo "$sp" >> {output.statistics} 
					continue
				else
					ln -s $link {params.wd}/results/assemblies/"$sp".fna.gz
				fi
			fi
		done	
		for spe in {params.local_species}; do
			sparr=(${{spe//,/ }})
			if [[ -L {params.wd}/results/assemblies/"${{sparr[0]}}".fna ]]; then
				echo "${{sparr[0]}}" >> {output.statistics_local}
				continue
			else
				echo "${{sparr[0]}}" >> {output.statistics_local}
				ln -s {params.wd}/"${{sparr[1]}}" {params.wd}/results/assemblies/"${{sparr[0]}}".fna
			fi
		done
		if [[ ! -f {output.statistics} ]]; then touch {output.statistics}; fi
		if [[ ! -f {output.statistics_local} ]]; then touch {output.statistics_local}; fi
		touch {output.checkpoint}
		"""



rule download_busco_set:
	output:
		busco_set = directory("results/orthology/busco/busco_set"),
		checkpoint = "results/checkpoints/download_busco_set.done"
	params:
		set = config["busco"]["set"]
	benchmark:
		"results/statistics/benchmarks/setup/dowload_busco_set.txt"
	shell:
		"""
		echo {params.set}
		if [ -d {output.busco_set} ]; then rm -rf {output.busco_set}; fi
		if [ ! -d {output.busco_set} ]; then mkdir {output.busco_set}; fi
		wget -c {params.set} -O - | tar -xz --strip-components 1 -C {output.busco_set}/
		touch {output.checkpoint}
		"""

rule prepare_augustus:
	output:
		augustus_config_path = directory("results/augustus_config_path"),
		checkpoint = "results/checkpoints/prepare_augustus.done"
	singularity:
		"docker://reslp/busco:3.0.2"
	benchmark:
		"results/statistics/benchmarks/setup/prepare_augustus.txt"
	shell:
		"""
		cp -r /opt/conda/config {output.augustus_config_path}
		touch {output.checkpoint}
		"""

rule get_genome_download_statistics:
	input:
		checkpoint = "results/checkpoints/download_genomes.done"
	output:
		statistics = "results/statistics/downloaded_genomes_statistics.txt"
	benchmark:
		"results/statistics/setup/get_genome_download_statistics.txt"
	shell:
		"""
		echo "name\tid\tassembly_accession\tbioproject\tbiosample\twgs_master\trefseq_category\ttaxid\tspecies_taxid\torganism_name\tinfraspecific_name\tisolate\tversion_status\tassembly_level\trelease_type\tgenome_rep\tseq_rel_date\tasm_name\tsubmitter\tgbrs_paired_asm\tpaired_asm_comp\tftp_path\texcluded_from_refseq\trelation_to_type_material" > {output.statistics}
			for file in $(ls results/downloaded_genomes/*_db_genbank.tsv); 
				do
					species=$(echo $file | awk -F '_db_' '{{print $1}}' | awk -F '/' '{{print $(NF)}}')
					echo $species
					if [[ -f "results/downloaded_genomes/"$species"_genomic_genbank.fna.gz" ]]; then 
						output=$(sed -n 2p $file)
						echo $species"\t""$output" >> {output.statistics}
					fi
				done
		"""
rule setup:
	input:
		"results/checkpoints/download_genomes.done",
		"results/checkpoints/download_busco_set.done",
		"results/checkpoints/prepare_augustus.done",
		"results/checkpoints/rename_assemblies.done",
		"results/statistics/downloaded_genomes_statistics.txt"
		#expand("results/assemblies/{species}.fna", species=samples)
	output:
		".phylogenomics_setup.done"

	shell:
		"""
		touch {output}
		mkdir -p results/statistics
		touch "results/statistics/runlog.txt"
		echo "$(date) - Pipeline setup done." >> results/statistics/runlog.txt
		"""

rule add_genomes:
	input:
		"results/checkpoints/download_genomes.done",
		"results/checkpoints/rename_assemblies.done",
		"results/statistics/downloaded_genomes_statistics.txt"
	output:
		".add_genomes.done"
	shell:
		"""
		touch {output}
		"""
