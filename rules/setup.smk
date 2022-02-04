import pandas as pd
import numpy as np
import yaml

sample_data = pd.read_csv(config["species"]).set_index("species", drop=False)
samples = [sample.replace(" ", "_") for sample in sample_data["species"].tolist()]

# get list of containers to use:
with open("data/containers.yaml", "r") as yaml_stream:
    containers = yaml.safe_load(yaml_stream)

def determine_concurrency_limit():
	names = [name.replace(" ", "_") for name in sample_data.loc[sample_data["web_local"].str.startswith("web"), "species"].to_list()]
	if len(names) < config["concurrency"]:
		return range(1, len(names)+1)
	else:
		return range(1, config["concurrency"]+1)

batches = determine_concurrency_limit()

def get_species_names(wildcards):
	names = [name.replace(" ", "_") for name in sample_data.loc[sample_data["web_local"].str.startswith("web"), "species"].to_list()]
	acc = sample_data.loc[sample_data["web_local"].str.startswith("web"), "web_local"].to_list()
	nbatches = config["concurrency"]
	#final_names = [names[i:i+nbatches] for i in range(0, len(names), nbatches)]
	#final_acc = [acc[i:i+nbatches] for i in range(0, len(acc), nbatches)]
	final_names = [list(n) for n in np.array_split(names, nbatches)]
	final_acc = [list(n) for n in np.array_split(acc, nbatches)]
	full = []
	final_names = final_names[int(wildcards.batch)-1]
	final_acc = final_acc[int(wildcards.batch)-1]
	for i in range(len(final_names)):
		if len(final_acc[i].split("=")) > 1:
			full.append(final_names[i]+"="+final_acc[i].split("=")[1])
		else:
			full.append(final_names[i])
	full= ",".join(full)
	return full

def get_species_names_rename(wildcards):
	names = [name.replace(" ", "_") for name in sample_data.loc[sample_data["web_local"] == "web", "species"].to_list()]
	names= " ".join(names)
	return names


def get_local_species_names_rename(wildcards):
	names_with_accession = [name for name in sample_data.loc[sample_data["web_local"].str.startswith("web="), "species"].to_list()]
	names_with_web = [name for name in sample_data.loc[sample_data["web_local"] == "web", "species"].to_list()]
	all_names = [name for name in sample_data["species"].to_list()]
	
	# get all names of web and web= species
	all_online_names = names_with_accession + names_with_web
	# select only the ones which do not have web or web=
	names = [x for x in all_names if x not in all_online_names]
	# now get assembly paths for those
	assembly_paths = sample_data[sample_data["species"].isin(names)]["web_local"].to_list()
	
	# get them into output format:	
	name_asse_pair = [name.replace(" ", "_") + "," + ass for name, ass in zip(names,assembly_paths)]
	names= " ".join(name_asse_pair)
	print("Local names: ", names)
	return names
        
def get_all_species_names(wildcards):
	print("All species names:")
	names = [name for name in sample_data.loc["species"].to_list()]
	#print(names)
	return names


rule download_genome_overview:
	output:
		overview = "results/downloaded_genomes/assembly_summary_genbank.txt"
	singularity: containers["biopython"]
	log:
		"log/download_genomes_overview.log"
	threads: 2
	params:
		wd = os.getcwd(),
		email = config["email"]
	shell:
		"""
		python bin/genome_download.py --entrez_email {params.email} --outdir {params.wd}/results/downloaded_genomes/ --overview-only &> {log}
		"""

rule download_genomes:
	input:
		config["species"],
		rules.download_genome_overview.output.overview
	output:
		checkpoint = "results/checkpoints/genome_download/download_genomes_{batch}-"+str(config["concurrency"])+".done",
		download_overview = "results/downloaded_genomes/download_overview_{batch}.txt",
		success = "results/downloaded_genomes/successfully_downloaded_{batch}.txt",
		failed = "results/downloaded_genomes/not_downloaded_{batch}.txt"
#	benchmark:
#		"results/statistics/benchmarks/setup/download_genomes.txt"
	singularity: containers["biopython"]
	log:
		"log/download_genomes_{batch}-"+str(config["concurrency"])+".log"
	params:
		species = get_species_names,
		wd = os.getcwd(),
		email = config["email"],
		batch = "{batch}"
	shell:
		"""
		if [[ ! -f results/statistics/runlog.txt ]]; then if [[ ! -d results/statistics ]]; then mkdir -p results/statistics; fi; touch results/statistics/runlog.txt; fi
		if [[ "{params.species}" != "" ]]; then
			echo "$(date) - Setup: Will download species now - batch: {params.batch}" >> results/statistics/runlog.txt
			python bin/genome_download.py --entrez_email {params.email} --outdir {params.wd}/results/downloaded_genomes/ --genomes {params.species} --batch {params.batch} 2>&1 | tee {log}
		else
			# need to touch these files, since they are usually produced by the python script.
			touch {output.success}
			touch {output.download_overview}
			touch {output.failed}
			echo "$(date) - Setup: No species to download." >> results/statistics/runlog.txt
		fi
		if [[ ! -d results/checkpoints/genome_download ]]; then mkdir -p results/checkpoints/genome_download; fi # sometimes this folder is not created, better do it to be safe.
		touch {output.checkpoint}
		"""

rule get_genome_download_statistics:
	input:
		checkpoints = expand(rules.download_genomes.output.checkpoint, batch=batches)
	output:
		statistics = "results/statistics/downloaded_genomes_statistics.txt",
		success = "results/downloaded_genomes/successfully_downloaded.txt",
		failed = "results/downloaded_genomes/not_downloaded.txt",
		overview = "results/downloaded_genomes/download_overview.txt"
	params:
		checkpoints = expand("results/checkpoints/genome_download/download_genomes_{batch}-"+str(config["concurrency"])+".done", batch=batches),
		success = expand("results/downloaded_genomes/successfully_downloaded_{b}.txt", b=batches),
		overview = expand("results/downloaded_genomes/download_overview_{b}.txt", b=batches),
		failed = expand("results/downloaded_genomes/not_downloaded_{b}.txt", b=batches)
		
	benchmark:
		"results/statistics/setup/get_genome_download_statistics.txt"
	shell:
		"""
		if [[ -z "{input}" ]]; then
			echo "All genomes local. Will touch empty output files."
			echo "name\tid\tassembly_accession\tbioproject\tbiosample\twgs_master\trefseq_category\ttaxid\tspecies_taxid\torganism_name\tinfraspecific_name\tisolate\tversion_status\tassembly_level\trelease_type\tgenome_rep\tseq_rel_date\tasm_name\tsubmitter\tgbrs_paired_asm\tpaired_asm_comp\tftp_path\texcluded_from_refseq\trelation_to_type_material" > {output.statistics}
			touch {output.success}
			touch {output.failed}
			touch {output.overview}
		else
			cat {params.success} > {output.success}
			cat {params.failed} > {output.failed}	
			cat {params.overview} > {output.overview}
			echo "name\tid\tassembly_accession\tbioproject\tbiosample\twgs_master\trefseq_category\ttaxid\tspecies_taxid\torganism_name\tinfraspecific_name\tisolate\tversion_status\tassembly_level\trelease_type\tgenome_rep\tseq_rel_date\tasm_name\tsubmitter\tgbrs_paired_asm\tpaired_asm_comp\tftp_path\texcluded_from_refseq\trelation_to_type_material" > {output.statistics}
			for file in $(ls results/downloaded_genomes/*_db_genbank.tsv); 
				do
					echo "---- Adding info ----"
					echo "file: "$file
					species=$(echo $file | awk -F '_db_' '{{print $1}}' | awk -F '/' '{{print $(NF)}}')
					
					echo "species: "$species
					if [[ -f "results/downloaded_genomes/"$species"_genomic_genbank.fna.gz" ]]; then 
						output=$(sed -n 2p $file)
						echo $species"\t""$output" >> {output.statistics}
					fi
				done
		fi
		"""


rule rename_assemblies:
	input:
		overview=rules.get_genome_download_statistics.output.overview
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
			if [[ -f {params.wd}/results/assemblies/"$sp".fasta.gz ]]; then
				continue
			else
				link="{params.wd}/results/downloaded_genomes/"$sp"_genomic_genbank.fna.gz"
				echo $link
				if [[ ! -f "$link" ]]; then
					echo "$sp" >> {output.statistics} 
					continue
				else
					ln -s $link {params.wd}/results/assemblies/"$sp".fasta.gz
				fi
			fi
		done	
		for spe in {params.local_species}; do
			sparr=(${{spe//,/ }})
			if [[ -L {params.wd}/results/assemblies/"${{sparr[0]}}".fasta ]]; then
				echo "${{sparr[0]}}" >> {output.statistics_local}
				continue
			else
				echo "${{sparr[0]}}" >> {output.statistics_local}
				ln -s {params.wd}/"${{sparr[1]}}" {params.wd}/results/assemblies/"${{sparr[0]}}".fasta
			fi
		done
		if [[ ! -f {output.statistics} ]]; then touch {output.statistics}; fi
		if [[ ! -f {output.statistics_local} ]]; then touch {output.statistics_local}; fi
		touch {output.checkpoint}
		"""



rule download_busco_set:
	output:
		busco_set = directory("results/orthology/busco/busco_set/"+config["busco"]["set"]),
		checkpoint = "results/checkpoints/download_busco_set.done"
	params:
		set = config["busco"]["set"],
		busco_version = config["busco"]["version"]
	log:
		"log/download_busco_set.log"
	benchmark:
		"results/statistics/benchmarks/setup/dowload_busco_set.txt"
	shell:
		"""
		echo -e "[$(date)]\\tBUSCO set specified: {params.set}" 2>&1 | tee {log}
		if [ -d {output.busco_set} ]; then rm -rf {output.busco_set}; fi
		mkdir {output.busco_set}

		if [ "{params.busco_version}" == "3.0.2" ]
		then
			base_url="https://busco.ezlab.org/v3/datasets"
			echo -e "[$(date)]\\tDownloading .." 2>&1 | tee -a {log}
			wget -q -c $base_url/{params.set}.tar.gz -O - --no-check-certificate | tar -xz --strip-components 1 -C {output.busco_set}/
		elif [ "{params.busco_version}" == "5.2.1" ]
		then
			base_url="https://busco-data.ezlab.org/v5/data/lineages"
			current=$(curl -s $base_url/ | grep "{params.set}" | cut -d ">" -f 2 | sed 's/<.*//')
			echo -e "[$(date)]\\tCurrent version is: $current" 2>&1 | tee -a {log}
			echo -e "[$(date)]\\tDownloading .." 2>&1 | tee -a {log}
			wget -q -c $base_url/$current -O - --no-check-certificate | tar -xz --strip-components 1 -C {output.busco_set}/
		else
			echo -e "\\n######\\nPlease specify a valid BUSCO version in your config file - supported are '3.0.2' and '5.0.2'\\n######" 2>&1 | tee -a {log}
			exit 1
		fi

		echo -ne "[$(date)]\\tDone!\\n" 2>&1 | tee -a {log}
		touch {output.checkpoint}
		"""


rule setup:
	input:
		expand("results/checkpoints/genome_download/download_genomes_{batch}-"+str(config["concurrency"])+".done", batch=batches),
		"results/checkpoints/download_busco_set.done",
		"results/checkpoints/rename_assemblies.done",
		"results/statistics/downloaded_genomes_statistics.txt"
		#expand("results/assemblies/{species}.fna", species=samples)
	output:
		"results/checkpoints/modes/phylogenomics_setup.done"

	shell:
		"""
		touch {output}
		mkdir -p results/statistics
		touch "results/statistics/runlog.txt"
		echo "$(date) - phylociraptor setup done." >> results/statistics/runlog.txt
		"""

rule add_genomes:
	input:
		"results/checkpoints/download_genomes.done",
		"results/checkpoints/rename_assemblies.done",
		"results/statistics/downloaded_genomes_statistics.txt"
	output:
		"results/checkpoints/modes/add_genomes.done"
	shell:
		"""
		touch {output}
		"""

