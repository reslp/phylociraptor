# needs to be run first before other rules can be run.
rule download_genomes:
	output:
		checkpoint = "results/checkpoints/download_genomes.done",
		download_overview = "results/downloaded_genomes/download_overview.txt",
		success = "results/downloaded_genomes/successfully_downloaded.txt",
		runlog = "results/statistics/runlog.txt"
	benchmark:
		"results/statistics/benchmarks/setup/download_genomes.txt"
	singularity:
		"docker://reslp/biomartr:0.9.2_exp"
	log:
		"log/download_genomes.log"
	params:
		species = get_species_names,
                wd = os.getcwd()
	shell:
		"""
		if [[ ! -f results/statistics/runlog.txt ]]; then touch results/statistics/runlog.txt; fi
		if [[ "{params.species}" != "" ]]; then
			echo "(date) - Setup: Will download species now" >> results/statistics/runlog.txt
			Rscript bin/genome_download.R {params.species} {params.wd} 2>&1 | tee {log}
		else
			echo "(date) - Setup: Now species to download." >> results/statistics/runlog.txt
		fi
		touch {output}
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
		for spe in $(cat {input.success}); do
			echo $spe
			if [[ -f {params.wd}/results/assemblies/"$spe".fna ]]; then
				continue
			else
				link=$(tail -n +2 "{input.overview}" | grep "^$spe," | awk -F',' '{{print $2}}')
				if [[ ! -f "$link" ]]; then
					echo "$spe" >> {output.statistics} 
					continue
				else
					ln -s $link {params.wd}/results/assemblies/"$spe".fna
				fi
			fi
		done	
		for spe in {params.local_species}; do
			sparr=(${{spe//,/ }})
			if [[ -f {params.wd}/results/assemblies/"${{sparr[0]}}".fna ]]; then
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
		busco_set = directory("results/busco_set"),
		checkpoint = "results/checkpoints/download_busco_set.done"
	params:
		set = config["busco"]["set"]
	benchmark:
		"results/statistics/benchmarks/setup/dowload_busco_set.txt"
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
		echo "file_name\torganism\turl\tdatabase\tpath\trefseq_category\tassembly_accession\tbioproject\tbiosample\ttaxid\tinfraspecific_name\tversion_status\trelease_type\tgenome_rep\tseq_rel_date\tsubmitter" > {output.statistics}	
			for file in $(ls results/downloaded_genomes/*_db_genbank.tsv); 
				do
					species=$(echo $file | awk -F 'doc_' '{{print $2}}' | awk -F '_db' '{{print $1}}')
					if [[ -f "results/downloaded_genomes/"$species"_genomic_genbank.fna" ]]; then 
						sed -n 2p $file >> {output.statistics}
					fi
				done
		"""

