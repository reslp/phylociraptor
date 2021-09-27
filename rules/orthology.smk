configfile: "data/config.yaml"

import os
import glob

def get_assemblies(wildcards):
	sp = "{wildcards.species}".format(wildcards=wildcards)
	sp.replace(" ", "_")
	#print(sp)
	if os.path.isfile("results/assemblies/" + sp + ".fna"):
		return ["results/assemblies/" + sp + ".fna"]
	elif os.path.isfile("results/assemblies/" + sp + ".fna.gz"):
		return ["results/assemblies/" + sp + ".fna.gz"]
	else:
		return []

def select_species(dir="results/assemblies"):
	sps = [sp.split("/")[-1].split(".fna")[0] for sp in glob.glob(dir+"/*fna*")]
#	print("Species ("+str(len(sps))+"):"+str(sps))
	if config["exclude_orthology"]:
		blacklist = []
		with open(config["exclude_orthology"]) as file:
			blacklist = [line.strip() for line in file]
#		print("blacklist: "+str(blacklist))
		for i in reversed(range(len(sps))):
			if sps[i] in blacklist:
				del(sps[i])
#		print("Species ("+str(len(sps))+"): "+str(sps))
	return sps

species = select_species()

if config["busco"]["version"] == "3.0.2":
	rule run_busco:
		input:
			assembly = get_assemblies,
			busco_set = "results/orthology/busco/busco_set/"+config["busco"]["set"]
		output:
			checkpoint = "results/checkpoints/busco/busco_{species}.done",
			augustus_output = "results/orthology/busco/busco_runs/{species}/run_busco/augustus_output.tar.gz",
			blast_output = "results/orthology/busco/busco_runs/{species}/run_busco/blast_output.tar.gz",
			hmmer_output = "results/orthology/busco/busco_runs/{species}/run_busco/hmmer_output.tar.gz",
			full_table = "results/orthology/busco/busco_runs/{species}/run_busco/full_table_busco.tsv",
			short_summary ="results/orthology/busco/busco_runs/{species}/run_busco/short_summary_busco.txt",
			missing_busco_list ="results/orthology/busco/busco_runs/{species}/run_busco/missing_busco_list_busco.tsv",
			single_copy_buscos = "results/orthology/busco/busco_runs/{species}/run_busco/single_copy_busco_sequences.tar",
			single_copy_buscos_tarlist = "results/orthology/busco/busco_runs/{species}/run_busco/single_copy_busco_sequences.txt"

		benchmark: "results/statistics/benchmarks/busco/run_busco_{species}.txt"
		threads: int(config["busco"]["threads"])
		shadow: "shallow"
		log: "log/{species}_busco.log"
		params:
			wd = os.getcwd(),
			sp = config["busco"]["augustus_species"],
			additional_params = config["busco"]["additional_parameters"],
			species = lambda wildcards: wildcards.species,
			mode = config["busco"]["mode"],
			augustus_config_in_container = "/opt/conda/config",
			set = config["busco"]["set"]
		singularity:
			"docker://reslp/busco:3.0.2"
		shell:
			"""
			mkdir -p log
			dir=results/busco/{params.species}
			# prepare stripped down version auf augustus config path.
			# this is introduced to lower the number of files.
			mkdir augustus
			cp -R {params.augustus_config_in_container}/cgp augustus
			cp {params.augustus_config_in_container}/config.ini augustus
			cp -R {params.augustus_config_in_container}/extrinsic augustus
			cp -R {params.augustus_config_in_container}/model augustus
			cp -R {params.augustus_config_in_container}/profile augustus
			mkdir augustus/species
			
			if [ -d {params.augustus_config_in_container}/species/{params.sp} ]
			then
				cp -R {params.augustus_config_in_container}/species/{params.sp} augustus/species
			fi		

			export AUGUSTUS_CONFIG_PATH=$(pwd)/augustus
			echo $AUGUSTUS_CONFIG_PATH
			
			# handle gzipped and other assemblies differently
			if echo $(file $(readlink -f "{input.assembly}")) | grep compressed ;
			then
				fullname=$(basename "{input.assembly}")
				filename="${{fullname%.*}}"
				gunzip -c $(readlink -f "{input.assembly}") > "$filename"
			else
				filename="{input.assembly}"
			fi
			echo "Assembly used for BUSCO is "$filename 2>&1 | tee {log}
			run_busco -i $filename -f --out busco -c {threads} -sp {params.sp} --lineage_path {input.busco_set} -m {params.mode} {params.additional_params} 2>&1 | tee -a {log}

			# do some cleanup to save space
			echo -e "\\n[$(date)]\\tCleaning up after BUSCO to save space" 2>&1 | tee -a {log}
			bin/tar_folder.sh {output.blast_output} run_busco/blast_output 2>&1 | tee -a {log}
			bin/tar_folder.sh {output.hmmer_output} run_busco/hmmer_output 2>&1 | tee -a {log}
			bin/tar_folder.sh {output.augustus_output} run_busco/augustus_output 2>&1 | tee -a {log}
			tar -pcf {output.single_copy_buscos} -C run_busco/ single_copy_busco_sequences
			tar -tvf {output.single_copy_buscos} > {output.single_copy_buscos_tarlist} 2>&1 | tee -a {log}
			
			#move output files:
			mv run_busco/full_table_busco.tsv {output.full_table}
			mv run_busco/short_summary_busco.txt {output.short_summary}
			mv run_busco/missing_busco_list_busco.tsv {output.missing_busco_list}
			#mv run_busco/single_copy_busco_sequences {output.single_copy_buscos}
			
			buscos=$(tail -n +6 results/orthology/busco/busco_runs/{params.species}/run_busco/full_table_busco.tsv | cut -f 2 | sort | uniq -c | tr '\\n' ' ' | sed 's/ $/\\n/g')
			name="{params.species}"
			echo "$(date) $name $buscos" >> results/statistics/runlog.txt
			
			touch {output.checkpoint}
			"""

if config["busco"]["version"] == "5.2.1":
	rule run_busco:
		input:
			assembly = get_assemblies,
			busco_set = "results/orthology/busco/busco_set/"+config["busco"]["set"]
		output:
			checkpoint = "results/checkpoints/busco/busco_{species}.done",
			augustus_output = "results/orthology/busco/busco_runs/{species}/run_busco/augustus_output.tar.gz",
			blast_output = "results/orthology/busco/busco_runs/{species}/run_busco/blast_output.tar.gz",
			hmmer_output = "results/orthology/busco/busco_runs/{species}/run_busco/hmmer_output.tar.gz",
			logs = "results/orthology/busco/busco_runs/{species}/run_busco/logs.tar.gz",
			full_table = "results/orthology/busco/busco_runs/{species}/run_busco/full_table_busco.tsv",
			short_summary ="results/orthology/busco/busco_runs/{species}/run_busco/short_summary_busco.txt",
			missing_busco_list ="results/orthology/busco/busco_runs/{species}/run_busco/missing_busco_list_busco.tsv",
			single_copy_buscos = "results/orthology/busco/busco_runs/{species}/run_busco/single_copy_busco_sequences.tar",
			single_copy_buscos_tarlist = "results/orthology/busco/busco_runs/{species}/run_busco/single_copy_busco_sequences.txt"

		benchmark: "results/statistics/benchmarks/busco/run_busco_{species}.txt"
		threads: int(config["busco"]["threads"])
		shadow: "shallow"
		log: "log/{species}_busco.log"
		params:
			wd = os.getcwd(),
			sp = config["busco"]["augustus_species"],
			additional_params = config["busco"]["additional_parameters"],
			species = lambda wildcards: wildcards.species,
			mode = config["busco"]["mode"],
			augustus_config_in_container = "/usr/local/config",
			set = config["busco"]["set"]
		singularity:
			"docker://ezlabgva/busco:v5.2.1_cv1"
		shell:
			"""
			mkdir -p log
			dir=results/busco/{params.species}
			# prepare stripped down version auf augustus config path.
			# this is introduced to lower the number of files.
			mkdir augustus
			cp -R {params.augustus_config_in_container}/cgp augustus
			cp -R {params.augustus_config_in_container}/extrinsic augustus
			cp -R {params.augustus_config_in_container}/model augustus
			cp -R {params.augustus_config_in_container}/profile augustus
			mkdir augustus/species
			cp -R {params.augustus_config_in_container}/species/generic augustus/species/
			
			if [ -d {params.augustus_config_in_container}/species/{params.sp} ]
			then
				cp -R {params.augustus_config_in_container}/species/{params.sp} augustus/species
			fi		

			export AUGUSTUS_CONFIG_PATH=$(pwd)/augustus
			#export AUGUSTUS_SCRIPTS_PATH=/usr/local/bin/ #might be necessary in ezlabgva/busco:v5.2.1_cv1 image
			#export AUGUSTUS_BIN_PATH=/usr/local/bin/
			echo $AUGUSTUS_CONFIG_PATH
			
			# handle gzipped and other assemblies differently
			if [[ "{input.assembly}" =~ \.gz$ ]]
			then
				fullname=$(basename "{input.assembly}")
				filename="${{fullname%.*}}"
				gunzip -c $(readlink -f "{input.assembly}") > "$filename"
			else
				filename="{input.assembly}"
			fi
			echo "Assembly used for BUSCO is "$filename 2>&1 | tee {log}
			busco -i $filename -f --out {params.species} -c {threads} --augustus --augustus_species {params.sp} --lineage_dataset $(pwd)/{input.busco_set} -m {params.mode} {params.additional_params} 2>&1 | tee -a {log}
			# do some cleanup to save space
			echo -e "\\n[$(date)]\\tCleaning up after BUSCO to save space" 2>&1 | tee -a {log}
			basedir=$(pwd)
			cd {params.species}/run_{params.set}
			$basedir/bin/tar_folder.sh $basedir/{output.blast_output} blast_output 2>&1 | tee -a $basedir/{log}
			$basedir/bin/tar_folder.sh $basedir/{output.hmmer_output} hmmer_output 2>&1 | tee -a $basedir/{log}
			$basedir/bin/tar_folder.sh $basedir/{output.augustus_output} augustus_output 2>&1 | tee -a $basedir/{log}
			cd ..
			$basedir/bin/tar_folder.sh $basedir/{output.logs} logs 2>&1 | tee -a $basedir/{log}
			cd ..
			tar -pcf {output.single_copy_buscos} -C {params.species}/run_{params.set}/busco_sequences single_copy_busco_sequences 
			tar -tvf {output.single_copy_buscos} > {output.single_copy_buscos_tarlist} 2>&1 | tee -a $basedir/{log}

			#move output files:
			mv {params.species}/run_{params.set}/full_table.tsv {output.full_table}
			mv {params.species}/run_{params.set}/short_summary.txt {output.short_summary}
			mv {params.species}/run_{params.set}/missing_busco_list.tsv {output.missing_busco_list}
			
			buscos=$(tail -n +6 {output.full_table} | cut -f 2 | sort | uniq -c | tr '\\n' ' ' | sed 's/ $/\\n/g')
			name="{params.species}"
			echo "$(date) $name $buscos" >> results/statistics/runlog.txt
			
			#touch checkpoint
			touch {output.checkpoint}
			"""

rule busco:
	input:
		checks = expand("results/checkpoints/busco/busco_{species}.done", species=species)
	output:
		"results/checkpoints/busco.done"
	shell:
		"""
		touch {output}
		"""

rule extract_busco_table:
	input:
		busco_set = "results/orthology/busco/busco_set/"+config["busco"]["set"],
		busco = rules.busco.output
	output:
		busco_table = "results/orthology/busco/busco_table.txt",
		#checkpoint = "results/checkpoints/extract_busco_table.done"
	benchmark:
		"results/statistics/benchmarks/extract_busco_table.txt"
	singularity:
		"docker://reslp/biopython_plus:1.77"
	params:
		busco_dir = "results/orthology/busco/busco_runs/"
	shell:
		"""
		python bin/extract_busco_table.py --hmm {input.busco_set}/hmms --busco_results {params.busco_dir} -o {output.busco_table}
		echo "species\tcomplete\tsingle_copy\tduplicated\tfragmented\tmissing\ttotal" > results/statistics/busco_summary.txt
		for file in $(ls results/orthology/busco/busco_runs/*/run_busco/short_summary_busco.txt);
		do  
			name=$(echo $file | sed 's#results/busco/##' | sed 's#/run_busco/short_summary_busco.txt##'); 
			printf $name; cat $file | grep -P '\t\d' | awk -F "\t" '{{printf "\t"$2}}' | awk '{{print}}'; done >> results/statistics/busco_summary.txt
		"""
rule all_orthology:
	input:
		expand("results/checkpoints/busco/busco_{species}.done", species=species),
		#"results/checkpoints/busco.done",
		"results/orthology/busco/busco_table.txt",
		#"results/checkpoints/create_sequence_files.done",
		#"results/checkpoints/remove_duplicated_sequence_files.done"
	output:
		"results/checkpoints/modes/orthology.done"
	shell:
		"""
		touch {output}
		echo "$(date) - Pipeline part 1 (orthology) done." >> results/statistics/runlog.txt
		"""
