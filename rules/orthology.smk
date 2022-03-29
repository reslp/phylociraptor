import os
import glob
import pandas as pd
import yaml

# read in CSV file again. This is needed to determine the mode for BUSCO
sample_data = pd.read_csv(config["species"])
sample_data["species"] = sample_data["species"].str.replace(" ","_")
sample_data.set_index("species", drop=False)
samples_from_csv = sample_data["species"].to_list() # needed for crosscheck if taxa are removed or renamed in csv file. 
print(sample_data["species"].to_list())

# get list of containers to use:
with open("data/containers.yaml", "r") as yaml_stream:
    containers = yaml.safe_load(yaml_stream)

def get_busco_mode(wildcards):
	sp = "{wildcards.species}".format(wildcards=wildcards)
	sp.replace(" ", "_")
	#first check of mode column exists:
	if "mode" in list(sample_data.columns):
		print("(phylociraptor): INFO: Will use custom provided BUSCO mode given in .csv file for species:", sp)
		if sample_data.loc[sample_data["species"] == sp, "mode"].isnull().values.any(): #check if value was actually provided
			print("(phylocirpator): No BUSCO mode value was provided for" , sp , ", will use the default:",config["busco"]["mode"]) 
			if not config["busco"]["version"].startswith("5.") and config["busco"]["mode"] == "transcriptome":
				print("phylociraptor: ERROR: Incompatible parameters. Transcriptome mode only workes with BUSCO 5. Please check your config files. Exiting.")
				sys.exit(1)
			elif not config["busco"]["version"].startswith("5.") and config["busco"]["mode"] == "protein":
				print("phylociraptor: ERROR: Incompatible parameters. Protein mode only workes with BUSCO 5. Please check your config files. Exiting.")
				sys.exit(1)
			else:	
				return config["busco"]["mode"]
		else:
			mode = sample_data.loc[sample_data["species"] == sp, "mode"].to_string(index=False)
			if mode == "transcriptome" and not config["busco"]["version"].startswith("5."):
				print("phylociraptor: ERROR: Incompatible parameters. Transcriptome mode only workes with BUSCO 5. Please check your config files. Exiting.")
				sys.exit(1)
			elif mode == "protein" and not config["busco"]["version"].startswith("5."):
				print("phylociraptor: ERROR: Incompatible parameters. Protein mode only workes with BUSCO 5. Please check your config files. Exiting.")
				sys.exit(1)
			else:		
				return mode 
	else:
		print("phylocirpator: INFO: mode column not found, will use globally set option from config file.")		
		return config["busco"]["mode"]
	

def get_assemblies(wildcards):
	sp = "{wildcards.species}".format(wildcards=wildcards)
	sp.replace(" ", "_")
	if os.path.isfile("results/assemblies/" + sp + ".fasta"):
		return ["results/assemblies/" + sp + ".fasta"]
	elif os.path.isfile("results/assemblies/" + sp + ".fasta.gz"):
		return ["results/assemblies/" + sp + ".fasta.gz"]
	else:
		return []

def select_species(dir="results/assemblies"):
	sps = [sp.split("/")[-1].split(".fasta")[0] for sp in glob.glob(dir+"/*fasta*")]
#	print("Species ("+str(len(sps))+"):"+str(sps))
	sps = list(set(sps).intersection(samples_from_csv)) # crosscheck with CSV file to see if taxa have been removed
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
print(species)
if config["busco"]["version"] == "3.0.2":
	rule busco:
		input:
			assembly = get_assemblies,
			busco_set = "results/orthology/busco/busco_set/"+config["busco"]["set"]
		output:
			checkpoint = "results/checkpoints/busco/busco_{species}.done",
			output = "results/orthology/busco/busco_runs/{species}/run_busco/software_output.tar.gz",
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
			mode = get_busco_mode,
			augustus_config_in_container = "/opt/conda/config",
			set = config["busco"]["set"]
		singularity:
			containers["busco3"]
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
			#echo $AUGUSTUS_CONFIG_PATH
			
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
			basedir=$(pwd)
			cd run_busco
			mkdir software_outputs
			mv *_output software_outputs
			$basedir/bin/tar_folder.sh $basedir/{output.output} software_outputs 2>&1 | tee -a $basedir/{log}
			cd ..
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
	rule busco:
		input:
			assembly = get_assemblies,
			busco_set = "results/orthology/busco/busco_set/"+config["busco"]["set"]
		output:
			checkpoint = "results/checkpoints/busco/busco_{species}.done",
			output = "results/orthology/busco/busco_runs/{species}/run_busco/software_output.tar.gz",
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
			mode = get_busco_mode, 
			augustus_config_in_container = "/usr/local/config",
			set = config["busco"]["set"]
		singularity:
			containers["busco5"]
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
			#echo $AUGUSTUS_CONFIG_PATH
			
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
			busco -i $filename -f --out {params.species} -c {threads} $(if [[ "{params.sp}" != "None" ]]; then echo "--augustus --augustus_species {params.sp}"; fi) --lineage_dataset $(pwd)/{input.busco_set} -m {params.mode} {params.additional_params} 2>&1 | tee -a {log}
			# do some cleanup to save space
			echo -e "\\n[$(date)]\\tCleaning up after BUSCO to save space" 2>&1 | tee -a {log}
			basedir=$(pwd)
			cd {params.species}/run_{params.set}
			mkdir software_outputs
			mv *_output software_outputs
			$basedir/bin/tar_folder.sh $basedir/{output.output} software_outputs 2>&1 | tee -a $basedir/{log}
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

rule aggregate_orthology:
	input:
		checks = expand("results/checkpoints/busco/busco_{species}.done", species=species)
	output:
		"results/checkpoints/busco.done"
	shell:
		"""
		touch {output}
		"""

rule extract_orthology_table:
	input:
		busco_set = "results/orthology/busco/busco_set/"+config["busco"]["set"],
		busco = rules.aggregate_orthology.output
	output:
		busco_table = "results/orthology/busco/busco_table.txt",
		#checkpoint = "results/checkpoints/extract_busco_table.done"
	benchmark:
		"results/statistics/benchmarks/extract_busco_table.txt"
	singularity:
		containers["biopython"]
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
rule orthology:
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
