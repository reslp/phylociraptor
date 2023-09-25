import os
import glob
import pandas as pd
import yaml

include: "functions.smk"

# read in CSV file again. This is needed to determine the mode for BUSCO
sample_data = pd.read_csv(config["species"])
sample_data["species"] = sample_data["species"].str.replace(" ","_")
sample_data.set_index("species", drop=False)
samples_from_csv = sample_data["species"].to_list() # needed for crosscheck if taxa are removed or renamed in csv file. 
print(sample_data["species"].to_list())

# get list of containers to use:
with open("data/containers.yaml", "r") as yaml_stream:
    containers = yaml.safe_load(yaml_stream)

#determine whether some species in the csv file have specified a busco mode that is not the default according to the config file
def get_busco_mode(wildcards):
	sp = "{wildcards.species}".format(wildcards=wildcards)
	sp.replace(" ", "_")
	#first check of mode column exists:
	if "mode" in list(sample_data.columns):
		print("(phylociraptor): INFO: Will check for custom provided BUSCO mode given in .csv file for species:", sp)
		if sample_data.loc[sample_data["species"] == sp, "mode"].isnull().values.any(): #check if value was actually provided
			print("(phylocirpator): No BUSCO mode value was provided for" , sp , "- will use the default:",config["orthology"]["busco_options"]["mode"]) 
			if not config["orthology"]["busco_options"]["version"].startswith("5.") and config["orthology"]["busco_options"]["mode"] == "transcriptome":
				print("phylociraptor: ERROR: Incompatible parameters. Transcriptome mode only workes with BUSCO 5. Please check your config files. Exiting.")
				sys.exit(1)
			elif not config["orthology"]["busco_options"]["version"].startswith("5.") and config["orthology"]["busco_options"]["mode"] == "protein":
				print("phylociraptor: ERROR: Incompatible parameters. Protein mode only workes with BUSCO 5. Please check your config files. Exiting.")
				sys.exit(1)
			else:	
				return config["orthology"]["busco_options"]["mode"]
		else:
			mode = sample_data.loc[sample_data["species"] == sp, "mode"].to_string(index=False)
			if mode == "transcriptome" and not config["orthology"]["busco_options"]["version"].startswith("5."):
				print("phylociraptor: ERROR: Incompatible parameters. Transcriptome mode only workes with BUSCO 5. Please check your config files. Exiting.")
				sys.exit(1)
			elif mode == "protein" and not config["orthology"]["busco_options"]["version"].startswith("5."):
				print("phylociraptor: ERROR: Incompatible parameters. Protein mode only workes with BUSCO 5. Please check your config files. Exiting.")
				sys.exit(1)
			else:		
				return mode 
	else:
		print("phylocirpator: INFO: mode column not found, will use default busco mode as specified config file for", sp)		
		return config["orthology"]["busco_options"]["mode"]
	

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
	if config["orthology"]["exclude"]:
		excludelist = []
		with open(config["orthology"]["exclude"]) as file:
			excludelist = [line.strip() for line in file]
#		print("excludelist: "+str(excludelist))
		for i in reversed(range(len(sps))):
			if sps[i] in excludelist:
				del(sps[i])
#		print("Species ("+str(len(sps))+"): "+str(sps))

	# get default settings for busco
	default = get_hash("", "orthology,method orthology,busco_options,set orthology,busco_options,version orthology,busco_options,mode orthology,busco_options,augustus_species orthology,busco_options,additional_parameters", configfi, returndict=True) 
	default_hash = get_hash(default)
	hashes = []
	if not "mode" in sample_data:
		hashes = [str(default_hash)] * len(sps)
	else:
		for sp in sps:
			if sample_data.loc[sample_data["species"] == sp, "mode"].isnull().values.any():
				print(sp," -> default mode: ",sample_data.loc[sample_data["species"] == sp, "mode"].to_string(index=False))
				hashes.append(str(default_hash))
			else:
				print(sp," -> mode from data file: ",sample_data.loc[sample_data["species"] == sp, "mode"].to_string(index=False))
				nondefault = default.copy()
				nondefault["orthology"]["busco_options"]["mode"] = sample_data.loc[sample_data["species"] == sp, "mode"].to_string(index=False)
				hashes.append(get_hash(nondefault))
	return [sps, hashes]

species = select_species()
	

#get hash for current step
hashes = collect_hashes("orthology", config, configfi, wd=os.getcwd())
current_hash = hashes["orthology"]["global"]

###############
rule read_params_global:
	input:
		compare("results/orthology/busco/parameters.orthology."+current_hash+".yaml", configfi, optional=[config['species']])
	output:
		"results/orthology/busco/parameters.orthology."+current_hash+".yaml"
	shell:
		"""
		bin/read_write_yaml.py {input[0]} {output} species orthology,method orthology,exclude orthology,busco_options,set orthology,busco_options,version orthology,busco_options,mode orthology,busco_options,augustus_species orthology,busco_options,additional_parameters
		"""

if config["orthology"]["busco_options"]["version"] == "3.0.2":
	rule busco:
		input:
			assembly = get_assemblies,
			busco_set = "results/orthology/busco/busco_set/"+config["orthology"]["busco_options"]["set"],
		output:
			checkpoint = "results/checkpoints/busco."+str(config["orthology"]["busco_options"]["set"])+"/busco_{species}.{hash}.done",
			output = "results/orthology/busco/busco_runs."+str(config["orthology"]["busco_options"]["set"])+"/{species}.{hash}/run_busco/software_output.tar.gz",
			full_table = "results/orthology/busco/busco_runs."+str(config["orthology"]["busco_options"]["set"])+"/{species}.{hash}/run_busco/full_table_busco.tsv",
			short_summary ="results/orthology/busco/busco_runs."+str(config["orthology"]["busco_options"]["set"])+"/{species}.{hash}/run_busco/short_summary_busco.txt",
			missing_busco_list ="results/orthology/busco/busco_runs."+str(config["orthology"]["busco_options"]["set"])+"/{species}.{hash}/run_busco/missing_busco_list_busco.tsv",
			single_copy_buscos = "results/orthology/busco/busco_runs."+str(config["orthology"]["busco_options"]["set"])+"/{species}.{hash}/run_busco/single_copy_busco_sequences.tar",
			single_copy_buscos_tarlist = "results/orthology/busco/busco_runs."+str(config["orthology"]["busco_options"]["set"])+"/{species}.{hash}/run_busco/single_copy_busco_sequences.txt",

		benchmark: "results/statistics/benchmarks/busco."+str(config["orthology"]["busco_options"]["set"])+"/run_busco_{species}.{hash}.txt"
		threads: int(config["orthology"]["threads"])
		shadow: "shallow"
		log: "log/orthology/{species}.{hash}_busco."+str(config["orthology"]["busco_options"]["set"])+".log"
		params:
			wd = os.getcwd(),
			sp = config["orthology"]["busco_options"]["augustus_species"],
			additional_params = config["orthology"]["busco_options"]["additional_parameters"],
			species = lambda wildcards: wildcards.species,
			mode = get_busco_mode,
			augustus_config_in_container = "/opt/conda/config",
			set = config["orthology"]["busco_options"]["set"]
		singularity:
			containers["busco3"]
		shell:
			"""
			# prepare stripped down version auf augustus config path.
			# this is introduced to lower the number of files.
			mkdir augustus
			cp -R {params.augustus_config_in_container}/cgp augustus
			cp {params.augustus_config_in_container}/config.ini augustus
			cp -R {params.augustus_config_in_container}/extrinsic augustus
			cp -R {params.augustus_config_in_container}/model augustus
			cp -R {params.augustus_config_in_container}/profile augustus
			mkdir augustus/species
			
			if [ -z "{params.sp}" ] #check if augustus training species parameter is provided correctly
			then
				echo "ERROR: augustus_species field in config YAML appears to be empty. Please check and provide a correct name."
				exit 1
			elif [ "{params.sp}" == "None" ]
			then	
				echo "ERROR: augustus_species field in config YAML appears to be empty. Please check and provide a correct name."
				exit 1
			elif [ -d {params.augustus_config_in_container}/species/{params.sp} ] # now check if directory is present
			then
				cp -R {params.augustus_config_in_container}/species/{params.sp} augustus/species
			else
				echo "ERROR: Species {params.sp} as set in yaml config file could not be found within pretrained Augustus species. Please recheck that the name is correct!"
				exit 1
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
			
			buscos=$(tail -n +6 {output.full_table} | cut -f 2 | sort | uniq -c | tr '\\n' ' ' | sed 's/ $/\\n/g')
			name="{params.species}"
			echo "$(date) $name $buscos" >> results/statistics/runlog.txt
			
			touch {output.checkpoint}
			"""

if config["orthology"]["busco_options"]["version"] == "5.2.1":
	rule busco:
		input:
			assembly = get_assemblies,
			busco_set = "results/orthology/busco/busco_set/"+config["orthology"]["busco_options"]["set"],
		output:
			checkpoint = "results/checkpoints/busco."+str(config["orthology"]["busco_options"]["set"])+"/busco_{species}.{hash}.done",
			output = "results/orthology/busco/busco_runs."+str(config["orthology"]["busco_options"]["set"])+"/{species}.{hash}/run_busco/software_output.tar.gz",
			logs = "results/orthology/busco/busco_runs."+str(config["orthology"]["busco_options"]["set"])+"/{species}.{hash}/run_busco/logs.tar.gz",
			full_table = "results/orthology/busco/busco_runs."+str(config["orthology"]["busco_options"]["set"])+"/{species}.{hash}/run_busco/full_table_busco.tsv",
			short_summary ="results/orthology/busco/busco_runs."+str(config["orthology"]["busco_options"]["set"])+"/{species}.{hash}/run_busco/short_summary_busco.txt",
			missing_busco_list ="results/orthology/busco/busco_runs."+str(config["orthology"]["busco_options"]["set"])+"/{species}.{hash}/run_busco/missing_busco_list_busco.tsv",
			single_copy_buscos = "results/orthology/busco/busco_runs."+str(config["orthology"]["busco_options"]["set"])+"/{species}.{hash}/run_busco/single_copy_busco_sequences.tar",
			multi_copy_buscos = "results/orthology/busco/busco_runs."+str(config["orthology"]["busco_options"]["set"])+"/{species}.{hash}/run_busco/multi_copy_busco_sequences.tar",
			single_copy_buscos_tarlist = "results/orthology/busco/busco_runs."+str(config["orthology"]["busco_options"]["set"])+"/{species}.{hash}/run_busco/single_copy_busco_sequences.txt",
			fragmented_buscos = "results/orthology/busco/busco_runs."+str(config["orthology"]["busco_options"]["set"])+"/{species}.{hash}/run_busco/fragmented_busco_sequences.tar"

		benchmark: "results/statistics/benchmarks/busco."+str(config["orthology"]["busco_options"]["set"])+"/run_busco_{species}.{hash}.txt"
		threads: int(config["orthology"]["threads"])
		shadow: "shallow"
		log: "log/orthology/{species}.{hash}_busco."+str(config["orthology"]["busco_options"]["set"])+".log"
		params:
			wd = os.getcwd(),
			sp = config["orthology"]["busco_options"]["augustus_species"],
			additional_params = config["orthology"]["busco_options"]["additional_parameters"],
			species = lambda wildcards: wildcards.species,
			mode = get_busco_mode, 
			augustus_config_in_container = "/usr/local/config",
			set = config["orthology"]["busco_options"]["set"]
		singularity:
			containers["busco5"]
		shell:
			"""
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
			
			# handle gzipped 
			if [[ "{input.assembly}" =~ \.gz$ ]]
			then
				fullname=$(basename "{input.assembly}")
				filename="${{fullname%.*}}"
				gunzip -c $(readlink -f "{input.assembly}") > "$filename"
			else
				filename="{input.assembly}"
			fi
			
			echo "Assembly used for BUSCO is "$filename 2>&1 | tee {log}
			if [[ "{params.sp}" ]] && [[ "{params.sp}" != "None" ]]
			then
				augustus="--augustus --augustus_species {params.sp}"
				echo -e "\\n[$(date)]\\tUsing Augustus with model '{params.sp}'" 2>&1 | tee -a {log}
			else
				augustus=""
			fi

			busco -i $filename -f --out {params.species} -c {threads} $augustus --lineage_dataset $(pwd)/{input.busco_set} -m {params.mode} {params.additional_params} 2>&1 | tee -a {log}

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
			ls -1 {params.species}/run_{params.set}/busco_sequences
			tar -pcf {output.single_copy_buscos} -C {params.species}/run_{params.set}/busco_sequences single_copy_busco_sequences 
			tar -tvf {output.single_copy_buscos} > {output.single_copy_buscos_tarlist} 2>&1 | tee -a $basedir/{log}
			tar -pcf {output.multi_copy_buscos} -C {params.species}/run_{params.set}/busco_sequences multi_copy_busco_sequences
			tar -pcf {output.fragmented_buscos} -C {params.species}/run_{params.set}/busco_sequences fragmented_busco_sequences

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
		checks = expand("results/checkpoints/busco."+str(config["orthology"]["busco_options"]["set"])+"/busco_{species}.{hash}.done", zip, species=species[0], hash=species[1]),
		params = rules.read_params_global.output
	output:
		check = "results/checkpoints/busco."+config["orthology"]["busco_options"]["set"]+"."+current_hash+".done",
		dir = directory("results/orthology/busco/busco_runs."+config["orthology"]["busco_options"]["set"]+"."+current_hash)
	params:
		set = config["orthology"]["busco_options"]["set"]
	shell:
		"""
		mkdir {output.dir}
		for c in {input.checks}
		do
			sp_mo=$(basename $c | sed 's/^busco_//' | sed 's/\.done$//')
			ln -s ../../../../results/orthology/busco/busco_runs.{params.set}/$sp_mo {output.dir}/
		done
		touch {output.check}
		"""

rule extract_orthology_table:
	input:
		busco_set = "results/orthology/busco/busco_set/"+config["orthology"]["busco_options"]["set"],
		busco_dir = rules.aggregate_orthology.output.dir
	output:
		busco_table = "results/orthology/orthology_table."+current_hash+".txt",
	benchmark:
		"results/statistics/benchmarks/extract_orthology_table."+current_hash+".txt"
	singularity:
		containers["biopython"]
	shell:
		"""
		python bin/extract_busco_table.py --hmm {input.busco_set}/hmms --busco_results {input.busco_dir}/ -o {output.busco_table}
		echo "species\tcomplete\tsingle_copy\tduplicated\tfragmented\tmissing\ttotal" > results/statistics/busco_summary.txt
		for file in $(ls {input.busco_dir}/*/run_busco/short_summary_busco.txt);
		do  
			name=$(echo $file | sed 's#results/busco/##' | sed 's#/run_busco/short_summary_busco.txt##'); 
			printf $name; cat $file | grep -P '\t\d' | awk -F "\t" '{{printf "\t"$2}}' | awk '{{print}}'; done >> results/statistics/busco_summary.txt
		"""
rule orthology:
	input:
		"results/orthology/orthology_table."+current_hash+".txt",
	output:
		"results/checkpoints/modes/orthology."+current_hash+".done"
	shell:
		"""
		touch {output}
		echo "$(date) - phylociraptor orthology done." >> results/statistics/runlog.txt
		"""

