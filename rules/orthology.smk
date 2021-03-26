configfile: "data/config.yaml"
rule run_busco:
	input:
		assembly = "results/assemblies/{species}.fna",
		augustus_config_path = "results/augustus_config_path",
		busco_set = "results/busco_set",
	output:
		checkpoint = "results/checkpoints/busco/busco_{species}.done",
		augustus_output = "results/busco/{species}/run_busco/augustus_output.tar.gz",
		blast_output = "results/busco/{species}/run_busco/blast_output.tar.gz",
		hmmer_output = "results/busco/{species}/run_busco/hmmer_output.tar.gz",
		full_table = "results/busco/{species}/run_busco/full_table_busco.tsv",
		short_summary ="results/busco/{species}/run_busco/short_summary_busco.txt",
		missing_busco_list ="results/busco/{species}/run_busco/missing_busco_list_busco.tsv",
		single_copy_buscos = directory("results/busco/{species}/run_busco/single_copy_busco_sequences")
	benchmark: "results/statistics/benchmarks/busco/run_busco_{species}.txt"
	threads: int(config["busco"]["threads"])
	shadow: "shallow"
	params:
		wd = os.getcwd(),
		sp = config["busco"]["augustus_species"],
		additional_params = config["busco"]["additional_parameters"],
		species = lambda wildcards: wildcards.species
	singularity:
		"docker://reslp/busco:3.0.2"
	shell:
		"""
		mkdir -p log
		dir=results/busco/{params.species}
		# prepare stripped down version auf augustus config path.
		# this is introduced to lower the number of files.
		mkdir augustus
		cp -R /opt/conda/config/cgp augustus
		cp /opt/conda/config/config.ini augustus
		cp -R /opt/conda/config/extrinsic augustus
		cp -R /opt/conda/config/model augustus
		cp -R /opt/conda/config/profile augustus
		mkdir augustus/species
		
		if [ -d /opt/conda/config/species/{params.sp} ]
		then
			cp -R /opt/conda/config/species/{params.sp} augustus/species
		fi		

		export AUGUSTUS_CONFIG_PATH=$(pwd)/augustus
		echo $AUGUSTUS_CONFIG_PATH
		run_busco -i {input.assembly} -f --out busco -c {threads} -sp {params.sp} --lineage_path {input.busco_set} -m genome {params.additional_params}
		# do some cleanup to save space
		bin/tar_folder.sh {output.blast_output} run_busco/blast_output
		bin/tar_folder.sh {output.hmmer_output} run_busco/hmmer_output
		bin/tar_folder.sh {output.augustus_output} run_busco/augustus_output

		
		#move output files:
		#mv run_busco/augustus_output.tar.gz {output.augustus_output}
		#mv run_busco/blast_output.tar.gz {output.blast_output}
		#mv run_busco/hmmer_output.tar.gz {output.hmmer_output}
		mv run_busco/full_table_busco.tsv {output.full_table}
		mv run_busco/short_summary_busco.txt {output.short_summary}
		mv run_busco/missing_busco_list_busco.tsv {output.missing_busco_list}
		mv run_busco/single_copy_busco_sequences {output.single_copy_buscos}
		
		buscos=$(tail -n +6 results/busco/{params.species}/run_busco/full_table_busco.tsv | cut -f 2 | sort | uniq -c | tr '\\n' ' ' | sed 's/ $/\\n/g')
		name="{params.species}"
		echo "$(date) $name $buscos" >> results/statistics/runlog.txt
		
		#touch checkpoint
		touch {output.checkpoint}
		"""

rule busco:
	input:
		checks = expand("results/checkpoints/busco/busco_{species}.done", species=glob_wildcards("results/assemblies/{species}.fna").species)
	output:
		"results/checkpoints/busco.done"
	shell:
		"""
		touch {output}
		"""

rule extract_busco_table:
	input:
		busco_set = "results/busco_set",
		busco = rules.busco.output
	output:
		busco_table = "results/busco_table/busco_table.txt",
		checkpoint = "results/checkpoints/extract_busco_table.done"
	benchmark:
		"results/statistics/benchmarks/busco/extract_busco_table.txt"
	singularity:
		"docker://reslp/biopython_plus:1.77"
	params:
		busco_dir = "results/busco/"
	shell:
		"""
		python bin/extract_busco_table.py --hmm {input.busco_set}/hmms --busco_results {params.busco_dir} -o {output.busco_table}
		echo "species\tcomplete\tsingle_copy\tduplicated\tfragmented\tmissing\ttotal" > results/statistics/busco_summary.txt
		for file in $(ls results/busco/*/run_busco/short_summary_busco.txt);do  name=$(echo $file | sed 's#results/busco/##' | sed 's#/run_busco/short_summary_busco.txt##'); printf $name; cat $file | grep -P '\t\d' | awk -F "\t" '{{printf "\t"$2}}' | awk '{{print}}'; done >> results/statistics/busco_summary.txt
		touch {output.checkpoint}
		"""
rule orthology:
	input:
		# expand("results/checkpoints/busco/busco_{species}.done", species=glob_wildcards("results/assemblies/{species}.fna").species),
		"results/checkpoints/busco.done",
		"results/checkpoints/extract_busco_table.done",
		#"results/checkpoints/create_sequence_files.done",
		#"results/checkpoints/remove_duplicated_sequence_files.done"
	output:
		"checkpoints/orthology.done"
	shell:
		"""
		touch {output}
		echo "$(date) - Pipeline part 1 (orthology) done." >> results/statistics/runlog.txt
		"""
