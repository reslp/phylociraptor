rule run_busco:
        input:
                assembly = "results/assemblies/{species}.fna",
                augustus_config_path = rules.prepare_augustus.output.augustus_config_path,
                busco_set = rules.download_busco_set.output.busco_set,
        output:
                checkpoint = "results/checkpoints/busco_{species}.done",
		augustus_output = "results/busco/{species}/run_busco/augustus_output.tar.gz",
		blast_output = "results/busco/{species}/run_busco/blast_output.tar.gz",
		hmmer_output = "results/busco/{species}/run_busco/hmmer_output.tar.gz",
		full_table = "results/busco/{species}/run_busco/full_table_busco.tsv",
		short_summary ="results/busco/{species}/run_busco/short_summary_busco.txt",
		missing_busco_list ="results/busco/{species}/run_busco/missing_busco_list_busco.tsv",
		single_copy_buscos = directory("results/busco/{species}/run_busco/single_copy_busco_sequences")
        threads: config["busco"]["threads"]
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
		echo "$(date) $name $buscos" >> results/report.txt
		
		#touch checkpoint
		touch {output.checkpoint}
                """

rule busco:
        input:
                checks = expand("results/checkpoints/busco_{sample}.done", sample=samples)
        output:
                "checkpoints/busco.done"
        shell:
                """
                touch {output}
                """


rule extract_busco_table:
        input:
                busco_set = rules.download_busco_set.output.busco_set,
                busco = rules.busco.output
        output:
                busco_table = "results/busco_table/busco_table.txt",
                checkpoint = "results/checkpoints/extract_busco_table.done"
        singularity:
                "docker://continuumio/miniconda3:4.7.10"
        params:
                busco_dir = "results/busco/"
        shell:
                """
                python bin/extract_busco_table.py --hmm {input.busco_set}/hmms --busco_results {params.busco_dir} > {output.busco_table}
                touch {output.checkpoint}
                """

rule create_sequence_files:
        input:
                busco_table = rules.extract_busco_table.output.busco_table,
                checkpoint = rules.busco.output
        output:
                sequence_dir = directory("results/busco_sequences"),
                checkpoint = "results/checkpoints/create_sequence_files.done"
        params:
                cutoff=config["filtering"]["cutoff"],
                minsp=config["filtering"]["minsp"],
		busco_dir = "results/busco"
        singularity:
                "docker://continuumio/miniconda3:4.7.10"
        conda:
                "../envs/biopython.yml"
        shell:
                """
		mkdir -p {output.sequence_dir}
                python bin/create_sequence_files.py --busco_table {input.busco_table} --busco_results {params.busco_dir} --cutoff {params.cutoff} --outdir {output.sequence_dir} --minsp {params.minsp}
		touch {output.checkpoint}
                
		"""

rule remove_duplicated_sequence_files:
        input:
                checkpoint = rules.create_sequence_files.output.checkpoint
        output:
                checkpoint = "results/checkpoints/remove_duplicated_sequence_files.done"
        singularity:
                "docker://continuumio/miniconda3:4.7.10"
        conda:
                "../envs/biopython.yml"
        params:
                wd = os.getcwd(),
		dupseq=config["filtering"]["dupseq"]
        shell:
                """
                if [[ -d results/busco_sequences_deduplicated ]]; then
                        rm -rf results/busco_sequences_deduplicated
                fi
                mkdir results/busco_sequences_deduplicated

		if [[ {params.dupseq} == "persample" ]];
                	then
                		echo "$(date) - BUSCO files will be filtered on a per-sample basis. This could lower the number of species in the final tree." >> results/report.txt
			else
				echo "$(date) - BUSCO files will be filtered on a per BUSCO gene basis. This could lower the number of genes used to calculate the final tree." >> results/report.txt
		fi

                for file in results/busco_sequences/*.fas;
                do
			if [[ {params.dupseq} == "persample" ]];
			then
				# per sequence filtering
				python bin/filter_alignments.py --alignments {params.wd}/$file --outdir "{params.wd}/results/busco_sequences_deduplicated" --per_sample
			else
				# whole alignment filtering
                        	python bin/filter_alignments.py --alignments {params.wd}/$file --outdir "{params.wd}/results/busco_sequences_deduplicated"
                	fi
		done
                echo "$(date) - Number of BUSCO sequence files: $(ls results/busco_sequences/*.fas | wc -l)" >> results/report.txt
                echo "$(date) - Number of deduplicated BUSCO sequence files: $(ls results/busco_sequences_deduplicated/*.fas | wc -l)" >> results/report.txt
                touch {output.checkpoint}
                """
