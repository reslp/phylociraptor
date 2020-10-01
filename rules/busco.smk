rule run_busco:
        input:
                assembly = "results/assemblies/{species}.fna",
                augustus_config_path = rules.prepare_augustus.output.augustus_config_path,
                busco_set = rules.download_busco_set.output.busco_set,
        output:
                busco_dir = directory("results/busco/{species}"),
                checkpoint = "results/checkpoints/busco_{species}.done"
        threads: config["busco"]["threads"]
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
		if [ ! -d {output.busco_dir} ]; then mkdir {output.busco_dir}; fi
                cd {output.busco_dir}
                export AUGUSTUS_CONFIG_PATH={params.wd}/{input.augustus_config_path}
                echo $AUGUSTUS_CONFIG_PATH
                run_busco -i ../../../{input.assembly} --out busco -c {threads} -sp {params.sp} --lineage_path ../../../{input.busco_set} -m genome {params.additional_params}
                cd ../../../
                
		buscos=$(tail -n +6 results/busco/{params.species}/run_busco/full_table_busco.tsv | cut -f 2 | sort | uniq -c | tr '\\n' ' ' | sed 's/ $/\\n/g')
		name="{params.species}"
		echo "$(date) $name $buscos" >> results/report.txt

		touch {output.checkpoint}
                """

rule busco:
        input:
                checks = expand("results/checkpoints/busco_{sample}.done", sample=samples),
		dirs = expand("results/busco/{sample}", sample=samples)
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
                cutoff=config["extract_sequences"]["cutoff"],
                minsp=config["extract_sequences"]["minsp"],
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
                wd = os.getcwd()
        shell:
                """
                if [[ -d results/busco_sequences_deduplicated ]]; then
                        rm -rf results/busco_sequences_deduplicated
                fi
                mkdir results/busco_sequences_deduplicated

                for file in results/busco_sequences/*.fas;
                do
                        python bin/filter_alignments.py --alignments {params.wd}/$file --outdir "{params.wd}/results/busco_sequences_deduplicated"
                done
                echo "$(date) - Number of BUSCO sequence files: $(ls results/busco_sequences/*.fas | wc -l)" >> results/report.txt
                echo "$(date) - Number of deduplicated BUSCO sequence files: $(ls results/busco_sequences_deduplicated/*.fas | wc -l)" >> results/report.txt
                touch {output.checkpoint}
                """
