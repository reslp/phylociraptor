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

rule busco:
        input:
                expand("results/checkpoints/busco_{sample}.done", sample=samples)
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
                checkpoint = "results/checkpoints/busco_table.done"
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
                busco_dir = rules.busco.output
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

