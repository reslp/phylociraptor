# needs to be run first before other rules can be run.
rule download_genomes:
        output:
                "results/checkpoints/download_genomes.done"
        singularity:
                "docker://reslp/biomartr:0.9.2"
        params:
                species = get_species_names,
                wd = os.getcwd()
        shell:
                """
                Rscript bin/genome_download.R {params.species} {params.wd}
                touch {output}
                """

rule download_busco_set:
        output:
                busco_set = directory("results/busco_set"),
                checkpoint = "results/checkpoints/download_busco_set.done"
        params:
                set = config["busco"]["set"]
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
        shell:
                """
                cp -r /opt/conda/config {output.augustus_config_path}
                touch {output.checkpoint}
                """
