import subprocess

BUSCOS, = glob_wildcards("results/busco_sequences/{busco}_all.fas")
BUSCOS = BUSCOS[1:10]
rule align:
        input:
                checkpoint = rules.create_sequence_files.output.checkpoint,
                sequence_file = "results/busco_sequences/{busco}_all.fas"
        output:
                alignment = "results/alignments/{busco}_aligned.fas",
                checkpoint = "results/checkpoints/mafft/{busco}_aligned.done"
        singularity:
                "docker://continuumio/miniconda3:4.7.10"
        conda:
                "../envs/mafft.yml"
        threads:
                config["mafft"]["threads"]
        params:
                config["mafft"]["parameters"]
        shell:
                """
                mafft {params} {input.sequence_file} > {output.alignment}
                touch {output.checkpoint}
                """

rule trim:
        input:
                rules.align.output.alignment
        output:
                trimmed_alignment = "results/trimmed_alignments/{busco}_aligned_trimmed.fas",
                checkpoint = "results/checkpoints/trimmed/{busco}_trimmed.done"
        singularity:
                "docker://continuumio/miniconda3:4.7.10"
        conda:
                "../envs/trimal.yml"
        shell:
                """
                trimal -gappyout -in {input} -out {output.trimmed_alignment}
                touch {output.checkpoint}
                """

#def aggregate_trimmed_alignments(wildcards):
#    checkpoint_output = checkpoints.create_sequence_files.get(**wildcards).output[0]
#    file_names = expand("results/trimmed_alignments/{busco}_aligned_trimmed.fas", busco = glob_wildcards(os.path.join(checkpoint_output, "{busco}_all.fas")).busco)
#    return file_names

rule get_all_trimmed_files:
        input:
                expand("results/trimmed_alignments/{bus}_aligned_trimmed.fas", bus=BUSCOS)
        output:
                "results/trimmed_alignments/busco_list.txt"
	shell:
                """
		for (file in {input}); do
                	echo $(basename $file) >> {output}
		done
                """

