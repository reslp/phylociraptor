import subprocess

BUSCOS, = glob_wildcards("results/busco_sequences_deduplicated/{busco}_all.fas")

rule align:
        input:
                checkpoint = rules.remove_duplicated_sequence_files.output.checkpoint,
                sequence_file = "results/busco_sequences_deduplicated/{busco}_all.fas"
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

if config["trimming"]["method"] == "trimal":
	rule trim:
		input:
			rules.align.output.alignment
		output:
			trimmed_alignment = "results/trimmed_alignments/{busco}_aligned_trimmed.fas",
			checkpoint = "results/checkpoints/trimmed/{busco}_trimmed.done"
		params:
			trimmer = config["trimming"]["parameters"]
		singularity:
			"docker://continuumio/miniconda3:4.7.10"
		conda:
			"../envs/trimal.yml"
		shell:
			"""
			trimal {params.trimmer} -in {input} -out {output.trimmed_alignment}
			touch {output.checkpoint}
			"""
elif config["trimming"]["method"] == "aliscore":
	rule trim:
		input:
			rules.align.output.alignment
		output:
			trimmed_alignment = "results/trimmed_alignments/{busco}_aligned_trimmed.fas",
			checkpoint = "results/checkpoints/trimmed/{busco}_trimmed.done"
		params:
			trimmer = config["trimming"]["parameters"],
			busco = "{busco}",
			wd = os.getcwd()
		singularity:
			"docker://chrishah/alicut-aliscore-docker:2.31"
		shell:
			"""
			mkdir -p results/trimmed_alignments/{params.busco}
			cd "results/trimmed_alignments/{params.busco}"
			ln -s -f  {params.wd}/{input} {params.busco}_aligned.fas 
			Aliscore.pl {params.trimmer} {params.busco}_aligned.fas &> aliscore_{params.busco}.log || true
			
			if [[ -f {params.busco}_aligned.fas_List_random.txt ]]; then
				echo "$(date) - The aliscore output file does not exist. Check results for BUSCO: {params.busco}" >> {params.wd}/results/report.txt
				if [[ $(cat {params.busco}_aligned.fas_List_random.txt | head -n 1 | grep '[0-9]' -c) != 0 ]]; then
					echo "$(date) - The aliscore output appears to be empty. Check results for BUSCO: {params.busco}" >> {params.wd}/results/report.txt
					ALICUT.pl -s &> alicut_{params.busco}.log
				fi
			fi
			
			# this check is included because alicut very rarely does not  produce an output file.
			# in this case an empty file will be touched. This is necessary so the rule does not fail
			# The empty file will later be exluded again in the next rule.
			if [[ ! -f {params.wd}/results/trimmed_alignments/{params.busco}/ALICUT_{params.busco}_aligned.fas ]]; then
				echo "$(date) - The ALICUT output appears to be empty. Will touch an empty file so the pipeline will continue. Check results for BUSCO: {params.busco}" >> {params.wd}/results/report.txt
				touch {params.wd}/{output.trimmed_alignment}
			else
				if [[ "$(cat ALICUT_{params.busco}_aligned.fas | grep -v ">" | sed 's/-//g' | grep "^$" | wc -l)" -gt 0 ]]; then
					echo "$(date) - Alignment of BUSCO: {params.busco} contains empty sequence after aliscore/alicut. This sequence will be removed." >> {params.wd}/results/report.txt
					cp ALICUT_{params.busco}_aligned.fas ALICUT_{params.busco}_aligned.fas_tmp
					cat ALICUT_{params.busco}_aligned.fas_tmp | perl -ne 'chomp; $h=$_; $s=<>; chomp($s); $check=$s; $check=~s/-//g; if (length($check) > 0){{print "$h\n$s\n"}}' > ALICUT_{params.busco}_aligned.fas
				fi
				mv ALICUT_{params.busco}_aligned.fas {params.wd}/{output.trimmed_alignment}
			fi
			touch {params.wd}/{output.checkpoint}
			"""
			
#def aggregate_trimmed_alignments(wildcards):
#    checkpoint_output = checkpoints.create_sequence_files.get(**wildcards).output[0]
#    file_names = expand("results/trimmed_alignments/{busco}_aligned_trimmed.fas", busco = glob_wildcards(os.path.join(checkpoint_output, "{busco}_all.fas")).busco)
#    return file_names

rule get_all_trimmed_files:
        input:
                expand("results/trimmed_alignments/{bus}_aligned_trimmed.fas", bus=BUSCOS)
        output:
                checkpoint = "results/checkpoints/get_all_trimmed_files.done"
	singularity:
		"docker://continuumio/miniconda3:4.7.10"
	conda:
		"../envs/biopython.yml"
	params:
		wd = os.getcwd()
	shell:
                """
		if [[ -d results/filtered_alignments ]]; then
			rm -rf results/filtered_alignments
		fi
		mkdir results/filtered_alignments
		
		for file in results/trimmed_alignments/*.fas;
		do
			if [[ -s {params.wd}/$file ]]; then 
				python bin/filter_alignments.py --alignments {params.wd}/$file --outdir "{params.wd}/results/filtered_alignments"
			else #do nothing if file is empty (happens rarely when ALICUT fails)
				continue
			fi
		done
		echo "$(date) - Number of alignments: $(ls results/alignments/*.fas | wc -l)" >> results/report.txt
		echo "$(date) - Number of trimmed alignments: $(ls results/filtered_alignments/*.fas | wc -l)" >> results/report.txt
		echo "$(date) - Number of alignments after filtering: $(ls results/filtered_alignments/*.fas | wc -l)" >> results/report.txt
		touch {output.checkpoint}
		"""

