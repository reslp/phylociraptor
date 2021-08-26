configfile: "data/config.yaml"

import subprocess
import glob

BUSCOS, = glob_wildcards("results/orthology/busco/busco_sequences_deduplicated/{busco}_all.fas")

rule trim_trimal:
		input:
			checkpoint = "checkpoints/align.done",
			alignment = "results/alignments/full/{aligner}/{busco}_aligned.fas"
		output:
			trimmed_alignment = "results/alignments/trimmed/{aligner}-trimal/{busco}_aligned_trimmed.fas",
		benchmark:
			"results/statistics/benchmarks/align/{aligner}_trimal_{busco}.txt"
		params:
			trimmer = config["trimming"]["trimal_parameters"]
		singularity:
			"docker://reslp/trimal:1.4.1"
		shell:
			"""
			trimal {params.trimmer} -in {input.alignment} -out {output.trimmed_alignment}
			"""
rule trim_aliscotri:
		input:
			checkpoint = "checkpoints/align.done",
			alignment = "results/alignments/full/{aligner}/{busco}_aligned.fas"
		output:
			trimmed_alignment = "results/alignments/trimmed/{aligner}-aliscore/{busco}_aligned_trimmed.fas",
		benchmark:
			"results/statistics/benchmarks/align/{aligner}_aliscore_{busco}.txt"
		params:
			trimmer = config["trimming"]["aliscore_parameters"],
			busco = "{busco}",
			wd = os.getcwd()
		singularity:
			"docker://chrishah/alicut-aliscore-docker:2.31"
		shell:
			"""
			mkdir -p results/alignments/trimmed/{wildcards.aligner}-aliscore/{params.busco}
			cd "results/alignments/trimmed/{wildcards.aligner}-aliscore/{params.busco}"
			ln -s -f  {params.wd}/{input.alignment} {params.busco}_aligned.fas 
			Aliscore.pl {params.trimmer} -i {params.busco}_aligned.fas &> aliscore_{params.busco}.log || true
			
			if [[ -f {params.busco}_aligned.fas_List_random.txt ]]; then
				echo "$(date) - The aliscore output file does not exist. Check results for BUSCO: {params.busco}" >> {params.wd}/results/statistics/runlog.txt
				if [[ $(cat {params.busco}_aligned.fas_List_random.txt | head -n 1 | grep '[0-9]' -c) != 0 ]]; then
					echo "$(date) - The aliscore output appears to be empty. Check results for BUSCO: {params.busco}" >> {params.wd}/results/statistics/runlog.txt
					ALICUT.pl -s &> alicut_{params.busco}.log
				fi
			fi
			
			# this check is included because alicut very rarely does not  produce an output file.
			# in this case an empty file will be touched. This is necessary so the rule does not fail
			# The empty file will later be exluded again in the next rule.
			if [[ ! -f {params.wd}/results/alignments/trimmed/{wildcards.aligner}-aliscore/{params.busco}/ALICUT_{params.busco}_aligned.fas ]]; then
				echo "$(date) - The ALICUT output appears to be empty. Will touch an empty file so the pipeline will continue. Check results for BUSCO: {params.busco}" >> {params.wd}/results/statistics/runlog.txt
				touch {params.wd}/{output.trimmed_alignment}
			else
				if [[ "$(cat ALICUT_{params.busco}_aligned.fas | grep -v ">" | sed 's/-//g' | grep "^$" | wc -l)" -gt 0 ]]; then
					echo "$(date) - Alignment of BUSCO: {params.busco} contains empty sequence after aliscore/alicut. This sequence will be removed." >> {params.wd}/results/statistics/runlog.txt
					cp ALICUT_{params.busco}_aligned.fas ALICUT_{params.busco}_aligned.fas_tmp
					cat ALICUT_{params.busco}_aligned.fas_tmp | perl -ne 'chomp; $h=$_; $s=<>; chomp($s); $check=$s; $check=~s/-//g; if (length($check) > 0){{print "$h\n$s\n"}}' > ALICUT_{params.busco}_aligned.fas
				fi
				mv ALICUT_{params.busco}_aligned.fas {params.wd}/{output.trimmed_alignment}
			fi
			"""

rule get_trimmed_statistics:
	input:
		alignments = expand("results/alignments/trimmed/{{aligner}}-{{alitrim}}/{bus}_aligned_trimmed.fas", bus=BUSCOS)
	output:
		statistics_trimmed = "results/statistics/statistics_trimmed_{alitrim}_{aligner}-{batch}-"+str(config["filtering"]["concurrency"])+".txt",
		checkpoint = "results/checkpoints/get_trim_statistics_{alitrim}_{aligner}-{batch}-"+str(config["filtering"]["concurrency"])+".done"
	params:
		wd = os.getcwd(),
		ids = config["species"],
		datatype = config["filtering"]["seq_type"],
		nbatches = config["filtering"]["concurrency"],
	singularity: "docker://reslp/concat:0.21"
	shadow: "minimal"
	shell:
		"""
		# here the ids for the alignments need to be filtered as well first. maybe this can be changed in the concat.py script, so that an id file is not needed anymore.
		concat.py -i $(ls -1 {params.wd}/results/alignments/trimmed/{wildcards.aligner}-{wildcards.alitrim}/*.fas | sed -n '{wildcards.batch}~{params.nbatches}p' | tr '\\n' ' ') -t <(ls -1 {params.wd}/results/orthology/busco/busco_runs) --runmode concat -o results/statistics/ --biopython --statistics --seqtype {params.datatype} --noseq
		mv results/statistics/statistics.txt {output.statistics_trimmed}
		touch {output.checkpoint}
		"""

rule filter_alignments:
	input:
		expand("results/checkpoints/get_trim_statistics_{{alitrim}}_{{aligner}}-{batch}-"+str(config["filtering"]["concurrency"])+".done", aligner=config["alignment"]["method"], alitrim=config["alignment"]["method"], batch=range(1 , config["filtering"]["concurrency"]+1))
	output:
		checkpoint = "results/checkpoints/filter_alignments_{alitrim}_{aligner}.done",
		filter_info = "results/statistics/alignment_filter_information_{alitrim}_{aligner}.txt"
	benchmark:
		"results/statistics/benchmarks/align/filter_alignments_{alitrim}_{aligner}.txt"
	singularity:
		"docker://reslp/biopython_plus:1.77"
	params:
		wd = os.getcwd(),
		minsp=config["filtering"]["minsp"],
		trimming_method = config["trimming"]["method"],
		min_pars_sites = config["filtering"]["min_parsimony_sites"],
	shell:
		"""
		if [[ -d results/alignments/filtered/{wildcards.aligner}-{wildcards.alitrim} ]]; then
			rm -rf results/alignments/filtered/{wildcards.aligner}-{wildcards.alitrim}
		fi
		mkdir -p results/alignments/filtered/{wildcards.aligner}-{wildcards.alitrim}

		# concatenate the statistics files from the individual batches (for some reason snakemake complained if I did it all in one step, so this looks a bit ugly now, but it runs)
		cat results/statistics/statistics_trimmed_{wildcards.alitrim}_{wildcards.aligner}-* > results/statistics/statistics_trimmed_temp-{wildcards.alitrim}_{wildcards.aligner}.txt
		head -n 1 results/statistics/statistics_trimmed_temp-{wildcards.alitrim}_{wildcards.aligner}.txt > results/statistics/statistics_trimmed_{wildcards.alitrim}_{wildcards.aligner}.txt
		cat results/statistics/statistics_trimmed_temp-{wildcards.alitrim}_{wildcards.aligner}.txt | grep -v "^alignment" >> results/statistics/statistics_trimmed_{wildcards.alitrim}_{wildcards.aligner}.txt
		rm results/statistics/statistics_trimmed_temp-{wildcards.alitrim}_{wildcards.aligner}.txt

		for file in results/alignments/trimmed/{wildcards.aligner}-{wildcards.alitrim}/*.fas;
		do
			if [[ "$(cat {params.wd}/$file | grep ">" -c)" -lt {params.minsp} ]]; then
                                echo "$(date) - File $file contains less than {params.minsp} sequences after trimming with {params.trimming_method}. This file will not be used for tree reconstruction." >> {params.wd}/results/statistics/runlog.txt
                                        continue
                        fi
                        if [[ -s {params.wd}/$file ]]; then
                                python bin/filter_alignments.py --alignments {params.wd}/$file --outdir "{params.wd}/results/alignments/filtered/{wildcards.aligner}-{wildcards.alitrim}" --statistics-file results/statistics/statistics_trimmed_{wildcards.alitrim}_{wildcards.aligner}.txt --min-parsimony {params.min_pars_sites} --minsp {params.minsp} >> {output.filter_info}
                        else #do nothing if file is empty (happens rarely when ALICUT fails)
                                continue
                        fi
		done
		echo "$(date) - Number of alignments ({wildcards.aligner}): $(ls results/alignments/full/{wildcards.aligner}/*.fas | wc -l)" >> results/statistics/runlog.txt
		echo "$(date) - Number of trimmed alignments ({wildcards.aligner} - {wildcards.alitrim}): $(ls results/alignments/trimmed/{wildcards.aligner}-{wildcards.alitrim}/*.fas | wc -l)" >> results/statistics/runlog.txt
		echo "$(date) - Number of alignments ({wildcards.aligner} - {wildcards.alitrim}) after filtering: $(ls results/alignments/filtered/{wildcards.aligner}-{wildcards.alitrim}/*.fas | wc -l)" >> results/statistics/runlog.txt
		touch {output.checkpoint}
		"""

rule get_filter_statistics:
	input:
		rules.filter_alignments.output.checkpoint
	output:
		statistics_filtered = "results/statistics/statistics_filtered_{alitrim}_{aligner}-{batch}-"+str(config["filtering"]["concurrency"])+".txt",	
	params:
		wd = os.getcwd(),
		ids = config["species"],
		datatype = config["filtering"]["seq_type"],
		nbatches = config["filtering"]["concurrency"],
	singularity: "docker://reslp/concat:0.21"
	shadow: "minimal"
	shell:
		"""
		# here the ids for the alignments need to be filtered as well first. maybe this can be changed in the concat.py script, so that an id file is not needed anymore.
		concat.py -i $(ls -1 {params.wd}/results/alignments/filtered/{wildcards.aligner}-{wildcards.alitrim}/*.fas | sed -n '{wildcards.batch}~{params.nbatches}p' | tr '\\n' ' ') -t <(ls -1 {params.wd}/results/orthology/busco/busco_runs) --runmode concat -o results/statistics/ --biopython --statistics --seqtype {params.datatype} --noseq
		mv results/statistics/statistics.txt {output.statistics_filtered}
		"""

rule all_filter_align:
	input:
		#"results/checkpoints/get_all_trimmed_alignments.done",
		expand("results/checkpoints/{aligner}_aggregate_align.done", aligner=config["alignment"]["method"]),
		expand("results/checkpoints/filter_alignments_{alitrim}_{aligner}.done", aligner=config["alignment"]["method"], alitrim=config["trimming"]["method"]),
		expand("results/statistics/statistics_filtered_{alitrim}_{aligner}-{batch}-"+str(config["filtering"]["concurrency"])+".txt", aligner=config["alignment"]["method"], alitrim=config["trimming"]["method"], batch=range(1 , config["filtering"]["concurrency"]+1))
	output:
		"checkpoints/filter_align.done"
	shell:
		"""
		echo "$(date) - Pipeline part filter_align done." >> results/statistics/runlog.txt
		touch {output}
		"""	
