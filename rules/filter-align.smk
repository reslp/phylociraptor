configfile: "data/config.yaml"

import subprocess
import glob

BUSCOS, = glob_wildcards("results/alignments/{busco}_aligned.fas")

#def get_all_trimmed_alignments(wildcards):
#	return glob.glob("results/trimmed_alignments/*.fas")

if config["trimming"]["method"] == "trimal":
	rule trim:
		input:
			"results/alignments/{busco}_aligned.fas"
		output:
			trimmed_alignment = "results/trimmed_alignments/{busco}_aligned_trimmed.fas",
		benchmark:
			"results/statistics/benchmarks/align/trim_{busco}.txt"
		params:
			trimmer = config["trimming"]["parameters"]
		singularity:
			"docker://reslp/trimal:1.4.1"
		shell:
			"""
			trimal {params.trimmer} -in {input} -out {output.trimmed_alignment}
			"""
elif config["trimming"]["method"] == "aliscore":
	rule trim:
		input:
			"results/alignments/{busco}_aligned.fas"
		output:
			trimmed_alignment = "results/trimmed_alignments/{busco}_aligned_trimmed.fas",
		benchmark:
			"results/statistics/benchmarks/align/trim_{busco}.txt"
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
			if [[ ! -f {params.wd}/results/trimmed_alignments/{params.busco}/ALICUT_{params.busco}_aligned.fas ]]; then
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

#rule get_all_trimmed_alignments:
#	input:
#		#get_all_trimmed_alignments
#		expand("results/trimmed_alignments/{bus}_aligned_trimmed.fas", bus=BUSCOS)
#	output:
#		checkpoint = "results/checkpoints/get_all_trimmed_alignments.done"
#	benchmark:
#		"results/statistics/benchmarks/align/get_all_trimmed_alignments.txt"
#	singularity:
#		"docker://reslp/biopython_plus:1.77"
#	params:
#		wd = os.getcwd(),
#		trimming_method = config["trimming"]["method"]
#	shell:
#		"""
#		if [[ -d results/filtered_alignments ]]; then
#			rm -rf results/filtered_alignments
#		fi
#		mkdir results/filtered_alignments
#		
#		for file in results/trimmed_alignments/*.fas;
#		do
#			if [[ "$(cat {params.wd}/$file | grep ">" -c)" -lt 3 ]]; then
#				echo "$(date) - File $file contains less than 3 sequences after trimming with {params.trimming_method}. This file will not be used for tree reconstruction." >> {params.wd}/results/statistics/runlog.txt
#					continue
#			fi	
#			if [[ -s {params.wd}/$file ]]; then 
#				python bin/filter_alignments.py --alignments {params.wd}/$file --outdir "{params.wd}/results/filtered_alignments"
#			else #do nothing if file is empty (happens rarely when ALICUT fails)
#				continue
#			fi
#		done
#		echo "$(date) - Number of alignments: $(ls results/alignments/*.fas | wc -l)" >> results/statistics/runlog.txt
#		echo "$(date) - Number of trimmed alignments: $(ls results/filtered_alignments/*.fas | wc -l)" >> results/statistics/runlog.txt
#		echo "$(date) - Number of alignments after filtering: $(ls results/filtered_alignments/*.fas | wc -l)" >> results/statistics/runlog.txt
#		touch {output.checkpoint}
#		"""

rule get_trimmed_statistics:
	input:
		alignments = expand("results/trimmed_alignments/{bus}_aligned_trimmed.fas", bus=BUSCOS)
	output:
		statistics_trimmed = "results/statistics/statistics_trimmed.txt",
		checkpoint = "results/checkpoints/get_trim_statistics.done"
	params:
		ids = config["species"],
		datatype = config["filtering"]["seq_type"],
	singularity: "docker://reslp/concat:0.21"
	shell:
		"""
		tail -n +2 {params.ids} | awk -F "," '{{print $1;}}' | sed 's/^"//g' | sed 's/"$//g' | sed 's/ /_/g' > results/statistics/ids_alignments.txt
		# here the ids for the alignments need to be filtered as well first. maybe this can be changed in the concat.py script, so that an id file is not needed anymore.
		concat.py -d results/trimmed_alignments/ -t results/statistics/ids_alignments.txt --runmode concat -o results/statistics/ --biopython --statistics --seqtype aa --noseq
		mv results/statistics/statistics.txt {output.statistics_trimmed}
		touch {output.checkpoint}
		"""

rule filter_alignments:
	input:
		#alignments = expand("results/trimmed_alignments/{bus}_aligned_trimmed.fas", bus=BUSCOS),
		rules.get_trimmed_statistics.output.checkpoint
		#stats = rules.get_alignment_statistics.output.statistics_trimmed
	output:
		checkpoint = "results/checkpoints/filter_alignments.done",
		filter_info = "results/statistics/alignment_filter_information.txt"
	benchmark:
		"results/statistics/benchmarks/align/filter_alignments.txt"
	singularity:
		"docker://reslp/biopython_plus:1.77"
	params:
		wd = os.getcwd(),
		trimming_method = config["trimming"]["method"],
		min_pars_sites = config["filtering"]["min_parsimony_sites"],
		stats = "results/statistics/statistics_trimmed.txt"
	shell:
		"""
		if [[ -d results/filtered_alignments ]]; then
			rm -rf results/filtered_alignments
		fi
		mkdir results/filtered_alignments
		
		for file in results/trimmed_alignments/*.fas;
		do
			if [[ "$(cat {params.wd}/$file | grep ">" -c)" -lt 3 ]]; then
                                echo "$(date) - File $file contains less than 3 sequences after trimming with {params.trimming_method}. This file will not be used for tree reconstruction." >> {params.wd}/results/statistics/runlog.txt
                                        continue
                        fi
                        if [[ -s {params.wd}/$file ]]; then
                                python bin/filter_alignments.py --alignments {params.wd}/$file --outdir "{params.wd}/results/filtered_alignments" --statistics-file {params.stats} --min-parsimony {params.min_pars_sites} >> {output.filter_info}
                        else #do nothing if file is empty (happens rarely when ALICUT fails)
                                continue
                        fi
		done
		echo "$(date) - Number of alignments: $(ls results/alignments/*.fas | wc -l)" >> results/statistics/runlog.txt
		echo "$(date) - Number of trimmed alignments: $(ls results/filtered_alignments/*.fas | wc -l)" >> results/statistics/runlog.txt
		echo "$(date) - Number of alignments after filtering: $(ls results/filtered_alignments/*.fas | wc -l)" >> results/statistics/runlog.txt
		touch {output.checkpoint}
		"""

rule get_filter_statistics:
	input:
		rules.filter_alignments.output.checkpoint
	output:
		statistics_filtered = "results/statistics/statistics_filtered.txt",	
		statistics_trimmed = "results/statistics/statistics_trimmed.txt"
	params:
		ids = config["species"],
		datatype = config["filtering"]["seq_type"],
	singularity: "docker://reslp/concat:0.21"
	shell:
		"""
		tail -n +2 {params.ids} | awk -F "," '{{print $1;}}' | sed 's/^"//g' | sed 's/"$//g' | sed 's/ /_/g' > results/statistics/ids_alignments.txt
		# here the ids for the alignments need to be filtered as well first. maybe this can be changed in the concat.py script, so that an id file is not needed anymore.
		concat.py -d results/filtered_alignments/ -t results/statistics/ids_alignments.txt --runmode concat -o results/statistics/ --biopython --statistics --seqtype {params.datatype} --noseq
		mv results/statistics/statistics.txt {output.statistics_filtered}
		concat.py -d results/trimmed_alignments/ -t results/statistics/ids_alignments.txt --runmode concat -o results/statistics/ --biopython --statistics --seqtype aa --noseq
		mv results/statistics/statistics.txt {output.statistics_trimmed}
		"""

rule part_filter_align:
	input:
		#"results/checkpoints/get_all_trimmed_alignments.done",
		"results/checkpoints/aggregate_align.done",
		"results/checkpoints/filter_alignments.done",
		"results/statistics/statistics_filtered.txt"
	output:
		"checkpoints/filter_align.done"
	shell:
		"""
		echo "$(date) - Pipeline part filter_align done." >> results/statistics/runlog.txt
		touch {output}
		"""	
