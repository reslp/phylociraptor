rule aliscore:
		input:
			checkpoint = return_aligner_checkpoint,
			params = return_aliscore_params,
			alignment = per_aligner
		output:
			trimmed_alignment = "results/alignments/trimmed/{aligner}-aliscore.{hash}/{busco}_aligned_trimmed.fas",
		benchmark:
			"results/statistics/benchmarks/align/{aligner}_aliscore_{busco}.{hash}.txt"
		params:
			trimmer = config["trimming"]["options"]["aliscore"],
			busco = "{busco}",
			target_dir = "results/alignments/trimmed/{aligner}-aliscore.{hash}/{busco}",
			wd = os.getcwd()
		log:
			"log/filter-align/aliscore/aliscore_alicut_{aligner}_{busco}.{hash}.txt"
		singularity:
			containers["aliscore"]	
		shell:
			"""
			if [[ -d {params.target_dir} ]]; then rm -rf {params.target_dir}; fi
			mkdir -p {params.target_dir}
			cd {params.target_dir}
			ln -s -f  {params.wd}/{input.alignment} {params.busco}_aligned.fas 
			
			echo "ALISCORE output:\n" > {params.wd}/{log}	
			Aliscore.pl $(if [[ "{params.trimmer}" != "None" ]]; then echo "{params.trimmer}"; fi) -i {params.busco}_aligned.fas 2>&1 | tee {params.wd}/{log} aliscore_{params.busco}.log || true
			
			echo "\n\n ALICUT output:\n" >> {params.wd}/{log}	
			if [[ -f {params.busco}_aligned.fas_List_random.txt ]]; then
				if [[ $(cat {params.busco}_aligned.fas_List_random.txt | head -n 1 | grep '[0-9]' -c) != 0 ]]; then
					ALICUT.pl -s 2>&1 | tee -a {params.wd}/{log}
				else
					echo "$(date) - The aliscore output {params.wd}/{params.target_dir}/{params.busco}_aligned.fas_List_random.txt appears to be empty. Check results for BUSCO: {params.busco}" >> {params.wd}/results/statistics/runlog.txt
				fi
			else
				echo "$(date) - The aliscore output file {params.wd}/{params.target_dir}/{params.busco}_aligned.fas_List_random.txt does not exist. Check results for BUSCO: {params.busco}" >> {params.wd}/results/statistics/runlog.txt
				
			fi
			
			# this check is included because alicut very rarely does not  produce an output file.
			# in this case an empty file will be touched. This is necessary so the rule does not fail
			# The empty file will later be exluded again in the next rule.
			if [[ ! -f {params.wd}/{params.target_dir}/ALICUT_{params.busco}_aligned.fas ]]; then
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
