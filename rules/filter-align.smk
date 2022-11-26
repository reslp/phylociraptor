include: "functions.smk"
import subprocess
import glob
import yaml

ruleorder: read_params_global > read_params_per 

# get list of containers to use:
with open("data/containers.yaml", "r") as yaml_stream:
    containers = yaml.safe_load(yaml_stream)

#create new hashes for current stage (alignment) by combining the previuos hash with a newly generated
hashes = collect_hashes("filter-align")

filter_orthology_hash = hashes['filter-orthology']
aligner_hashes = hashes['align']["per"]
previous_hash = hashes['align']["global"]
current_hash = hashes['filter-align']["global"]
trimmer_hashes = hashes['filter-align']["per"]

BUSCOS, = glob_wildcards("results/orthology/busco/busco_sequences_deduplicated."+filter_orthology_hash+"/{busco}_all.fas")

def previous_params_global(wildcards):
	return "results/alignments/full/parameters.align."+previous_hash+".yaml"
def previous_params_per(wildcards):
	return "results/alignments/full/"+wildcards.aligner+"."+aligner_hashes[wildcards.aligner]+"/parameters.align."+wildcards.aligner+"."+aligner_hashes[wildcards.aligner]+".yaml"

def compare_filter_align(wildcards):
	return [trigger("results/alignments/trimmed/{aligner}-{alitrim}.{hash}/parameters.filter-align.{aligner}-{alitrim}.{hash}.yaml".format(aligner=wildcards.aligner, alitrim=wildcards.alitrim, hash=wildcards.hash), configfi)]	

def per_aligner(wildcards):
#	print(wildcards)
	return "results/alignments/full/"+wildcards.aligner+"."+aligner_hashes[wildcards.aligner]+"/"+wildcards.busco+"_aligned.fas"

def per_trimmer(wildcards):
	lis = []
#	print(wildcards)
	for b in BUSCOS:
		lis.append("results/alignments/trimmed/"+wildcards.aligner+"-"+wildcards.alitrim+"."+trimmer_hashes[wildcards.alitrim][wildcards.aligner]+"/"+str(b)+"_aligned_trimmed.fas")
	return lis

def determine_concurrency_limit():
	fname = "results/alignments/full"
	if os.path.isdir(fname):
		ngenes = glob.glob(fname+"/*/*.fas")
		ngenes = len(ngenes)
		if ngenes < config["concurrency"]:
			return range(1, ngenes + 1)
		else:
			return range(1, config["concurrency"] + 1)
	else:
		return

batches = determine_concurrency_limit()


aligners = get_aligners()		
trimmers = get_trimmers()		

rule read_params_per:
	input:
		trigger = compare_filter_align,
		previous = previous_params_per
	output:
		"results/alignments/trimmed/{aligner}-{alitrim}.{hash}/parameters.filter-align.{aligner}-{alitrim}.{hash}.yaml"
	shell:
		"""
		bin/read_write_yaml.py {input.trigger} {output} trimming,options,{wildcards.alitrim}
		cat {input.previous} >> {output}
		"""

rule read_params_global:
	input:
		trigger = compare("results/alignments/trimmed/parameters.filter-align."+current_hash+".yaml", configfi),
		previous = previous_params_global
	output:
		"results/alignments/trimmed/parameters.filter-align."+current_hash+".yaml"
	shell:
		"""
		bin/read_write_yaml.py {input.trigger} {output} trimming,method trimming,options trimming,min_parsimony_sites
		cat {input.previous} >> {output}
		"""
def return_aligner_checkpoint(wildcards):
	return "results/checkpoints/"+wildcards.aligner+"_aggregate_align."+aligner_hashes[wildcards.aligner]+".done"

def return_bmge_params(wildcards):
	return "results/alignments/trimmed/"+wildcards.aligner+"-bmge."+trimmer_hashes["bmge"][wildcards.aligner]+"/parameters.filter-align."+wildcards.aligner+"-bmge."+trimmer_hashes["bmge"][wildcards.aligner]+".yaml"
def return_trimal_params(wildcards):
	return "results/alignments/trimmed/"+wildcards.aligner+"-trimal."+trimmer_hashes["trimal"][wildcards.aligner]+"/parameters.filter-align."+wildcards.aligner+"-trimal."+trimmer_hashes["trimal"][wildcards.aligner]+".yaml"
def return_aliscore_params(wildcards):
	return "results/alignments/trimmed/"+wildcards.aligner+"-aliscore."+trimmer_hashes["aliscore"][wildcards.aligner]+"/parameters.filter-align."+wildcards.aligner+"-aliscore."+trimmer_hashes["aliscore"][wildcards.aligner]+".yaml"

rule bmge:
		input:
			checkpoint = return_aligner_checkpoint,
			params = return_bmge_params,
			alignment = per_aligner
		output:
			trimmed_alignment = "results/alignments/trimmed/{aligner}-bmge.{hash}/{busco}_aligned_trimmed.fas",
		benchmark:
			"results/statistics/benchmarks/align/{aligner}_bmge_{busco}.{hash}.txt"
		params:
			trimmer = config["trimming"]["options"]["bmge"],
			type = config["filtering"]["seq_type"] 
		log:
			"log/filter-align/bmge/bmge_{aligner}_{busco}.{hash}.txt"
		singularity:
			containers["bmge"]	
		shell:
			"""
			if [[ "{params.type}" == "aa" ]]
			then
				echo "Setting sequence type to AA" 2>&1 | tee {log}
				seqtype="AA"
			elif [[ "{params.type}" == "nu" ]]
			then
				seqtype="DNA"
				echo "Setting sequence type to DNA" 2>&1 | tee {log}
			fi
			bmge -i {input.alignment} -t $seqtype -of {output.trimmed_alignment} 2>&1 | tee -a {log}
			"""

rule trimal:
		input:
			checkpoint = return_aligner_checkpoint,
			params = return_trimal_params,
			alignment = per_aligner
		output:
			trimmed_alignment = "results/alignments/trimmed/{aligner}-trimal.{hash}/{busco}_aligned_trimmed.fas",
		benchmark:
			"results/statistics/benchmarks/align/{aligner}_trimal_{busco}.{hash}.txt"
		params:
			trimmer = config["trimming"]["options"]["trimal"]
		log:
			"log/filter-align/trimal/trimal_{aligner}_{busco}.{hash}.txt"
		singularity:
			containers["trimal"]	
		shell:
			"""
			trimal {params.trimmer} -in {input.alignment} -out {output.trimmed_alignment} 2>&1 | tee {log}
			#it can happen that the alignment is empty after trimming. In this case trimal doesn't produce an error, but does not write the file so we do that manually
			if [[ ! -f "{output.trimmed_alignment}" ]]; then echo "Making dummy file" 2>&1 | tee {log}; touch {output.trimmed_alignment}; fi
			"""
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
				echo "$(date) - The aliscore output file does not exist. Check results for BUSCO: {params.busco}" >> {params.wd}/results/statistics/runlog.txt
				if [[ $(cat {params.busco}_aligned.fas_List_random.txt | head -n 1 | grep '[0-9]' -c) != 0 ]]; then
					echo "$(date) - The aliscore output appears to be empty. Check results for BUSCO: {params.busco}" >> {params.wd}/results/statistics/runlog.txt
					ALICUT.pl -s 2>&1 | tee -a {params.wd}/{log}
				fi
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

rule get_trimmed_statistics:
	input:
		alignments = per_trimmer
	output:
		statistics_trimmed = "results/statistics/trim-{aligner}-{alitrim}.{hash}/stats_trimmed_{alitrim}_{aligner}-{batch}-"+str(config["concurrency"])+".txt"
		#checkpoint = "results/checkpoints/trim/get_trim_statistics_{alitrim}_{aligner}-{batch}-"+str(config["concurrency"])+".done"
	params:
		wd = os.getcwd(),
		ids = config["species"],
		datatype = config["filtering"]["seq_type"],
		nbatches = config["concurrency"],
		set = config["orthology"]["busco_options"]["set"],
		orthology_hash = hashes['orthology']
	singularity: containers["concat"] 
	shadow: "minimal"
	shell:
		"""
		# here the ids for the alignments need to be filtered as well first. maybe this can be changed in the concat.py script, so that an id file is not needed anymore.
		concat.py -i $(ls -1 {params.wd}/results/alignments/trimmed/{wildcards.aligner}-{wildcards.alitrim}.{wildcards.hash}/*.fas | sed -n '{wildcards.batch}~{params.nbatches}p' | tr '\\n' ' ') -t <(for name in $(ls -1 {params.wd}/results/orthology/busco/busco_runs.{params.set}.{params.orthology_hash}); do echo "${{name%.*}}"; done) --runmode concat -o results/statistics/ --biopython --statistics --seqtype {params.datatype} --noseq
		mv results/statistics/statistics.txt {output.statistics_trimmed}
		"""

def pull_trimmer(wildcards):
	lis = []
	for i in range(1, config["concurrency"] + 1):
		lis.append("results/statistics/trim-"+wildcards.aligner+"-"+wildcards.alitrim+"."+trimmer_hashes[wildcards.alitrim][wildcards.aligner]+"/stats_trimmed_"+wildcards.alitrim+"_"+wildcards.aligner+"-"+str(i)+"-"+str(config["concurrency"])+".txt")
	return lis

def get_aligner_hash(wildcards):
	return aligner_hashes[wildcards.aligner]

rule filter_alignments:
	input:
#		expand("results/statistics/trim-{{aligner}}-{{alitrim}}/stats_trimmed_{{alitrim}}_{{aligner}}-{batch}-"+str(config["concurrency"])+".txt", aligner=config["alignment"]["method"], alitrim=config["trimming"]["method"], batch=batches)
		pull_trimmer
	output:
		filter_info = "results/statistics/filter-{aligner}-{alitrim}.{hash}/alignment_filter_information_{alitrim}_{aligner}.txt",
		trim_info = "results/statistics/trim-{aligner}-{alitrim}.{hash}/statistics_trimmed_{alitrim}_{aligner}.txt"
	benchmark:
		"results/statistics/benchmarks/align/filter_alignments_{alitrim}_{aligner}.{hash}.txt"
	singularity:
		containers["biopython"]	
	params:
		wd = os.getcwd(),
		minsp=config["filtering"]["minsp"],
		trimming_method = config["trimming"]["method"],
		min_pars_sites = config["trimming"]["min_parsimony_sites"],
		aligner_hash = get_aligner_hash,
		target_dir = "results/alignments/filtered/{aligner}-{alitrim}.{hash}"
	shell:
		"""
		if [[ -d {params.target_dir} ]]; then
			rm -rf {params.target_dir}
		fi
		mkdir -p {params.target_dir}

		# concatenate the statistics files from the individual batches (for some reason snakemake complained if I did it all in one step, so this looks a bit ugly now, but it runs)
		cat results/statistics/trim-{wildcards.aligner}-{wildcards.alitrim}.{wildcards.hash}/stats_trimmed_{wildcards.alitrim}_{wildcards.aligner}-* > results/statistics/trim-{wildcards.aligner}-{wildcards.alitrim}.{wildcards.hash}/statistics_trimmed_temp-{wildcards.alitrim}_{wildcards.aligner}.txt
		head -n 1 results/statistics/trim-{wildcards.aligner}-{wildcards.alitrim}.{wildcards.hash}/statistics_trimmed_temp-{wildcards.alitrim}_{wildcards.aligner}.txt > results/statistics/trim-{wildcards.aligner}-{wildcards.alitrim}.{wildcards.hash}/statistics_trimmed_{wildcards.alitrim}_{wildcards.aligner}.txt
		

		echo "\\n ######## BEFORE #######"
		grep -v "^alignment" results/statistics/trim-{wildcards.aligner}-{wildcards.alitrim}.{wildcards.hash}/statistics_trimmed_temp-{wildcards.alitrim}_{wildcards.aligner}.txt >> {output.trim_info}
		echo "\\n ######## AFTER #######"
		rm results/statistics/trim-{wildcards.aligner}-{wildcards.alitrim}.{wildcards.hash}/statistics_trimmed_temp-{wildcards.alitrim}_{wildcards.aligner}.txt
		for file in results/alignments/trimmed/{wildcards.aligner}-{wildcards.alitrim}.{wildcards.hash}/*.fas;
		do
                        if [[ -s {params.wd}/$file ]]; then
				if [[ "$(cat {params.wd}/$file | grep ">" -c)" -lt {params.minsp} ]]; then
					echo -e "$file\tTOO_FEW_SEQUENCES" 2>&1 | tee -a {output.filter_info}
                                	echo "$(date) - File $file contains less than {params.minsp} sequences after trimming with {params.trimming_method}. This file will not be used for tree reconstruction." >> {params.wd}/results/statistics/runlog.txt
                                        	continue
                        	fi
				if [[ $(grep "$(basename $file)" {output.trim_info} | cut -f 6) -ge {params.min_pars_sites} ]]
				then
					echo -e "$file\tPASS" 2>&1 | tee -a {output.filter_info}
					ln -s {params.wd}/$file {params.wd}/{params.target_dir}/
				else
					echo -e "$file\tNOT_INFORMATIVE" 2>&1 | tee -a {output.filter_info}
				fi

#                                python bin/filter_alignments.py --alignments {params.wd}/$file --outdir "{params.wd}/{params.target_dir}" --statistics-file results/statistics/trim-{wildcards.aligner}-{wildcards.alitrim}.{wildcards.hash}/statistics_trimmed_{wildcards.alitrim}_{wildcards.aligner}.txt --min-parsimony {params.min_pars_sites} --minsp {params.minsp} >> {output.filter_info}
                        else #do nothing if file is empty (happens rarely when ALICUT fails)
                                continue
                        fi
		done
		# remove unnecessary statistics file
#		rm results/statistics/trim-{wildcards.aligner}-{wildcards.alitrim}.{wildcards.hash}/statistics_trimmed_{wildcards.alitrim}_{wildcards.aligner}.txt

		echo "$(date) - Number of alignments ({wildcards.aligner}): $(ls results/alignments/full/{wildcards.aligner}.{params.aligner_hash}/*.fas | wc -l)" >> results/statistics/runlog.txt
		echo "$(date) - Number of trimmed alignments ({wildcards.aligner} - {wildcards.alitrim}): $(ls results/alignments/trimmed/{wildcards.aligner}-{wildcards.alitrim}.{wildcards.hash}/*.fas | wc -l)" >> results/statistics/runlog.txt
		echo "$(date) - Number of alignments ({wildcards.aligner} - {wildcards.alitrim}) after filtering: $(ls {params.target_dir}/*.fas | wc -l)" >> results/statistics/runlog.txt
		"""

rule get_filter_statistics:
	input:
		rules.filter_alignments.output.filter_info
	output:
		statistics_filtered = "results/statistics/filter-{aligner}-{alitrim}.{hash}/statistics_filtered_{alitrim}_{aligner}-{batch}-"+str(config["concurrency"])+".txt",	
	params:
		wd = os.getcwd(),
		ids = config["species"],
		datatype = config["filtering"]["seq_type"],
		nbatches = config["concurrency"],
		set = config["orthology"]["busco_options"]["set"],
		orthology_hash = hashes['orthology']
	singularity: containers["concat"] 
	shadow: "minimal"
	shell:
		"""
		# here the ids for the alignments need to be filtered as well first. maybe this can be changed in the concat.py script, so that an id file is not needed anymore.
		concat.py -i $(ls -1 {params.wd}/results/alignments/filtered/{wildcards.aligner}-{wildcards.alitrim}.{wildcards.hash}/*.fas | sed -n '{wildcards.batch}~{params.nbatches}p' | tr '\\n' ' ') -t <(for name in $(ls -1 {params.wd}/results/orthology/busco/busco_runs.{params.set}.{params.orthology_hash}); do echo "${{name%.*}}"; done) --runmode concat -o results/statistics/ --biopython --statistics --seqtype {params.datatype} --noseq
		mv results/statistics/statistics.txt {output.statistics_filtered}
		"""

def pull_filtered_stats(wildcards):
	lis = []
	for t in trimmer_hashes.keys():
		for a in trimmer_hashes[t].keys():
			for i in range(1, config["concurrency"] + 1):
				lis.append("results/statistics/filter-"+a+"-"+t+"."+trimmer_hashes[t][a]+"/statistics_filtered_"+t+"_"+a+"-"+str(i)+"-"+str(config["concurrency"])+".txt")
	return lis

def pull_algn_info(wildcards):
	lis = []
	for t in trimmer_hashes.keys():
		for a in trimmer_hashes[t].keys():
			for i in range(1, config["concurrency"] + 1):
				lis.append("results/statistics/filter-"+a+"-"+t+"."+trimmer_hashes[t][a]+"/alignment_filter_information_"+t+"_"+a+".txt")
	return lis

rule filter_align:
	input:
		filtered = pull_filtered_stats,
		algntrim = pull_algn_info,
		params = rules.read_params_global.output,
	output:
		checkpoint = "results/checkpoints/modes/filter_align."+current_hash+".done",
		stats = "results/statistics/filtered-alignments-stats."+current_hash+".txt"
	params:
		dupseq = config["filtering"]["dupseq"],
		cutoff = config["filtering"]["cutoff"],
		minsp = config["filtering"]["minsp"],
		seq_type = config["filtering"]["seq_type"],
		min_parsimony_sites = config["trimming"]["min_parsimony_sites"],
		algn_trim = expand("{al}__{tr}", al=config["alignment"]["method"], tr=config["trimming"]["method"]),
		trimal_params = config["trimming"]["options"]["trimal"],
		aliscore_params = config["trimming"]["options"]["aliscore"],
		bmge_params = config["trimming"]["options"]["bmge"]
	shell:
		"""
		echo "combo\t$(head -n 1 {input.filtered[0]} | cut -f 1,4-)" > {output.stats}
		for f in {input.filtered}; do combo=$(echo $f | cut -d "/" -f 3 | sed 's/filter-//'); tail -n +2 $f | sed "s/^/$combo\t/"; done | cut -f 1,2,4- >> {output.stats}

		# gather aligner and trimming combinations and corresponding settings for statistics:	
		combs="{params.algn_trim}"
		echo "aligner\ttrimmer\ttrimming_params\tdupseq\tcutoff\tminsp\tseq_type\tmin_parsimony_sites" > results/statistics/trim-filter-overview.txt
		for combination in ${{combs}}; do
			if [[ $combination == *"trimal"* ]]; then
				out=$combination"__{params.trimal_params}"
			elif [[ $combination == *"aliscore"* ]]; then
				out=$combination"__{params.aliscore_params}"
			elif [[ $combination == *"bmge"* ]]; then
				out=$combination"__{params.bmge_params}"
			else
				out=$combination"__no parameters"
			fi
			out=$out"__{params.dupseq}__{params.cutoff}__{params.minsp}__{params.seq_type}__{params.min_parsimony_sites}"
			echo $out | sed 's/__/\t/g' >> results/statistics/trim-filter-overview.txt
		done  
		echo "$(date) - phylociraptor filter-align done." >> results/statistics/runlog.txt
		touch {output.checkpoint}
		"""	
