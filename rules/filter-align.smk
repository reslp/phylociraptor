include: "functions.smk"
import subprocess
import glob
import yaml


ruleorder: read_params_global > read_params_per 

# get list of containers to use:
with open("data/containers.yaml", "r") as yaml_stream:
    containers = yaml.safe_load(yaml_stream)

#create new hashes for current stage (alignment) by combining the previuos hash with a newly generated
hashes = collect_hashes("filter-align", config, configfi, wd=os.getcwd())

filter_orthology_hash = hashes['filter-orthology']["global"]
aligner_hashes = hashes['align']["per"]
previous_hash = hashes['align']["global"]
current_hash = hashes['filter-align']["global"]
trimmed_hashes = hashes['filter-align']["per-trimming"]
filtered_hashes = hashes['filter-align']["per"]

BUSCOS, = glob_wildcards("results/orthology/single-copy-orthologs_deduplicated."+filter_orthology_hash+"/{busco}_all.fas")

def previous_params_global(wildcards):
	return "results/alignments/full/parameters.align."+previous_hash+".yaml"

def previous_params_per(wildcards):
	return "results/alignments/full/"+wildcards.aligner+"."+aligner_hashes[wildcards.aligner]+"/parameters.align."+wildcards.aligner+"."+aligner_hashes[wildcards.aligner]+".yaml"

def compare_filter_align(wildcards):
	return [trigger("results/alignments/filtered/{aligner}-{alitrim}.{hash}/parameters.filter-align.{aligner}-{alitrim}.{hash}.yaml".format(aligner=wildcards.aligner, alitrim=wildcards.alitrim, hash=wildcards.hash), configfi)]	
def compare_trimmed_align(wildcards):
	return [trigger("results/alignments/trimmed/{aligner}-{alitrim}.{hash}/parameters.filter-align.{aligner}-{alitrim}.{hash}.yaml".format(aligner=wildcards.aligner, alitrim=wildcards.alitrim, hash=wildcards.trimhash), configfi)]	

def per_aligner(wildcards):
#	print(wildcards)
	return "results/alignments/full/"+wildcards.aligner+"."+aligner_hashes[wildcards.aligner]+"/"+wildcards.busco+"_aligned.fas"

def per_trimmer(wildcards):
	lis = []
#	print(wildcards)
	for b in BUSCOS:
		lis.append("results/alignments/trimmed/"+wildcards.aligner+"-"+wildcards.alitrim+"."+trimmed_hashes[wildcards.alitrim][wildcards.aligner]+"/"+str(b)+"_aligned_trimmed.fas")
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

rule read_params_per_trimmer:
	input:
		trigger = compare_trimmed_align,
		previous = previous_params_per
	output:
		"results/alignments/trimmed/{aligner}-{alitrim}.{trimmhash}/parameters.filter-align.{aligner}-{alitrim}.{trimhash}.yaml"
	params:
		configfile = configfi
	shell:
		"""
		bin/read_write_yaml.py {input.trigger} {output} trimming,options,{wildcards.alitrim}
		cat {input.previous} >> {output}
		"""


rule read_params_per:
	input:
		trigger = compare_filter_align,
		previous = previous_params_per
	output:
		"results/alignments/filtered/{aligner}-{alitrim}.{hash}/parameters.filter-align.{aligner}-{alitrim}.{hash}.yaml"
	shell:
		"""
		bin/read_write_yaml.py {input.trigger} {output} trimming,options,{wildcards.alitrim} trimming,min_parsimony_sites trimming,max_rcv_score
		cat {input.previous} >> {output}
		"""

rule read_params_global:
	input:
		trigger = compare("results/alignments/filtered/parameters.filter-align."+current_hash+".yaml", configfi),
		previous = previous_params_global
	output:
		"results/alignments/trimmed/parameters.filter-align."+current_hash+".yaml"
	shell:
		"""
		bin/read_write_yaml.py {input.trigger} {output} trimming,method trimming,options trimming,min_parsimony_sites trimming,max_rcv_score
		cat {input.previous} >> {output}
		"""
def return_aligner_checkpoint(wildcards):
	return "results/checkpoints/"+wildcards.aligner+"_aggregate_align."+aligner_hashes[wildcards.aligner]+".done"

def return_bmge_params(wildcards):
	return "results/alignments/trimmed/"+wildcards.aligner+"-bmge."+trimmed_hashes["bmge"][wildcards.aligner]+"/parameters.filter-align."+wildcards.aligner+"-bmge."+trimmed_hashes["bmge"][wildcards.aligner]+".yaml"

def return_trimal_params(wildcards):
	return "results/alignments/trimmed/"+wildcards.aligner+"-trimal."+trimmed_hashes["trimal"][wildcards.aligner]+"/parameters.filter-align."+wildcards.aligner+"-trimal."+trimmed_hashes["trimal"][wildcards.aligner]+".yaml"

def return_aliscore_params(wildcards):
	return "results/alignments/trimmed/"+wildcards.aligner+"-aliscore."+trimmed_hashes["aliscore"][wildcards.aligner]+"/parameters.filter-align."+wildcards.aligner+"-aliscore."+trimmed_hashes["aliscore"][wildcards.aligner]+".yaml"



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
		orthology_hash = hashes['orthology']["global"],
		mode = config["orthology"]["method"]
	singularity: containers["concat"] 
	shadow: "minimal"
	shell:
		"""
		# here the ids for the alignments need to be filtered as well first. maybe this can be changed in the concat.py script, so that an id file is not needed anymore.
		if [[ "{params.mode}" == "orthofinder" ]]; then
			concat.py -i $(ls -1 {params.wd}/results/alignments/trimmed/{wildcards.aligner}-{wildcards.alitrim}.{wildcards.hash}/*.fas | sed -n '{wildcards.batch}~{params.nbatches}p' | tr '\\n' ' ') --runmode concat -o results/statistics/ --biopython --statistics --seqtype {params.datatype} --noseq
		else
			concat.py -i $(ls -1 {params.wd}/results/alignments/trimmed/{wildcards.aligner}-{wildcards.alitrim}.{wildcards.hash}/*.fas | sed -n '{wildcards.batch}~{params.nbatches}p' | tr '\\n' ' ') -t <(for name in $(ls -1 {params.wd}/results/orthology/busco/busco_runs.{params.set}.{params.orthology_hash}); do echo "${{name%.*}}"; done) --runmode concat -o results/statistics/ --biopython --statistics --seqtype {params.datatype} --noseq
		fi

		mv results/statistics/statistics.txt {output.statistics_trimmed}
		"""

def pull_trimmer(wildcards):
	lis = []
	for i in range(1, config["concurrency"] + 1):
		lis.append("results/statistics/trim-"+wildcards.aligner+"-"+wildcards.alitrim+"."+trimmed_hashes[wildcards.alitrim][wildcards.aligner]+"/stats_trimmed_"+wildcards.alitrim+"_"+wildcards.aligner+"-"+str(i)+"-"+str(config["concurrency"])+".txt")
	return lis

def get_aligner_hash(wildcards):
	return aligner_hashes[wildcards.aligner]

def get_trimmer_hash(wildcards):
	return trimmed_hashes[wildcards.alitrim][wildcards.aligner]

rule filter_alignments:
	input:
		pull_trimmer,
		"results/alignments/filtered/{aligner}-{alitrim}.{hash}/parameters.filter-align.{aligner}-{alitrim}.{hash}.yaml"	
	output:
		filter_info = "results/statistics/filter-{aligner}-{alitrim}.{hash}/alignment_filter_information_{alitrim}_{aligner}.txt",
	benchmark:
		"results/statistics/benchmarks/align/filter_alignments_{alitrim}_{aligner}.{hash}.txt"
	singularity:
		containers["biopython"]	
	params:
		wd = os.getcwd(),
		minsp=config["trimming"]["minsp"],
		trimming_method = config["trimming"]["method"],
		min_pars_sites = config["trimming"]["min_parsimony_sites"],
		max_rcv_score = config["trimming"]["max_rcv_score"],
		aligner_hash = get_aligner_hash,
		trimmer_hash = get_trimmer_hash,
		target_dir = "results/alignments/filtered/{aligner}-{alitrim}.{hash}"
	shell:
		"""

		# concatenate the statistics files from the individual batches (for some reason snakemake complained if I did it all in one step, so this looks a bit ugly now, but it runs)
		cat results/statistics/trim-{wildcards.aligner}-{wildcards.alitrim}.{params.trimmer_hash}/stats_trimmed_{wildcards.alitrim}_{wildcards.aligner}-* > results/statistics/trim-{wildcards.aligner}-{wildcards.alitrim}.{params.trimmer_hash}/statistics_trimmed_temp-{wildcards.alitrim}_{wildcards.aligner}.txt
		head -n 1 results/statistics/trim-{wildcards.aligner}-{wildcards.alitrim}.{params.trimmer_hash}/statistics_trimmed_temp-{wildcards.alitrim}_{wildcards.aligner}.txt > results/statistics/trim-{wildcards.aligner}-{wildcards.alitrim}.{params.trimmer_hash}/statistics_trimmed_{wildcards.alitrim}_{wildcards.aligner}.txt
		

		echo "\\n ######## BEFORE #######"
		grep -v "^alignment" results/statistics/trim-{wildcards.aligner}-{wildcards.alitrim}.{params.trimmer_hash}/statistics_trimmed_temp-{wildcards.alitrim}_{wildcards.aligner}.txt >> results/statistics/trim-{wildcards.aligner}-{wildcards.alitrim}.{params.trimmer_hash}/statistics_trimmed_{wildcards.alitrim}_{wildcards.aligner}.txt
		echo "\\n ######## AFTER #######"
		rm results/statistics/trim-{wildcards.aligner}-{wildcards.alitrim}.{params.trimmer_hash}/statistics_trimmed_temp-{wildcards.alitrim}_{wildcards.aligner}.txt
		for file in results/alignments/trimmed/{wildcards.aligner}-{wildcards.alitrim}.{params.trimmer_hash}/*.fas;
		do
                        if [[ -s {params.wd}/$file ]]; then
				if [[ "$(cat {params.wd}/$file | grep ">" -c)" -lt {params.minsp} ]]; then
					echo -e "$file\tTOO_FEW_SEQUENCES" 2>&1 | tee -a {output.filter_info}
                                	echo "$(date) - File $file contains less than {params.minsp} sequences after trimming with {params.trimming_method}. This file will not be used for tree reconstruction." >> {params.wd}/results/statistics/runlog.txt
                                        continue
				fi
				if [[ $(grep "$(basename $file)" results/statistics/trim-{wildcards.aligner}-{wildcards.alitrim}.{params.trimmer_hash}/statistics_trimmed_{wildcards.alitrim}_{wildcards.aligner}.txt | cut -f 6) -ge {params.min_pars_sites} && $(grep "$(basename $file)" results/statistics/trim-{wildcards.aligner}-{wildcards.alitrim}.{params.trimmer_hash}/statistics_trimmed_{wildcards.alitrim}_{wildcards.aligner}.txt | cut -f 9 | awk '{{ if ($1 <= {params.max_rcv_score}) {{print "OK"}} }}') == "OK" ]]; then
					echo -e "$file\tPASS" 2>&1 | tee -a {output.filter_info}
					ln -s -f ../../../../$file {params.wd}/{params.target_dir}/
					continue
				elif [[ $(grep "$(basename $file)" results/statistics/trim-{wildcards.aligner}-{wildcards.alitrim}.{params.trimmer_hash}/statistics_trimmed_{wildcards.alitrim}_{wildcards.aligner}.txt | cut -f 6) -lt {params.min_pars_sites} ]]; then # case if alignment has too few pars inf sites
					echo -e "$file\tNOT_INFORMATIVE" 2>&1 | tee -a {output.filter_info}
					continue
				elif [[ $(grep "$(basename $file)" results/statistics/trim-{wildcards.aligner}-{wildcards.alitrim}.{params.trimmer_hash}/statistics_trimmed_{wildcards.alitrim}_{wildcards.aligner}.txt | cut -f 9 | awk '{{ if ($1 <= {params.max_rcv_score}) {{print "OK"}} }}') != "OK" ]]; then # case if rcv is too high
					echo -e "$file\tRCV_TOO_HIGH" 2>&1 | tee -a {output.filter_info}
					continue
				else # should not happen!
					echo -e "$file\tINVALID" 2>&1 | tee -a {output.filter_info}
					continue
				fi
				
                        else #do nothing if file is empty (happens rarely when ALICUT fails)
                                continue
                        fi
		done
		# remove unnecessary statistics file
#		rm results/statistics/trim-{wildcards.aligner}-{wildcards.alitrim}.{wildcards.hash}/statistics_trimmed_{wildcards.alitrim}_{wildcards.aligner}.txt

		echo "$(date) - Number of alignments ({wildcards.aligner}): $(ls results/alignments/full/{wildcards.aligner}.{params.aligner_hash}/*.fas | wc -l)" >> results/statistics/runlog.txt
		echo "$(date) - Number of trimmed alignments ({wildcards.aligner} - {wildcards.alitrim}): $(ls results/alignments/trimmed/{wildcards.aligner}-{wildcards.alitrim}.{params.trimmer_hash}/*.fas | wc -l)" >> results/statistics/runlog.txt
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
		orthology_hash = hashes['orthology']["global"],
		mode = config["orthology"]["method"]
	singularity: containers["concat"] 
	shadow: "minimal"
	shell:
		"""
		# here the ids for the alignments need to be filtered as well first. maybe this can be changed in the concat.py script, so that an id file is not needed anymore.
		if [[ "{params.mode}" == "orthofinder" ]]; then
			concat.py -i $(ls -1 {params.wd}/results/alignments/filtered/{wildcards.aligner}-{wildcards.alitrim}.{wildcards.hash}/*.fas | sed -n '{wildcards.batch}~{params.nbatches}p' | tr '\\n' ' ') --runmode concat -o results/statistics/ --biopython --statistics --seqtype {params.datatype} --noseq
		else
			concat.py -i $(ls -1 {params.wd}/results/alignments/filtered/{wildcards.aligner}-{wildcards.alitrim}.{wildcards.hash}/*.fas | sed -n '{wildcards.batch}~{params.nbatches}p' | tr '\\n' ' ') -t <(for name in $(ls -1 {params.wd}/results/orthology/busco/busco_runs.{params.set}.{params.orthology_hash}); do echo "${{name%.*}}"; done) --runmode concat -o results/statistics/ --biopython --statistics --seqtype {params.datatype} --noseq
		fi
		mv results/statistics/statistics.txt {output.statistics_filtered}
		"""

def pull_filtered_stats(wildcards):
	lis = []
	for t in filtered_hashes.keys():
		for a in filtered_hashes[t].keys():
			for i in range(1, config["concurrency"] + 1):
				lis.append("results/statistics/filter-"+a+"-"+t+"."+filtered_hashes[t][a]+"/statistics_filtered_"+t+"_"+a+"-"+str(i)+"-"+str(config["concurrency"])+".txt")
	return lis

def pull_algn_info(wildcards):
	lis = []
	for t in filtered_hashes.keys():
		for a in filtered_hashes[t].keys():
			for i in range(1, config["concurrency"] + 1):
				lis.append("results/statistics/filter-"+a+"-"+t+"."+filtered_hashes[t][a]+"/alignment_filter_information_"+t+"_"+a+".txt")
	return lis

def get_trimmer_params(wildcards):
	algn_trim = expand("{al}__{tr}", al=config["alignment"]["method"], tr=config["trimming"]["method"])
	out = []
	for comb in algn_trim:
		trimmer = comb.split("_")[-1]
		if trimmer == "untrimmed":
			config["trimming"]["options"][trimmer] = ""
			#out.append(comb + "__no§parameters") # this is a little weird hack because of the bash for loop below. but it works.
		if config["trimming"]["options"][trimmer] == None or config["trimming"]["options"][trimmer] == "":
			out.append(comb + "__no§parameters") # this is a little weird hack because of the bash for loop below. but it works.
		else:
			out.append(comb + "__" + config["trimming"]["options"][trimmer])
	return out

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
		minsp = config["trimming"]["minsp"],
		seq_type = config["filtering"]["seq_type"],
		min_parsimony_sites = config["trimming"]["min_parsimony_sites"],
		algn_trim = get_trimmer_params
	shell:
		"""
		echo "combo\t$(head -n 1 {input.filtered[0]} | cut -f 1,4-)" > {output.stats}
		for f in {input.filtered}; do combo=$(echo $f | cut -d "/" -f 3 | sed 's/filter-//'); tail -n +2 $f | sed "s/^/$combo\t/"; done | cut -f 1,2,4- >> {output.stats}

		# gather aligner and trimming combinations and corresponding settings for statistics:	
		combs="{params.algn_trim}"
		echo "aligner\ttrimmer\ttrimming_params\tdupseq\tcutoff\tminsp\tseq_type\tmin_parsimony_sites" > results/statistics/trim-filter-overview.txt
		for combination in ${{combs}}; do
			out=$combination"__{params.dupseq}__{params.cutoff}__{params.minsp}__{params.seq_type}__{params.min_parsimony_sites}"
			echo $out | sed 's/__/\t/g' | sed 's/§/ /g' >> results/statistics/trim-filter-overview.txt
		done  
		echo "$(date) - phylociraptor filter-align done." >> results/statistics/runlog.txt
		touch {output.checkpoint}
		"""	

for method in config["trimming"]["method"]:
	include: "trimmers/" + method + ".smk"
