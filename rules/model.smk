include: "functions.smk"
import os.path
import glob
import yaml

ruleorder: read_params_global > read_params_per

# get list of containers to use:
with open("data/containers.yaml", "r") as yaml_stream:
    containers = yaml.safe_load(yaml_stream)

#create new hashes for current stage 
hashes = collect_hashes("modeltest")

filter_orthology_hash = hashes['filter-orthology']
aligner_hashes = hashes['align']["per"]
trimmer_hashes = hashes['filter-align']["per"]
previous_hash = hashes['filter-align']["global"]

modeltest_hashes = hashes['modeltest']["per"]
current_hash = hashes['modeltest']["global"]

BUSCOS, = glob_wildcards("results/orthology/busco/busco_sequences_deduplicated."+filter_orthology_hash+"/{busco}_all.fas")

def previous_params_global(wildcards):
	return "results/alignments/trimmed/parameters.filter-align."+previous_hash+".yaml"
def previous_params_per(wildcards):
	return "results/alignments/trimmed/"+wildcards.aligner+"-"+wildcards.alitrim+"."+trimmer_hashes[wildcards.alitrim][wildcards.aligner]+"/parameters.filter-align."+wildcards.aligner+"-"+wildcards.alitrim+"."+trimmer_hashes[wildcards.alitrim][wildcards.aligner]+".yaml"

aligners = get_aligners()		
trimmers = get_trimmers()		
tree_methods = get_treemethods()
bscuts = get_bootstrap_cutoffs()

# keep old modeltest rule for later reference. The Alignment info part may be necessary for raxml later...
#rule modeltest:
#	input:
#		rules.part2.output
#	output:
#		models = "results/modeltest/best_models.txt",
#		checkpoint = "results/checkpoints/modeltest.done"
#	params:
#		wd = os.getcwd()
#	singularity:
#		containers["iqtree"]	
#	threads: 16
#	shell:
#		"""
#		mkdir -p results/modeltest
#		echo "name\tmodel\tnseq\tnsites\tconstsites\tinvsites\tparssites\tdistpatt" > {params.wd}/results/modeltest/best_models.txt
#		for file in $(ls results/filtered_alignments/*.fas);
#		do
#			outname=$(basename $file)
#			iqtree -m TESTONLY -s $file -msub nuclear -redo --prefix {params.wd}/results/modeltest/$outname -nt AUTO -ntmax {threads}
#			printf "$outname\t" >> {params.wd}/results/modeltest/best_models.txt
#			cat {params.wd}/results/modeltest/$outname.log | grep "Best-fit model:" | awk -F ":" '{{print $2}}' | awk -F " " '{{printf ("%s\t", $1)}}' >> {params.wd}/results/modeltest/best_models.txt
#			cat {params.wd}/results/modeltest/$outname.iqtree | grep "Input data" | awk -F ":" '{{print $2}}' | awk -F " " '{{printf ("%s\t", $1)}}' >> {params.wd}/results/modeltest/best_models.txt
#			cat {params.wd}/results/modeltest/$outname.iqtree | grep "Input data" | awk -F ":" '{{print $2}}' | awk -F " " '{{printf ("%s\t", $4)}}' >> {params.wd}/results/modeltest/best_models.txt
#			cat {params.wd}/results/modeltest/$outname.iqtree | grep "Number of constant sites" | awk -F ":" '{{print $2}}' | awk -F " " '{{printf ("%s\t", $1)}}' >> {params.wd}/results/modeltest/best_models.txt
#			cat {params.wd}/results/modeltest/$outname.iqtree | grep "Number of invariant" | awk -F ":" '{{print $2}}' | awk -F " " '{{printf ("%s\t", $1)}}' >> {params.wd}/results/modeltest/best_models.txt
#			cat {params.wd}/results/modeltest/$outname.iqtree | grep "Number of parsimony" | awk -F ":" '{{print $2}}' | awk -F " " '{{printf ("%s\t", $1)}}' >> {params.wd}/results/modeltest/best_models.txt
#			cat {params.wd}/results/modeltest/$outname.iqtree | grep "Number of distinct" | awk -F ":" '{{print $2}}' | awk -F " " '{{print $1}}' >> {params.wd}/results/modeltest/best_models.txt
#			# currently iqtree output will be removed and only the best model is saved.
#			rm {params.wd}/results/modeltest/$outname.*
#		done
#		touch {output.checkpoint}
#		"""

def return_target_modeltest_check(wildcards):
	lis = []
        for busco in BUSCOS:
		if os.path.isfile("results/alignments/filtered/"+wildcards.aligner+"-"+wildcards.alitrim+"."+trimmer_hashes[wildcards.alitrim][wildcards.aligner]+"/"+str(busco)+"_aligned_trimmed.fas"):
#			lis.append("results/modeltest/"+wildcards.aligner+"-"+wildcards.alitrim+"."+modeltest_hashes["iqtree"][wildcards.alitrim][wildcards.aligner]+"/"+busco+"/"+busco+".log")
			lis.append("results/checkpoints/modeltest/"+wildcards.aligner+"-"+wildcards.alitrim+"."+modeltest_hashes["iqtree"][wildcards.alitrim][wildcards.aligner]+"/"+str(busco)+"_modeltest_"+wildcards.aligner+"_"+wildcards.alitrim+".done")
	return lis

def return_genes(wildcards):
	lis = []
        for busco in BUSCOS:
		if os.path.isfile("results/alignments/filtered/"+wildcards.aligner+"-"+wildcards.alitrim+"."+trimmer_hashes[wildcards.alitrim][wildcards.aligner]+"/"+str(busco)+"_aligned_trimmed.fas"):
#			lis.append("results/modeltest/"+wildcards.aligner+"-"+wildcards.alitrim+"."+modeltest_hashes["iqtree"][wildcards.alitrim][wildcards.aligner]+"/"+busco+"/"+busco+".log")
			lis.append(busco)
	return lis

def compare_modeltest(wildcards):
	return [trigger("results/modeltest/parameters.modeltest.{aligner}-{alitrim}.{hash}.yaml".format(aligner=wildcards.aligner, alitrim=wildcards.alitrim, hash=wildcards.hash), configfi)]	

rule read_params_per:
	input:
		trigger = compare_modeltest,
		previous = previous_params_per
	output:
		"results/modeltest/parameters.modeltest.{aligner}-{alitrim}.{hash}.yaml"
	shell:
		"""
		bin/read_write_yaml.py {input.trigger} {output} seed modeltest,options,iqtree modeltest,bootstrap
		cat {input.previous} >> {output}
		"""
rule read_params_global:
	input:
		trigger = compare("results/modeltest/parameters.modeltest."+current_hash+".yaml", configfi),
		previous = previous_params_global
	output:
		"results/modeltest/parameters.modeltest."+current_hash+".yaml"
	shell:
		"""
		bin/read_write_yaml.py {input.trigger} {output} seed modeltest,method modeltest,options modeltest,bootstrap
		cat {input.previous} >> {output}
		"""

def alignments_in(wildcards):
	return "results/alignments/filtered/"+wildcards.aligner+"-"+wildcards.alitrim+"."+trimmer_hashes[wildcards.alitrim][wildcards.aligner]+"/"+wildcards.busco+"_aligned_trimmed.fas"

def return_trimmer_checkpoint(wildcards):
	return "results/statistics/filter-"+wildcards.aligner+"-"+wildcards.alitrim+"."+trimmer_hashes[wildcards.alitrim][wildcards.aligner]+"/alignment_filter_information_"+wildcards.alitrim+"_"+wildcards.aligner+".txt"

rule iqtree_mt:
	input:
		"results/modeltest/parameters.modeltest.{aligner}-{alitrim}.{hash}.yaml",
		alignment = alignments_in,
		checkpoint = return_trimmer_checkpoint,
	output:
		checkpoint = "results/checkpoints/modeltest/{aligner}-{alitrim}.{hash}/{busco}_modeltest_{aligner}_{alitrim}.done"
	benchmark:
		"results/statistics/benchmarks/model/modeltest_{busco}_{aligner}_{alitrim}.{hash}.txt"
	params:
		wd = os.getcwd(),
		busco = "{busco}",
		bb = config["modeltest"]["bootstrap"],
		additional_params = config["modeltest"]["options"]["iqtree"],
		target_dir = "results/modeltest/{aligner}-{alitrim}.{hash}/{busco}",
		seed = config["seed"]
	singularity:
		containers["iqtree"]	
	log:
		"log/modeltest/iqtree/iqtree_mt-{aligner}-{alitrim}-{busco}.{hash}.txt"
	threads: config["modeltest"]["threads"]
	shell:
		"""
		if [[ ! -d {params.target_dir} ]]; then mkdir -p {params.target_dir}; fi
		iqtree -s {input.alignment} -msub nuclear --prefix {params.wd}/{params.target_dir}/{params.busco} -nt {threads} -m MFP $(if [[ "{params.bb}" != "None" ]]; then echo "-bb {params.bb}"; fi) $(if [[ "{params.seed}" != "None" ]]; then echo "-seed {params.seed}"; fi) {params.additional_params} 2>&1 | tee {log}
		touch {output.checkpoint}
		"""

rule aggregate_best_models:
	input:
		checkfiles = return_target_modeltest_check
	output:
		best_models = "results/modeltest/best_models_{aligner}_{alitrim}.{hash}.txt",
		checkpoint = "results/checkpoints/modeltest/aggregate_best_models_{aligner}_{alitrim}.{hash}.done",
		genetree_filter_stats = "results/modeltest/genetree_filter_{aligner}_{alitrim}.{hash}.txt"
	benchmark:
		"results/statistics/benchmarks/model/aggregate_best_models_{aligner}_{alitrim}.{hash}.txt"
	params:
		wd = os.getcwd(),
		genes = return_genes
	log:
		"log/modeltest/aggregate_best_models-{aligner}-{alitrim}.{hash}.txt"
	shell:
		"""
		# echo "name\tmodel" > {output.best_models}
		for gene in {params.genes}
		do
			model=$(cat results/modeltest/{wildcards.aligner}-{wildcards.alitrim}.{wildcards.hash}/$gene/$gene.log | grep "Best-fit model:" | tail -n 1 | awk -F ":" '{{print $2}}' | awk -F " " '{{print $1}}')
			# This should not happen
			if [[ $(echo $model | wc -l) -eq 0 ]]
			then
				#this should not happen - just a failsafe
				echo "ERROR: No model detected for gene: $gene" >&2
			else
				printf "$gene\t" >> {output.best_models}
				echo $model >> {output.best_models}
			fi

			#now calculate mean bootsrap for each tree
			bootstrapvalues=$(grep -E '\)[0-9]+:' -o results/modeltest/{wildcards.aligner}-{wildcards.alitrim}.{wildcards.hash}/$gene/$gene.treefile | sed 's/)//' | sed 's/://' | tr '\n' '+' | sed 's/+$//')
#			bootstrapsum=$(echo "$bootstrapvalues" | bc)
			bootstrapsum=$(echo "$bootstrapvalues" | perl -ne 'chomp; @bs=split(/\+/); $sum=0; for (@bs){{$sum=$sum+$_}}; print "$sum"')
		
			totalbootstraps=$(echo "$bootstrapvalues" | sed 's/[0-9]//g' | wc -c)
#			meanbootstrap=$(echo "($bootstrapsum)/$totalbootstraps" | bc)
			meanbootstrap=$(( bootstrapsum / totalbootstraps ))
			echo -e "$gene\t{wildcards.aligner}\t{wildcards.alitrim}\t$bootstrapsum\t$totalbootstraps\t$meanbootstrap" | tee -a {output.genetree_filter_stats}
		done 2>&1 | tee {log}
		touch {output.checkpoint}

		"""

def pull(wildcards):
	lis = []
	for m in config["modeltest"]["method"]:
		for a in config["alignment"]["method"]:
			for t in config["trimming"]["method"]:
				lis.append("results/checkpoints/modeltest/aggregate_best_models_"+a+"_"+t+"."+modeltest_hashes[m][t][a]+".done")
				lis.append("results/modeltest/genetree_filter_"+a+"_"+t+"."+modeltest_hashes[m][t][a]+".txt")
				lis.append("results/modeltest/parameters.modeltest."+a+"-"+t+"."+modeltest_hashes[m][t][a]+".yaml")
	return lis

rule modeltest:
	input:
		pull,
		params = rules.read_params_global.output,
	output:
		"results/checkpoints/modes/modeltest."+current_hash+".done"
	shell:
		"""
		touch {output}
		echo "$(date) - phylociraptor modeltest done." >> results/statistics/runlog.txt
		"""
