include: "functions.smk"
import os.path
import glob
import yaml

# get list of containers to use:
with open("data/containers.yaml", "r") as yaml_stream:
    containers = yaml.safe_load(yaml_stream)

BUSCOS, = glob_wildcards("results/orthology/busco/busco_sequences_deduplicated/{busco}_all.fas")

aligners = get_aligners()		
trimmers = get_trimmers()		
tree_methods = get_treemethods()
bscuts = get_bootstrap_cutoffs()

#BUSCOS, = glob_wildcards("results/filtered_alignments/"+config["alignment"]["method"][0]+"/"+config["trimming"]["method"][0]+"/{busco}_aligned_trimmed.fas")

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
		if os.path.isfile("results/alignments/filtered/"+wildcards.aligner+"-"+wildcards.alitrim+"/"+str(busco)+"_aligned_trimmed.fas"):
#			lis.append("results/modeltest/"+wildcards.aligner+"-"+wildcards.alitrim+"/"+busco+"/"+busco+".log")
			lis.append("results/checkpoints/modeltest/"+wildcards.aligner+"-"+wildcards.alitrim+"/"+str(busco)+"_modeltest_"+wildcards.aligner+"_"+wildcards.alitrim+".done")
	return lis

def return_genes(wildcards):
	lis = []
        for busco in BUSCOS:
		if os.path.isfile("results/alignments/filtered/"+wildcards.aligner+"-"+wildcards.alitrim+"/"+str(busco)+"_aligned_trimmed.fas"):
#			lis.append("results/modeltest/"+wildcards.aligner+"-"+wildcards.alitrim+"/"+busco+"/"+busco+".log")
			lis.append(busco)
	return lis

rule iqtree_mt:
	input:
#		alignment = get_alignment
		alignment = "results/alignments/filtered/{aligner}-{alitrim}/{busco}_aligned_trimmed.fas"
	output:
#		logfile = "results/modeltest/{aligner}-{alitrim}/{busco}/{busco}.log",
		checkpoint = "results/checkpoints/modeltest/{aligner}-{alitrim}/{busco}_modeltest_{aligner}_{alitrim}.done"
	benchmark:
		"results/statistics/benchmarks/model/modeltest_{busco}_{aligner}_{alitrim}.txt"
	params:
		wd = os.getcwd(),
		busco = "{busco}",
		bb = config["genetree"]["bootstrap"],
		additional_params = config["iqtree"]["additional_params"],
		seed = config["seed"]
	singularity:
		containers["iqtree"]	
	threads: config["modeltest"]["threads"]
	shell:
		"""
		if [[ ! -d "results/modeltest/{wildcards.aligner}-{wildcards.alitrim}/{wildcards.busco}" ]]; then mkdir -p results/modeltest/{wildcards.aligner}-{wildcards.alitrim}/{wildcards.busco}; fi
		iqtree -s {input.alignment} -msub nuclear --prefix {params.wd}/results/modeltest/{wildcards.aligner}-{wildcards.alitrim}/{params.busco}/{params.busco} -nt {threads} -m MFP $(if [[ "{params.bb}" != "None" ]]; then echo "-bb {params.bb}"; fi) $(if [[ "{params.seed}" != "None" ]]; then echo "-seed {params.seed}"; fi) {params.additional_params}
		touch {output.checkpoint}
		"""

rule aggregate_best_models:
	input:
		checkfiles = return_target_modeltest_check
	output:
		best_models = "results/modeltest/best_models_{aligner}_{alitrim}.txt",
		checkpoint = "results/checkpoints/modeltest/aggregate_best_models_{aligner}_{alitrim}.done",
		genetree_filter_stats = "results/modeltest/genetree_filter_{aligner}_{alitrim}.txt"
	benchmark:
		"results/statistics/benchmarks/model/aggregate_best_models_{aligner}_{alitrim}.txt"
	params:
		wd = os.getcwd(),
		genes = return_genes
	shell:
		"""
		# echo "name\tmodel" > {output.best_models}
		for gene in {params.genes}
		do
			model=$(cat results/modeltest/{wildcards.aligner}-{wildcards.alitrim}/$gene/$gene.log | grep "Best-fit model:" | tail -n 1 | awk -F ":" '{{print $2}}' | awk -F " " '{{print $1}}')
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
			bootstrapvalues=$(grep -E '\)[0-9]+:' -o results/modeltest/{wildcards.aligner}-{wildcards.alitrim}/$gene/$gene.treefile | sed 's/)//' | sed 's/://' | tr '\n' '+' | sed 's/+$//')
#			bootstrapsum=$(echo "$bootstrapvalues" | bc)
			bootstrapsum=$(echo "$bootstrapvalues" | perl -ne 'chomp; @bs=split(/\+/); $sum=0; for (@bs){{$sum=$sum+$_}}; print "$sum"')
		
			totalbootstraps=$(echo "$bootstrapvalues" | sed 's/[0-9]//g' | wc -c)
#			meanbootstrap=$(echo "($bootstrapsum)/$totalbootstraps" | bc)
			meanbootstrap=$(( bootstrapsum / totalbootstraps ))
			echo -e "$gene\t{wildcards.aligner}\t{wildcards.alitrim}\t$bootstrapsum\t$totalbootstraps\t$meanbootstrap" | tee -a {output.genetree_filter_stats}
		done
		touch {output.checkpoint}

		"""

rule modeltest:
	input:
		expand("results/checkpoints/modeltest/aggregate_best_models_{aligner}_{alitrim}.done", aligner=aligners, alitrim=trimmers),
		expand("results/modeltest/genetree_filter_{aligner}_{alitrim}.txt", aligner=aligners, alitrim=trimmers)
	output:
		"results/checkpoints/modes/modeltest.done"
	shell:
		"""
		touch {output}
		echo "$(date) - Pipeline part modeltest (model) done." >> results/statistics/runlog.txt
		"""
