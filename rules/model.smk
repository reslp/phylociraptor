import os.path
import glob
import yaml

configfile: "data/config.yaml"

# get list of containers to use:
with open("data/containers.yaml", "r") as yaml_stream:
    containers = yaml.safe_load(yaml_stream)

BUSCOS, = glob_wildcards("results/orthology/busco/busco_sequences_deduplicated/{busco}_all.fas")

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

def return_target_modeltest_log(wildcards):
	lis = []
        for busco in BUSCOS:
		if os.path.isfile("results/alignments/filtered/"+wildcards.aligner+"-"+wildcards.alitrim+"/"+str(busco)+"_aligned_trimmed.fas"):
			lis.append("results/modeltest/"+wildcards.aligner+"-"+wildcards.alitrim+"/"+busco+"/"+busco+".log")
#	input_files = glob.glob("results/alignments/filtered/"+wildcards.aligner+"-"+wildcards.alitrim+"/*_aligned_trimmed.fas")
#	for f in input_files:
#		busco = f.split("/")[-1].replace("_aligned_trimmed.fas","")
#		lis.append("results/modeltest/"+busco+"/"+wildcards.aligner+"-"+wildcards.alitrim+"/"+busco+".log")
	return lis

rule modeltest:
	input:
#		alignment = get_alignment
		alignment = "results/alignments/filtered/{aligner}-{alitrim}/{busco}_aligned_trimmed.fas"
	output:
		logfile = "results/modeltest/{aligner}-{alitrim}/{busco}/{busco}.log",
		checkpoint = "results/checkpoints/modeltest/{aligner}-{alitrim}/{busco}_modeltest_{aligner}_{alitrim}.done"
	benchmark:
		"results/statistics/benchmarks/model/modeltest_{busco}_{aligner}_{alitrim}.txt"
	params:
		wd = os.getcwd(),
		busco = "{busco}",
		maxmem = config["iqtree"]["maxmem"],
		bb = config["genetree"]["boostrap"],
		additional_params = config["iqtree"]["additional_params"]
	singularity:
		containers["iqtree"]	
	threads: config["modeltest"]["threads"]
	shell:
		"""
		iqtree -s {input.alignment} -msub nuclear -redo --prefix {params.wd}/results/modeltest/{wildcards.aligner}-{wildcards.alitrim}/{params.busco}/{params.busco} -nt {threads} -m MFP $(if [[ "{params.bb}" != "None" ]]; then echo "-bb {params.bb}"; fi) {params.additional_params}
		touch {output.checkpoint}
		"""

rule aggregate_best_models:
	input:
		logfiles = return_target_modeltest_log
	output:
		best_models = "results/modeltest/best_models_{aligner}_{alitrim}.txt",
		checkpoint = "results/checkpoints/modeltest/aggregate_best_models_{aligner}_{alitrim}.done",
		genetree_filter_stats = "results/modeltest/genetree_filter_{aligner}_{alitrim}.txt"
	benchmark:
		"results/statistics/benchmarks/model/aggregate_best_models_{aligner}_{alitrim}.txt"
	params:
		wd = os.getcwd(),
	shell:
		"""
		# echo "name\tmodel" > {output.best_models}
		for file in $(ls -1 results/modeltest/{wildcards.aligner}-{wildcards.alitrim}/*/*.log)
		do
			outname=$(basename $file | awk -F "." '{{print $1}}')
			printf "$outname\t" >> {output.best_models}
			cat $file | grep "Best-fit model:" | awk -F ":" '{{print $2}}' | awk -F " " '{{print $1}}' >> {output.best_models}
		done
		# now calculate mean bootstrap for each tree
		for gene in $(ls -d results/modeltest/{wildcards.aligner}-{wildcards.alitrim}/*)
				do
					bootstrapvalues=$(grep -E '\)[0-9]+:' -o $gene/*.treefile | sed 's/)//' | sed 's/://' | tr '\n' '+' | sed 's/+$//')
					bootstrapsum=$(echo "$bootstrapvalues" | bc)
					totalbootstraps=$(grep -E '\)[0-9]+:' -o $gene/*.treefile | wc -l )
					meanbootstrap=$(echo "($bootstrapsum)/$totalbootstraps" | bc)
					genename=$(basename $gene)
					echo -e "$genename\t{wildcards.aligner}\t{wildcards.alitrim}\t$bootstrapsum\t$totalbootstraps\t$meanbootstrap" | tee -a {output.genetree_filter_stats}
				done
		touch {output.checkpoint}
		"""

rule all_modeltest:
	input:
		expand("results/checkpoints/modeltest/aggregate_best_models_{aligner}_{alitrim}.done", aligner=config["alignment"]["method"], alitrim=config["trimming"]["method"]),
		expand("results/modeltest/genetree_filter_{aligner}_{alitrim}.txt", aligner=config["alignment"]["method"], alitrim=config["trimming"]["method"])
	output:
		"results/checkpoints/modes/modeltest.done"
	shell:
		"""
		touch {output}
		echo "$(date) - Pipeline part modeltest (model) done." >> results/statistics/runlog.txt
		"""
