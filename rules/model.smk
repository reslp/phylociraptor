BUSCOS, = glob_wildcards("results/filtered_alignments/{busco}_aligned_trimmed.fas")

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
#		"docker://reslp/iqtree:2.0.7"
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


rule modeltest:
	input:
		rules.part2.output,
		alignment = "results/filtered_alignments/{busco}_aligned_trimmed.fas"
	output:
		logfile = "results/modeltest/{busco}/{busco}.log",
		checkpoint = "results/checkpoints/modeltest/{busco}_modeltest.done"
	benchmark:
		"results/statistics/benchmarks/model/modeltest_{busco}.txt"
	params:
		wd = os.getcwd(),
		busco = "{busco}"
	singularity:
		"docker://reslp/iqtree:2.0.7"
	threads: config["iqtree"]["threads"]
	shell:
		"""
		iqtree -m TESTONLY -s {input.alignment} -msub nuclear -redo --prefix {params.wd}/results/modeltest/{params.busco}/{params.busco} -nt AUTO -ntmax {threads}
		touch {output.checkpoint}
		"""

rule aggregate_best_models:
	input:
		checkpoints = expand("results/checkpoints/modeltest/{b}_modeltest.done", b = BUSCOS),
		logfiles = expand("results/modeltest/{b}/{b}.log", b = BUSCOS)
	output:
		best_models = "results/modeltest/best_models.txt",
		checkpoint = "results/checkpoints/modeltest/aggregate_best_models.done"
	benchmark:
		"results/statistics/benchmarks/model/aggregate_best_models.txt"
	params:
		wd = os.getcwd()
	shell:
		"""
		#echo "name\tmodel" > {output.best_models}
		for file in {input.logfiles};
		do
		outname=$(basename $file | awk -F "." '{{print $1}}')
		printf "$outname\t" >> {output.best_models}
		cat $file | grep "Best-fit model:" | awk -F ":" '{{print $2}}' | awk -F " " '{{print $1}}' >> {output.best_models}
		done
		touch {output.checkpoint}
		"""

