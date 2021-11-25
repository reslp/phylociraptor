configfile:  "data/config.yaml"

include: "concatenate.smk"

rule partition_alignment:
	input:
		rules.concatenate.output.statistics
	output:
		partitions = "results/phylogeny-{bootstrap}/concatenate/{aligner}-{alitrim}/partitions.txt"
	params:
		wd = os.getcwd(),
		models = "results/modeltest/best_models_{aligner}_{alitrim}.txt",
		datatype = config["filtering"]["seq_type"]
	shell:
		"""
		if [[ -f {params.wd}/{params.models} && {params.wd}/checkpoints/modeltest.done ]]; then
			echo "$(date) - 'phylociraptor model' finished successfully before. Will run raxml with best models." >> {params.wd}/results/statistics/runlog.txt
			awk 'FNR==NR{{a[$1"_aligned_trimmed.fas"]=$2;next}}{{print $0"\\t"a[$1]}}' {params.models} results/phylogeny-{wildcards.bootstrap}/concatenate/{wildcards.aligner}-{wildcards.alitrim}/statistics.txt | awk -F"\\t" 'NR>1{{split($1,b,"_"); print $9", " b[1]"="$2"-"$3}}' > results/phylogeny-{wildcards.bootstrap}/concatenate/{wildcards.aligner}-{wildcards.alitrim}/partitions_unformated.txt
		else
			echo "$(date) - 'phylociraptor model' was NOT run before. Will run raxml with GTR or PROTGTR depending on input data type." >> {params.wd}/results/statistics/runlog.txt
			if [[ {params.datatype} == "aa" ]]; then
				awk '{{print $0"\\tPROTGTR"}}' {input} | awk -F"\\t" 'NR>1{{split($1,b,"_"); print $9", " b[1]"="$2"-"$3}}' > results/phylogeny-{wildcards.bootstrap}/concatenate/{wildcards.aligner}-{wildcards.alitrim}/partitions_unformated.txt
			else
				awk '{{print $0"\\tGTR"}}' {input} | awk -F"\\t" 'NR>1{{split($1,b,"_"); print $9", " b[1]"="$2"-"$3}}' > results/phylogeny-{wildcards.bootstrap}/concatenate/{wildcards.aligner}-{wildcards.alitrim}/partitions_unformated.txt
			fi
		fi
		# correct some model names to make them raxml compatible:
		# it is not quite clear which models are compatible. During more extensive tests additional problems should show up
		cat results/phylogeny-{wildcards.bootstrap}/concatenate/{wildcards.aligner}-{wildcards.alitrim}/partitions_unformated.txt | sed 's/JTTDCMut/JTT-DCMut/' > {output.partitions}
		touch {output.partitions}
		"""

rule raxmlng:
		input:
			alignment = rules.concatenate.output.alignment,
			partitions = rules.partition_alignment.output.partitions
		output:
			checkpoint = "results/checkpoints/raxml_{aligner}_{alitrim}_{bootstrap}.done",
			alignment = "results/phylogeny-{bootstrap}/raxml/{aligner}-{alitrim}/concat.fas",
			partitions = "results/phylogeny-{bootstrap}/raxml/{aligner}-{alitrim}/partitions.txt",
			statistics = "results/statistics/mltree/mltree_raxml_{aligner}_{alitrim}_{bootstrap}_statistics.txt"
		benchmark:
			"results/statistics/benchmarks/tree/raxmlng_{aligner}_{alitrim}_{bootstrap}.txt"
		singularity:
			"docker://reslp/raxml-ng:1.0.0"
		threads: config["raxmlng"]["threads"]
		params:
			bs = config["raxmlng"]["bootstrap"],
			wd = os.getcwd(),
			additional_params = config["raxmlng"]["additional_params"]
		shell:
			"""
			cp {input.alignment} {output.alignment}
			cp {input.partitions} {output.partitions}
			cd results/phylogeny-{wildcards.bootstrap}/raxml/{wildcards.aligner}-{wildcards.alitrim}
			raxml-ng --msa {params.wd}/{output.alignment} --prefix raxmlng -threads {threads} --bs-trees {params.bs} --model {params.wd}/{output.partitions} --all {params.additional_params}
			statistics_string="raxmlng\t{wildcards.aligner}\t{wildcards.alitrim}\t{params.bs}\t{wildcards.bootstrap}\t$(cat {params.wd}/{output.partitions} | wc -l)\t$(cat raxmlng.raxml.bestTree)"
			echo -e $statistics_string > {params.wd}/{output.statistics}
			touch {params.wd}/{output.checkpoint}
			"""
rule iqtree:
		input:
#			"results/statistics/filter-{aligner}-{alitrim}/alignment_filter_information_{alitrim}_{aligner}.txt" 
			"results/checkpoints/modes/modeltest.done"
		output:
			checkpoint = "results/checkpoints/iqtree_{aligner}_{alitrim}_{bootstrap}.done",
			statistics = "results/statistics/mltree/mltree_iqtree_{aligner}_{alitrim}_{bootstrap}_statistics.txt"
		benchmark:
			"results/statistics/benchmarks/tree/iqtree_{aligner}_{alitrim}_{bootstrap}.txt"
		singularity:
			"docker://reslp/iqtree:2.0.7"
		params:
			wd = os.getcwd(),
			models = "results/modeltest/best_models_{aligner}_{alitrim}.txt",
			nt = "AUTO",
			bb = config["iqtree"]["bootstrap"],
			m = config["iqtree"]["model"],
			maxmem = config["iqtree"]["maxmem"],
			additional_params = config["iqtree"]["additional_params"],
#			bootstrap_cutoff_file = "results/statistics/genetree_filter_{aligner}_{alitrim}.txt",
			genes = get_input_genes
		threads:
			config["iqtree"]["threads"]
		shell:
			"""
			rm -rf results/phylogeny-{wildcards.bootstrap}/iqtree/{wildcards.aligner}-{wildcards.alitrim}/algn
			if [[ ! -d results/phylogeny-{wildcards.bootstrap}/iqtree/{wildcards.aligner}-{wildcards.alitrim} ]]; then mkdir -p results/phylogeny-{wildcards.bootstrap}/iqtree/{wildcards.aligner}-{wildcards.alitrim}; fi
			cd results/phylogeny-{wildcards.bootstrap}/iqtree/{wildcards.aligner}-{wildcards.alitrim}/
			mkdir algn
			echo "$(date) - iqtree {wildcards.aligner}-{wildcards.alitrim}: Will use bootstrap cutoff ({wildcards.bootstrap}) before creating concatenated alignment" >> {params.wd}/results/statistics/runlog.txt
			for gene in $(echo "{params.genes}")
			do
				cp {params.wd}/results/alignments/filtered/{wildcards.aligner}-{wildcards.alitrim}/"$gene"_aligned_trimmed.fas algn
			done
			# here we decide how iqtree should be run. In case modeltesting was run before, this will not be repeated here.				
			echo "Will create NEXUS partition file with model information now."
			echo "#nexus" > concat.nex
			echo "begin sets;" >> concat.nex 
			i=1
			charpart=""
			for gene in $(echo "{params.genes}")
			do
				cat {params.wd}/{params.models} | grep $gene | awk -v i=$i '{{print "charset part"i" = algn/"$1"_aligned_trimmed.fas:*;"}}' >> concat.nex
				charpart=$charpart$(cat {params.wd}/{params.models} | grep $gene | awk -v i=$i '{{printf($2":part"i", ")}}' | sed 's/\\(.*\\), /\\1, /')
				i=$((++i))
			done
			echo "charpartition mine = "$charpart";" >> concat.nex
			echo "end;" >> concat.nex
			echo "$(date) - nexus file for iqtree written." >> {params.wd}/results/statistics/runlog.txt
			if [[ -z "{params.maxmem}" ]]; then
				iqtree -p concat.nex --prefix concat -bb {params.bb} -nt AUTO -ntmax {threads} -redo {params.additional_params}
			else
				iqtree -p concat.nex --prefix concat -bb {params.bb} -nt AUTO -ntmax {threads} -redo -mem {params.maxmem} {params.additional_params}
			fi
			#rm -r algn
			statistics_string="iqtree\t{wildcards.aligner}\t{wildcards.alitrim}\t{params.bb}\t{wildcards.bootstrap}\t$(ls algn | wc -l)\t$(cat concat.contree)"
			echo -e $statistics_string > {params.wd}/{output.statistics}	
			cd {params.wd}
			touch {output.checkpoint}
			"""

if "phylobayes" in config["tree"]["method"]:
	chains = [str(i) for i in range(1,int(config["phylobayes"]["nchains"])+1)]
	rule phylobayes:
		input:
			alignment = rules.concatenate.output.phylip_alignment
		output:
			checkpoint = "results/checkpoints/phylobayes_chain{chain}.done",
			trace = "results/phylogeny/phylobayes/phylobayes_chain{chain}.trace",
			treelist = "results/phylogeny/phylobayes/phylobayes_chain{chain}.treelist"
		singularity:
			"docker://reslp/phylobayes-mpi:git_dca7bdf"
		params:
			model = config["phylobayes"]["model"],
			threads = config["phylobayes"]["threads"],
			ngens = config["phylobayes"]["ngens"],
			additional_params = config["phylobayes"]["additional_params"],
			chain = "{chain}",
			wd = os.getcwd()
		shell:
			"""
			#mkdir -p results/phylogeny/phylobayes
			cp -n {input.alignment} results/phylogeny/phylobayes/concat.phy
			cd results/phylogeny/phylobayes
			unset PE_HOSTFILE
			mpirun -np {params.threads} pb_mpi -d concat.phy {params.model} -x 1 {params.ngens} {params.additional_params} phylobayes_chain{params.chain}
			cd {params.wd}
			touch {output.checkpoint}
			"""
	rule merge_phylobayes_chains:
		input:
			expand("results/checkpoints/phylobayes_chain{chain}.done", chain = chains)	
		output:
			checkpoint = "results/checkpoints/merge_phylobayes_chains.dome"
		shell:
			"""
			touch {output.checkpoint}
			"""
else:
	chains = [str(i) for i in range(1,int(config["phylobayes"]["nchains"])+1)]
	rule phylobayes:
		output:
			checkpoint = "results/checkpoints/phylobayes_chain{chain}.done"
		shell:
			"""
			touch {output.checkpoint}
			"""
	rule merge_phylobayes_chains:
		input:
			expand("results/checkpoints/phylobayes_chain{chain}.done", chain = chains)
		output:
			checkpoint = "results/checkpoints/merge_phylobayes_chains.dome"
		shell:
			"""
			touch {output.checkpoint}
			"""
rule all_trees:
	input:
		expand("results/checkpoints/{treeinfer}_{aligner}_{alitrim}_{bootstrap}.done", aligner=config["alignment"]["method"], alitrim=config["trimming"]["method"], treeinfer=config["tree"]["method"], bootstrap=config["filtering"]["bootstrap_cutoff"])
	output:
		"results/checkpoints/modes/trees.done"
	shell:
		"""
		touch {output}
		echo "$(date) - Pipeline part 3 (tree) done." >> results/statistics/runlog.txt
		"""
