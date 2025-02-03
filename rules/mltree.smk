include: "functions.smk"
include: "concatenate.smk"

os.environ["MODE"] = "mltree" # for communication with concatenate.smk

ruleorder: read_params_global > read_params_per

aligners = get_aligners()		
trimmers = get_trimmers()		
tree_methods = get_treemethods()
allbscuts = get_bootstrap_cutoffs()

#create new hashes for current stage 
hashes = collect_hashes("mltree", config, configfi, wd=os.getcwd())

filter_orthology_hash = hashes['filter-orthology']["global"]
aligner_hashes = hashes['align']["per"]
trimmer_hashes = hashes['filter-align']["per"]
modeltest_hashes = hashes['modeltest']["per"]
tinference_hashes = hashes['mltree']["per"]
current_hash = hashes['mltree']["global"]
previous_hash = hashes['modeltest']["global"]
previous_alitrim_hash = hashes['filter-align']["global"]

print(hashes["modeltest"])


############ functions specifically for this step
def compare_mltree(wildcards):
	return [trigger("results/phylogeny/{inference}/bootstrap-cutoff-{bootstrap}/parameters.mltree.{inference}-{aligner}-{alitrim}.{hash}.yaml".format(bootstrap=wildcards.bootstrap, inference=wildcards.inference, aligner=wildcards.aligner, alitrim=wildcards.alitrim, hash=wildcards.hash), configfi)]	

def get_iqtree_params(wildcards):
	return "results/phylogeny/iqtree/bootstrap-cutoff-" + wildcards.bootstrap + "/parameters.mltree.iqtree-"+wildcards.aligner+"-"+wildcards.alitrim+"."+wildcards.hash+".yaml"
############

rule read_params_per:
	input:
		trigger = compare_mltree,
		previous = previous_params_per
	output:
		"results/phylogeny/{inference}/bootstrap-cutoff-{bootstrap}/parameters.mltree.{inference}-{aligner}-{alitrim}.{hash}.yaml"
	shell:
		"""
		bin/read_write_yaml.py {input.trigger} {output} seed genetree_filtering,bootstrap_cutoff,{wildcards.bootstrap} mltree,options,{wildcards.inference} mltree,bootstrap,{wildcards.inference}
		cat {input.previous} >> {output}
		"""
rule read_params_global:
	input:
		trigger = compare("results/phylogeny/parameters.mltree."+current_hash+".yaml", configfi),
		previous = previous_params_global
	output:
		"results/phylogeny/parameters.mltree."+current_hash+".yaml"
	shell:
		"""
		bin/read_write_yaml.py {input.trigger} {output} seed genetree_filtering,bootstrap_cutoff mltree,method mltree,options mltree,bootstrap
		cat {input.previous} >> {output}
		"""


def get_best_models_filename(wildcards):
	hash = hashes["modeltest"]["per"]["iqtree"][wildcards.alitrim][wildcards.aligner]
	return "results/modeltest/best_models_" + wildcards.aligner + "_" + wildcards.alitrim + "." + hash + ".txt"

rule partition_alignment:
	input:
		rules.concatenate.output.statistics
	output:
		partitions = "results/phylogeny/concatenate/bootstrap-cutoff-{bootstrap}/{aligner}-{alitrim}.{hash}/partitions.txt"
	params:
		wd = os.getcwd(),
		models = get_best_models_filename,
		datatype = config["filtering"]["seq_type"]
	shell:
		"""
		if [[ -f {params.wd}/{params.models} ]]; then
			echo "$(date) - 'phylociraptor modeltest' finished successfully before. Will run raxml with best models." >> {params.wd}/results/statistics/runlog.txt
			awk 'FNR==NR{{a[$1"_aligned_trimmed.fas"]=$2;next}}{{print $0"\\t"a[$1]}}' {params.models} results/phylogeny/concatenate/bootstrap-cutoff-{wildcards.bootstrap}/{wildcards.aligner}-{wildcards.alitrim}.{wildcards.hash}/statistics.txt | awk -F"\\t" 'NR>1{{split($1,b,"_"); print $10", " b[1]"="$2"-"$3}}' > results/phylogeny/concatenate/bootstrap-cutoff-{wildcards.bootstrap}/{wildcards.aligner}-{wildcards.alitrim}.{wildcards.hash}/partitions_unformated.txt
		else
			echo "$(date) - 'phylociraptor modeltest' was NOT run before. Will run raxml with GTR or PROTGTR depending on input data type." >> {params.wd}/results/statistics/runlog.txt
			if [[ {params.datatype} == "aa" ]]; then
				awk '{{print $0"\\tPROTGTR"}}' {input} | awk -F"\\t" 'NR>1{{split($1,b,"_"); print $10", " b[1]"="$2"-"$3}}' > results/phylogeny/concatenate/bootstrap-cutoff-{wildcards.bootstrap}/{wildcards.aligner}-{wildcards.alitrim}.{wildcards.hash}/partitions_unformated.txt
			else
				awk '{{print $0"\\tGTR"}}' {input} | awk -F"\\t" 'NR>1{{split($1,b,"_"); print $10", " b[1]"="$2"-"$3}}' > results/phylogeny/concatenate/bootstrap-cutoff-{wildcards.bootstrap}/{wildcards.aligner}-{wildcards.alitrim}.{wildcards.hash}/partitions_unformated.txt
			fi
		fi
		# correct some model names to make them raxml compatible:
		# it is not quite clear which models are compatible. During more extensive tests additional problems should show up
		cat results/phylogeny/concatenate/bootstrap-cutoff-{wildcards.bootstrap}/{wildcards.aligner}-{wildcards.alitrim}.{wildcards.hash}/partitions_unformated.txt | sed 's/JTTDCMut/JTT-DCMut/' > {output.partitions}
		touch {output.partitions}
		"""

rule raxmlng:
		input:
			alignment = rules.concatenate.output.alignment,
			partitions = rules.partition_alignment.output.partitions
		output:
			checkpoint = "results/checkpoints/raxml_{aligner}_{alitrim}_{bootstrap}.{hash}.done",
			alignment = "results/phylogeny/raxml/bootstrap-cutoff-{bootstrap}/{aligner}-{alitrim}.{hash}/concat.fas",
			partitions = "results/phylogeny/raxml/bootstrap-cutoff-{bootstrap}/{aligner}-{alitrim}.{hash}/partitions.txt",
			statistics = "results/statistics/mltree/mltree_raxml_{aligner}_{alitrim}_{bootstrap}_statistics.{hash}.txt"
		benchmark:
			"results/statistics/benchmarks/tree/raxmlng_{aligner}_{alitrim}_{bootstrap}.{hash}.txt"
		singularity:
			containers["raxmlng"]	
		threads: config["mltree"]["threads"]["raxml"]
		params:
			bs = config["mltree"]["bootstrap"]["raxml"],
			wd = os.getcwd(),
			additional_params = config["mltree"]["options"]["raxml"],
			seed = config["seed"]
		log:
			"log/mltree/raxml/raxml-{aligner}-{alitrim}-{bootstrap}.{hash}.txt"
		shell:
			"""
			cp {input.alignment} {output.alignment}
			cp {input.partitions} {output.partitions}
			cd results/phylogeny/raxml/bootstrap-cutoff-{wildcards.bootstrap}/{wildcards.aligner}-{wildcards.alitrim}.{wildcards.hash}
			raxml-ng --msa {params.wd}/{output.alignment} $(if [[ "{params.seed}" != "None" ]]; then echo "--seed {params.seed}"; fi) --prefix raxmlng --threads auto{{{threads}}} --bs-trees {params.bs} --model {params.wd}/{output.partitions} --all {params.additional_params} 2>&1 | tee {params.wd}/{log}
			statistics_string="raxmlng\t{wildcards.aligner}\t{wildcards.alitrim}\t{params.bs}\t{wildcards.bootstrap}\t$(cat {params.wd}/{output.partitions} | wc -l)\t$(cat raxmlng.raxml.support)"
			echo -e $statistics_string > {params.wd}/{output.statistics}
			touch {params.wd}/{output.checkpoint}
			"""

rule prepare_iqtree:
		input:
			get_modeltest_checkpoint,
			get_iqtree_params
		output:
			algn = directory("results/phylogeny/iqtree/bootstrap-cutoff-{bootstrap}/{aligner}-{alitrim}.{hash}/algn"),
			nexus = "results/phylogeny/iqtree/bootstrap-cutoff-{bootstrap}/{aligner}-{alitrim}.{hash}/concat.nex"
		singularity:
			containers["iqtree"]
		params:
			wd = os.getcwd(),
			models = get_best_models,
			modelsfile = get_best_models_filename,
			alidir = get_alignment_dir, 
			genes = get_input_genes,
			iqtree_params = config["mltree"]["options"]["iqtree"]
		shell:
			"""
			mkdir {output.algn}
			echo "$(date) - prepare_iqtree {wildcards.aligner}-{wildcards.alitrim}: Will use bootstrap cutoff ({wildcards.bootstrap}) before creating concatenated alignment" >> {params.wd}/results/statistics/runlog.txt
			for gene in $(echo "{params.genes}")
			do
				cp {params.wd}/{params.alidir}/"$gene"_aligned_trimmed.fas {output.algn}/
			done
			
			echo "$(date) - Will create NEXUS partition file with model information now." >> {params.wd}/results/statistics/runlog.txt
			# check if -m is in iqtree_params

			
			if [[ -f {params.wd}/{params.models} ]]; then # modeltest was run
				if [[ $"{params.iqtree_params}" == *"-m "* ]]; then # model was specified (this overwrites best models from modeltest)
					model=$(echo "{params.iqtree_params}" | awk '{{ for(i=1;i<=NF;i++) {{ if($i=="-m") {{ if(i < NF) {{ print $(i+1); }}break; }}}}}}')
					if [[ "$model" == *"MFP"* || "$model" == *"TEST"* || "$model" == *"MERGE"* ]]; then
						model="None"
					fi
					echo "$(date) - 'phylociraptor modeltest' finished successfully before for {wildcards.aligner} and {wildcards.alitrim}, however a model was specified in the config file: $model. Will apply this model to all partitions." >> {params.wd}/results/statistics/runlog.txt
				else
					echo "$(date) - 'phylociraptor modeltest' finished successfully before for {wildcards.aligner} and {wildcards.alitrim}. Will run iqtree with best models for each gene/partition." >> {params.wd}/results/statistics/runlog.txt
					model="None"
				fi	 

				# write iqtree nexus partition file
				echo "#nexus" > {output.nexus}
				echo "begin sets;" >> {output.nexus}
				i=1
				charpart=""
					for gene in $(echo "{params.genes}")
					do
						cat {params.wd}/{params.models} | grep $gene | awk -v i=$i '{{print "charset part"i" = algn/"$1"_aligned_trimmed.fas:*;"}}' >> {output.nexus}
						if [[ "$model" == "None" ]]; then 
							charpart=$charpart$(cat {params.wd}/{params.models} | grep $gene | awk -v i=$i '{{printf($2":part"i", ")}}' | sed 's/\\(.*\\), /\\1, /')
						else
							charpart=$charpart$(cat {params.wd}/{params.models} | grep $gene | awk -v i=$i -v model=$model '{{printf(model":part"i", ")}}' | sed 's/\\(.*\\), /\\1, /')
						fi
						i=$((++i))
					done
				echo "charpartition mine = "$charpart";" >> {output.nexus}
				echo "end;" >> {output.nexus} #concat.nex

			else # modeltest was not run.
				if [[ $"{params.iqtree_params}" == *"-m "* ]]; then
					model=$(echo "{params.iqtree_params}" | awk '{{ for(i=1;i<=NF;i++) {{ if($i=="-m") {{ if(i < NF) {{ print $(i+1); }}break; }}}}}}')
					echo "$(date) - 'phylociraptor modeltest' was not run before for {wildcards.aligner} and {wildcards.alitrim}. Will run iqtree with the model specified in the config file for each partition/gene: $model" >> {params.wd}/results/statistics/runlog.txt
				else
					echo "$(date) - 'phylociraptor modeltest' was not run before for {wildcards.aligner} and {wildcards.alitrim}. Will run iqtree with -m MFP." >> {params.wd}/results/statistics/runlog.txt
					model="None"
				fi	 
				echo "#nexus" > {output.nexus}
				echo "begin sets;" >> {output.nexus}
				i=1
				charpart=""
					for gene in $(echo "{params.genes}")
					do
						echo "charset part$i = algn/$gene""_aligned_trimmed.fas:*;" >> {output.nexus}
						charpart="${{charpart}}$model:part$i,"
						i=$((++i))
					done

				if [[ "$model" != "None" ]]; then # only write model partition line when a model has been specified.
					echo "charpartition mine = "$charpart";" >> {output.nexus}
				fi
					
				echo "end;" >> {output.nexus} #concat.nex
				
			fi
				# insert code for specified model here
			echo "$(date) - nexus file for iqtree written." >> {params.wd}/results/statistics/runlog.txt
			"""
rule iqtree:
		input:
			rules.prepare_iqtree.output
		output:
			checkpoint = "results/checkpoints/iqtree_{aligner}_{alitrim}_{bootstrap}.{hash}.done",
			statistics = "results/statistics/mltree/mltree_iqtree_{aligner}_{alitrim}_{bootstrap}_statistics.{hash}.txt"
		benchmark:
			"results/statistics/benchmarks/tree/iqtree_{aligner}_{alitrim}_{bootstrap}.{hash}.txt"
		singularity:
			containers["iqtree"]
		params:
			wd = os.getcwd(),
			models = "results/modeltest/best_models_{aligner}_{alitrim}.{hash}.txt",
			nt = "AUTO",
			bb = config["mltree"]["bootstrap"]["iqtree"],
			additional_params = config["mltree"]["options"]["iqtree"],
			seed=config["seed"],
#			bootstrap_cutoff_file = "results/statistics/genetree_filter_{aligner}_{alitrim}.txt",
			genes = get_input_genes
		threads:
			config["mltree"]["threads"]["iqtree"]
		log:
			"log/mltree/iqtree/iqtree-{aligner}-{alitrim}-{bootstrap}.{hash}.txt"
		shell:
			"""
			cd results/phylogeny/iqtree/bootstrap-cutoff-{wildcards.bootstrap}/{wildcards.aligner}-{wildcards.alitrim}.{wildcards.hash}/
			if grep -q "charpartition" concat.nex; then
				# file contains model info
				iqtree --redo -p concat.nex --prefix concat -bb {params.bb} -nt {threads} $(if [[ "{params.seed}" != "None" ]]; then echo "-seed {params.seed}"; fi) {params.additional_params} 2>&1 | tee {params.wd}/{log}
				statistics_string="iqtree\t{wildcards.aligner}\t{wildcards.alitrim}\t{params.bb}\t{wildcards.bootstrap}\t$(ls algn | wc -l)\t$(cat concat.contree)"
			else # partition file does not contain models, also this means no models have been specified in config file.
				echo "Modeltest does not seem to have been run and no model has been specified in config file."
				echo "IQ-Tree will therefore run without the best models and partitioning scheme from phylociraptor modeltest and use -m MFP instead."
				iqtree --redo -p concat.nex -m MFP --prefix concat -bb {params.bb} -nt {threads} $(if [[ "{params.seed}" != "None" ]]; then echo "-seed {params.seed}"; fi) {params.additional_params} 2>&1 | tee {params.wd}/{log}
				statistics_string="iqtree\t{wildcards.aligner}\t{wildcards.alitrim}\t{params.bb}\t{wildcards.bootstrap}\t$(ls algn | wc -l)\t$(cat concat.contree)"
			fi
			echo -e $statistics_string > {params.wd}/{output.statistics}	
			touch {params.wd}/{output.checkpoint}
			"""

rule prepare_iqtree_unpartitioned:
		input:
			get_modeltest_checkpoint,
			get_iqtree_params,
			algn = rules.concatenate.output.alignment
		output:
			algn = directory("results/phylogeny/iqtree-unpartitioned/bootstrap-cutoff-{bootstrap}/{aligner}-{alitrim}.{hash}/algn"),
			concatalgn = "results/phylogeny/iqtree-unpartitioned/bootstrap-cutoff-{bootstrap}/{aligner}-{alitrim}.{hash}/concat.fas",
		singularity:
			containers["iqtree"]
		params:
			wd = os.getcwd(),
			models = get_best_models,
			modelsfile = get_best_models_filename,
			alidir = get_alignment_dir, 
			genes = get_input_genes
		shell:
			"""
			mkdir {output.algn}
			echo "$(date) - prepare_iqtree_unpartitioned {wildcards.aligner}-{wildcards.alitrim}: Will use bootstrap cutoff ({wildcards.bootstrap})" >> {params.wd}/results/statistics/runlog.txt
			echo "$(date) - prepare_iqtree_unpartitioned: Will not consider partitioning scheme for {wildcards.aligner} and {wildcards.alitrim}." >> {params.wd}/results/statistics/runlog.txt
			cp {input.algn} {output.concatalgn}
			# insert code for specified model here
			echo "$(date) - concat.fas file for iqtree-unpartitioned written." >> {params.wd}/results/statistics/runlog.txt
			"""

rule iqtree_unpartitioned:
		input:
			rules.prepare_iqtree_unpartitioned.output
		output:
			checkpoint = "results/checkpoints/iqtree-unpartitioned_{aligner}_{alitrim}_{bootstrap}.{hash}.done",
			statistics = "results/statistics/mltree/mltree_iqtree-unpartitioned_{aligner}_{alitrim}_{bootstrap}_statistics.{hash}.txt"
		benchmark:
			"results/statistics/benchmarks/tree/iqtreei-unpartitioned_{aligner}_{alitrim}_{bootstrap}.{hash}.txt"
		singularity:
			containers["iqtree"]
		params:
			wd = os.getcwd(),
			nt = "AUTO",
			bb = config["mltree"]["bootstrap"]["iqtree-unpartitioned"],
			additional_params = config["mltree"]["options"]["iqtree-unpartitioned"],
			seed=config["seed"],
			genes = get_input_genes
		threads:
			config["mltree"]["threads"]["iqtree-unpartitioned"]
		log:
			"log/mltree/iqtree/iqtree-unpartitioned-{aligner}-{alitrim}-{bootstrap}.{hash}.txt"
		shell:
			"""
			cd results/phylogeny/iqtree-unpartitioned/bootstrap-cutoff-{wildcards.bootstrap}/{wildcards.aligner}-{wildcards.alitrim}.{wildcards.hash}/

			echo "This analysis will be unpartitioned. IQ-Tree will run without the best models and partitioning scheme from phylociraptor modeltest."
			echo "Instead this will use modelsetting specified in config file: {params.additional_params}."
			if [[ $"{params.additional_params}" == *"-m "* ]]; then
				echo "Will use spcified model from config file"
				iqtree -s concat.fas {params.additional_params} --prefix concat -bb {params.bb} -nt {threads} --redo $(if [[ "{params.seed}" != "None" ]]; then echo "-seed {params.seed}"; fi) 2>&1 | tee {params.wd}/{log}
			else
				echo "No model has been specified in config file. Will use default: -m MFP"
				iqtree -s concat.fas --prefix concat -bb {params.bb} -nt {threads} $(if [[ "{params.seed}" != "None" ]]; then echo "-seed {params.seed}"; fi) -m MFP {params.additional_params} 2>&1 | tee {params.wd}/{log}
			fi
			statistics_string="iqtree\t{wildcards.aligner}\t{wildcards.alitrim}\t{params.bb}\t{wildcards.bootstrap}\t$(ls algn | wc -l)\t$(cat concat.contree)"
		
			echo -e $statistics_string > {params.wd}/{output.statistics}	
			touch {params.wd}/{output.checkpoint}
			"""

def pull(wildcards):
	lis = []
	bscuts = allbscuts
	if "NOMODELTEST" in os.environ.keys(): #modeltest was not run!
		bscuts = [0]
	for i in tree_methods:
		for m in config['modeltest']['method']:
			for a in aligners:
				for t in trimmers:
					for b in bscuts:
						if "NOMODELTEST" in os.environ.keys(): # modeltest was not run!
							lis.append("results/checkpoints/"+i+"_"+a+"_"+t+"_"+str(b)+"."+tinference_hashes[str(b)][i][m][t][a]+".done")
							lis.append("results/phylogeny/" + i + "/bootstrap-cutoff-"+str(b)+"/parameters.mltree."+i+"-"+a+"-"+t+"."+tinference_hashes[str(b)][i][m][t][a]+".yaml")
						else: # modeltest was run.
							lis.append("results/checkpoints/"+i+"_"+a+"_"+t+"_"+str(b)+"."+tinference_hashes[str(b)][i][m][t][a]+".done")
							lis.append("results/phylogeny/" + i + "/bootstrap-cutoff-"+str(b)+"/parameters.mltree."+i+"-"+a+"-"+t+"."+tinference_hashes[str(b)][i][m][t][a]+".yaml")
	return lis
	

rule mltree:
	input:
#		expand("results/checkpoints/{treeinfer}_{aligner}_{alitrim}_{bootstrap}.done", aligner=aligners, alitrim=trimmers, treeinfer=tree_methods, bootstrap=bscuts)
		pull,
		rules.read_params_global.output
	output:
		"results/checkpoints/modes/mltree."+current_hash+".done"
	shell:
		"""
		touch {output}
		echo "$(date) - phylociraptor mltree (tree) done." >> results/statistics/runlog.txt
		"""
