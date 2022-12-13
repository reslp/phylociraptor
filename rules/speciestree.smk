include: "functions.smk"

ruleorder: read_params_global > read_params_per

aligners = get_aligners()		
trimmers = get_trimmers()		
bscuts = get_bootstrap_cutoffs()

#create new hashes for current stage 
hashes = collect_hashes("speciestree", config, configfi)

filter_orthology_hash = hashes['filter-orthology']["global"]
aligner_hashes = hashes['align']["per"]
trimmer_hashes = hashes['filter-align']["per"]
modeltest_hashes = hashes['modeltest']["per"]
tinference_hashes = hashes['speciestree']["per"]
current_hash = hashes['speciestree']["global"]
previous_hash = hashes['modeltest']["global"]

######################## functions specifically for this step
def previous_params_per(wildcards):
	return "results/modeltest/parameters.modeltest."+wildcards.aligner+"-"+wildcards.alitrim+"."+modeltest_hashes["iqtree"][wildcards.alitrim][wildcards.aligner]+".yaml"

def previous_params_per(wildcards):
	return "results/modeltest/parameters.modeltest."+previous_hash+".yaml"

def compare_speciestree(wildcards):
	return [trigger("results/phylogeny/{inference}/bootstrap-cutoff-{bootstrap}/parameters.speciestree.{inference}-{aligner}-{alitrim}.{hash}.yaml".format(bootstrap=wildcards.bootstrap, inference=wildcards.inference, aligner=wildcards.aligner, alitrim=wildcards.alitrim, hash=wildcards.hash), configfi)]	
########################

BUSCOS, = glob_wildcards("results/orthology/busco/busco_sequences_deduplicated."+filter_orthology_hash+"/{busco}_all.fas")

def return_trees(wildcards):
	lis = []
        for busco in BUSCOS:
		if os.path.isfile("results/alignments/filtered/"+wildcards.aligner+"-"+wildcards.alitrim+"/"+str(busco)+"_aligned_trimmed.fas"):
			lis.append("results/phylogeny/gene_trees/"+wildcards.aligner+"-"+wildcards.alitrim+"/"+busco+"/"+busco+"_gt.treefile")
	return lis

rule read_params_per:
	input:
		trigger = compare_speciestree,
		previous = previous_params_per
	params:
		conf = configfi 
	output:
		"results/phylogeny/{inference}/bootstrap-cutoff-{bootstrap}/parameters.speciestree.{inference}-{aligner}-{alitrim}.{hash}.yaml"
	shell:
		"""
		bin/read_write_yaml.py {params.conf} {output} seed genetree_filtering,bootstrap_cutoff,{wildcards.bootstrap} speciestree,options,{wildcards.inference} speciestree,include
		cat {input.previous} >> {output}
		"""

rule read_params_global:
	input:
		trigger = compare("results/parameters.speciestree."+current_hash+".yaml", configfi),
		previous = "results/modeltest/parameters.modeltest."+previous_hash+".yaml"
	output:
		"results/parameters.speciestree."+current_hash+".yaml"
	shell:
		"""
		bin/read_write_yaml.py {input.trigger} {output} seed genetree_filtering,bootstrap_cutoff speciestree,method speciestree,options speciestree,include
		cat {input.previous} >> {output}
		"""

def get_modeltest_hash(wildcards):
	return modeltest_hashes["iqtree"][wildcards.alitrim][wildcards.aligner]

rule aggregate_gene_trees:
	input:
#		treefiles = return_trees,
		"results/phylogeny/astral/bootstrap-cutoff-{bootstrap}/parameters.speciestree.astral-{aligner}-{alitrim}.{hash}.yaml",
		get_modeltest_checkpoint,
	output:
		trees = "results/phylogeny/astral/bootstrap-cutoff-{bootstrap}/{aligner}-{alitrim}.{hash}/trees_{aligner}_{alitrim}.tre",
		checkpoint = "results/checkpoints/aggregate_gene_trees_{aligner}_{alitrim}_{bootstrap}.{hash}.done"
#		genetree_filter_stats = "results/statistics/genetree_filter_{aligner}_{alitrim}_{bootstrap}.{hash}.txt"
	params:
		include = config["speciestree"]["include"],
		wd = os.getcwd(),
		modeltest_hash = get_modeltest_hash,
		genes = get_input_genes
	shell:
		"""
		rm -rf {output.trees}

		if [[ ! -f "{params.include}" ]]
		then
			echo "$(date) - phylociraptor will use gene tree filtering based on average bootstrap support value of {wildcards.bootstrap}." >> {params.wd}/results/statistics/runlog.txt
			for gene in $(echo "{params.genes}")
			do
					cat results/modeltest/{wildcards.aligner}-{wildcards.alitrim}.{params.modeltest_hash}/$gene/*.treefile >> {output.trees}
			done
		else
			cat {params.include} > {output.trees}
		fi
		touch {output.checkpoint}
		"""

rule astral:
	input:
		trees = "results/phylogeny/astral/bootstrap-cutoff-{bootstrap}/{aligner}-{alitrim}.{hash}/trees_{aligner}_{alitrim}.tre",
#		trees = rules.aggregate_gene_trees.output.trees,
#		checkpoint = rules.aggregate_gene_trees.output.checkpoint
	output:
		species_tree = "results/phylogeny/astral/bootstrap-cutoff-{bootstrap}/{aligner}-{alitrim}.{hash}/species_tree.tre",
		statistics = "results/statistics/speciestree/astral_{aligner}-{alitrim}-{bootstrap}_statistics.{hash}.txt",
		checkpoint = "results/checkpoints/astral_{aligner}_{alitrim}_{bootstrap}.{hash}.done"
	benchmark:
		"results/statistics/benchmarks/speciestree/astral_species_tree_{aligner}_{alitrim}_{bootstrap}.{hash}.txt"
	params:
		wd = os.getcwd(),
		seed=config["seed"]
	singularity:
		containers["astral"]	
	log:	"log/speciestree/astral/astral-{aligner}-{alitrim}-{bootstrap}.{hash}.txt"
	shell:
		"""
		java -jar /ASTRAL-5.7.1/Astral/astral.5.7.1.jar -i {input.trees} -o {output.species_tree} $(if [[ "{params.seed}" != "None" ]]; then echo "-s {params.seed}"; fi) 2>&1 | tee {log}
		statistics_string="astral\t{wildcards.aligner}\t{wildcards.alitrim}\t{wildcards.bootstrap}\t$(cat {input.trees} | wc -l)\t$(cat {output.species_tree})"
		echo -e $statistics_string > {params.wd}/{output.statistics}	
		touch {output.checkpoint}
		"""

def pull(wildcards):
	lis = []
	for i in config["speciestree"]["method"]:
		for m in config['modeltest']['method']:
			for a in aligners:
				for t in trimmers:
					for b in bscuts:
						lis.append("results/checkpoints/"+i+"_"+a+"_"+t+"_"+str(b)+"."+tinference_hashes[str(b)][i][m][t][a]+".done")
						lis.append("results/phylogeny/"+i+"/bootstrap-cutoff-"+str(b)+"/parameters.speciestree."+i+"-"+a+"-"+t+"."+tinference_hashes[str(b)][i][m][t][a]+".yaml")
	return lis

rule speciestree:
	input:
		pull,
		rules.read_params_global.output
	output:
		"results/checkpoints/modes/speciestree."+current_hash+".done"
	shell:
		"""
		
		echo "$(date) - phylociraptor speciestree reconstruction done." >> results/statistics/runlog.txt
		touch {output}
		"""

