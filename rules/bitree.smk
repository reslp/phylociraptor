os.environ["MODE"] = "bitree" # for communication with concatenate.smk

include: "functions.smk"
include: "concatenate.smk"


ruleorder: read_params_global > read_params_per
if os.environ["MODELTEST"] == "no":
	print(now(), "INFO: Since modeltesting was not run previously will try to use only bootstrap cutoff 0.")

#create new hashes for current stage
hashes = collect_hashes("bitree", config, configfi, wd=os.getcwd())

aligners = get_aligners()		
trimmers = get_trimmers()		
tree_methods = get_treemethods()
chains = get_bichains()
bscuts = get_bootstrap_cutoffs()

filter_orthology_hash = hashes['filter-orthology']["global"]
aligner_hashes = hashes['align']["per"]
trimmer_hashes = hashes['filter-align']["per"]
modeltest_hashes = hashes['modeltest']["per"]
tinference_hashes = hashes['bitree']["per"]
current_hash = hashes['bitree']["global"]
previous_hash = hashes['modeltest']["global"]
previous_alitrim_hash = hashes['filter-align']["global"]


############ functions specifically for this step
def compare_bitree(wildcards):
	return [trigger("results/phylogeny/{inference}/bootstrap-cutoff-{bootstrap}/parameters.bitree.{inference}-{aligner}-{alitrim}.{hash}.yaml".format(bootstrap=wildcards.bootstrap, inference=wildcards.inference, aligner=wildcards.aligner, alitrim=wildcards.alitrim, hash=wildcards.hash), configfi)]	
############

rule read_params_per:
	input:
		trigger = compare_bitree,
		previous = previous_params_per
	output:
		"results/phylogeny/{inference}/bootstrap-cutoff-{bootstrap}/parameters.bitree.{inference}-{aligner}-{alitrim}.{hash}.yaml"
	shell:
		"""
		bin/read_write_yaml.py {input.trigger} {output} seed genetree_filtering,bootstrap_cutoff,{wildcards.bootstrap} bitree,options,{wildcards.inference} bitree,chains
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

rule phylobayes:
	input:
		alignment = rules.concatenate.output.phylip_alignment
	output:
		checkpoint = "results/checkpoints/phylobayes_{aligner}_{alitrim}_{bootstrap}_{chain}.{hash}.done",
		internalcheck = "results/phylogeny/phylobayes/bootstrap-cutoff-{bootstrap}/{aligner}-{alitrim}.{hash}/chain-{chain}/internal-check.done",
	singularity:
		containers["phylobayes"]	
	params:
		threads = config["bitree"]["threads"]["phylobayes"],
		options = config["bitree"]["options"]["phylobayes"],
		chains = config["bitree"]["chains"]["phylobayes"],
		generations = config["bitree"]["generations"]["phylobayes"],
		sampling = config["bitree"]["sampling"]["phylobayes"],
		wd = os.getcwd()
	threads:
		config["bitree"]["threads"]["phylobayes"]
	shell:
		"""
		# prepare things
		cp {input.alignment} results/phylogeny/phylobayes/bootstrap-cutoff-{wildcards.bootstrap}/{wildcards.aligner}-{wildcards.alitrim}.{wildcards.hash}/chain-{wildcards.chain}/
		cd results/phylogeny/phylobayes/bootstrap-cutoff-{wildcards.bootstrap}/{wildcards.aligner}-{wildcards.alitrim}.{wildcards.hash}/chain-{wildcards.chain}
		
		# get number of threads per chain:
		#THREADS=$(echo "{params.threads} {params.chains}" | awk '{{print int($1/$2)}}')
		#echo "Each chain of {params.chains} should run with $THREADS threads."
		#for i in {{1..{params.chains}}}
		#do
		CHAINNAME="chain-"{wildcards.chain}
		PARAMS="--oversubscribe"
		if [[ -f "$CHAINNAME.chain" ]]; then
			echo "Output seems to already exist. Will resume chain."
			mpirun -np {params.threads} $PARAMS pb_mpi -x {params.sampling} {params.generations} $CHAINNAME
		else
			echo "Output does not exist. Will start new chain..."
			mpirun -np {params.threads} $PARAMS pb_mpi {params.options} -x {params.sampling} {params.generations} -d $(basename {input.alignment}) $CHAINNAME
		fi
		#done
		touch {output.internalcheck}
		touch {output.checkpoint}
		"""

#	rule merge_phylobayes_chains:
#		input:
#			expand("results/checkpoints/phylobayes_chain{chain}.done", chain = chains)	
#		output:
#			checkpoint = "results/checkpoints/merge_phylobayes_chains.dome"
#		shell:
#			"""
#			touch {output.checkpoint}
#			"""
#else:
#	chains = [str(i) for i in range(1,int(config["phylobayes"]["nchains"])+1)]
#	rule phylobayes:
#		output:
#			checkpoint = "results/checkpoints/phylobayes_chain{chain}.done"
#		shell:
#			"""
#			touch {output.checkpoint}
#			"""
#	rule merge_phylobayes_chains:
#		input:
#			expand("results/checkpoints/phylobayes_chain{chain}.done", chain = chains)
#		output:
#			checkpoint = "results/checkpoints/merge_phylobayes_chains.dome"
#		shell:
#			"""
#			touch {output.checkpoint}
#			"""

def pull(wildcards):
	lis = []
	if "NOMODELTEST" in os.environ.keys():
		bscuts = [0]
	# decide what has been run and modify accordingly
	for i in config["bitree"]["method"]:
		for m in config['modeltest']['method']:
			for a in aligners:
				for t in trimmers:
					for b in bscuts:
						for c in range(1, int(config['bitree']['chains']['phylobayes'])+1):
							lis.append("results/checkpoints/"+i+"_"+a+"_"+t+"_"+str(b)+"_"+str(c)+"."+tinference_hashes[str(b)][i][m][t][a]+".done")
							lis.append("results/phylogeny/" + i + "/bootstrap-cutoff-"+str(b)+"/parameters.bitree."+i+"-"+a+"-"+t+"."+tinference_hashes[str(b)][i][m][t][a]+".yaml")
	print(lis)
	return lis

rule bitree:
	input:
#		expand("results/checkpoints/{treeinfer}_{aligner}_{alitrim}_{bootstrap}_{chain}.done", aligner=aligners, alitrim=trimmers, treeinfer=tree_methods, bootstrap=bscuts, chain=[1,2,3,4])
		pull,
		rules.read_params_global.output
	output:
		"results/checkpoints/modes/bitree."+current_hash+".done"
	shell:
		"""
		touch {output}
		echo "$(date) - phylociraptor bitree (tree) done." >> results/statistics/runlog.txt
		"""
