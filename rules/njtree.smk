include: "functions.smk"
include: "concatenate.smk"

os.environ["MODE"] = "njtree" #for communication with concatenate.smk

aligners = get_aligners()		
trimmers = get_trimmers()		
bscuts = get_bootstrap_cutoffs()

#create new hashes for current stage 
hashes = collect_hashes("njtree", config, configfi, wd=os.getcwd())

filter_orthology_hash = hashes['filter-orthology']["global"]
aligner_hashes = hashes['align']["per"]
trimmer_hashes = hashes['filter-align']["per"]
filter_align_hash = hashes['filter-align']['global']
modeltest_hashes = hashes['modeltest']["per"]
tinference_hashes = hashes['njtree']["per"]
current_hash = hashes['njtree']["global"]
previous_hash = hashes['modeltest']["global"]

########### functions specifically for this step
def compare_njtree(wildcards):
	return [trigger("results/phylogeny/quicktree/bootstrap-cutoff-{bootstrap}/parameters.njtree.{inference}-{aligner}-{alitrim}.{hash}.yaml".format(bootstrap=wildcards.bootstrap, inference=wildcards.inference, aligner=wildcards.aligner, alitrim=wildcards.alitrim, hash=wildcards.hash), configfi)]

def pull(wildcards):
	lis = []
	for i in config["njtree"]["method"]:
		for m in config['modeltest']['method']:
			for a in aligners:
				for t in trimmers:
					for b in bscuts:
						lis.append("results/checkpoints/"+i+"_"+a+"_"+t+"_"+str(b)+"."+tinference_hashes[str(b)][i][m][t][a]+".done")
						lis.append("results/phylogeny/" + i + "/bootstrap-cutoff-"+str(b)+"/parameters.njtree."+i+"-"+a+"-"+t+"."+tinference_hashes[str(b)][i][m][t][a]+".yaml")
	return lis
########### 

rule read_params_per:
	input:
		trigger = compare_njtree,
		previous = previous_params_per #from concatenate.smk
	output:
		"results/phylogeny/quicktree/bootstrap-cutoff-{bootstrap}/parameters.njtree.{inference}-{aligner}-{alitrim}.{hash}.yaml"
	shell:
		"""
		bin/read_write_yaml.py {input.trigger} {output} genetree_filtering,bootstrap_cutoff,{wildcards.bootstrap} njtree,options,{wildcards.inference}
		cat {input.previous} >> {output}
		"""

rule read_params_global:
	input:
		trigger = compare("results/parameters.njtree."+current_hash+".yaml", configfi),
		previous = "results/modeltest/parameters.modeltest."+previous_hash+".yaml"
	output:
		"results/phylogeny/quicktree/parameters.njtree."+current_hash+".yaml"
	shell:
		"""
		bin/read_write_yaml.py {input.trigger} {output} genetree_filtering,bootstrap_cutoff njtree,method njtree,options
		cat {input.previous} >> {output}
		"""

rule quicktree:
	input:
		rules.concatenate.output.stockholm_alignment
	output:
		tree = "results/phylogeny/quicktree/bootstrap-cutoff-{bootstrap}/{aligner}-{alitrim}.{hash}/njtree.tre",
		checkpoint = "results/checkpoints/quicktree_{aligner}_{alitrim}_{bootstrap}.{hash}.done"
	params:
		additional_params = config["njtree"]["options"]["quicktree"],
	singularity: containers["quicktree"] 
	shell:
		"""
		quicktree -in a $(if [[ "{params.additional_params}" != "None" ]]; then echo "{params.additional_params}"; fi) {input} > {output.tree}
		touch {output.checkpoint}
		"""

		
rule njtree:
	input:
		pull,
		rules.read_params_global.output
	output:
		"results/checkpoints/modes/njtree."+current_hash+".done"
	shell:
		"""
		echo "$(date) - Pipeline part njtre done." >> results/statistics/runlog.txt
		touch {output}
		"""
