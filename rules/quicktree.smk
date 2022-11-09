include: "functions.smk"
import yaml

# get list of containers to use:
with open("data/containers.yaml", "r") as yaml_stream:
    containers = yaml.safe_load(yaml_stream)

aligners = get_aligners()		
trimmers = get_trimmers()		
bscuts = get_bootstrap_cutoffs()

#create new hashes for current stage 
hashes = collect_hashes("njtree")

filter_orthology_hash = hashes['filter-orthology']
aligner_hashes = hashes['align']["per"]
trimmer_hashes = hashes['filter-align']["per"]
modeltest_hashes = hashes['modeltest']["per"]
tinference_hashes = hashes['tree_inference']["per"]
print("TEST",tinference_hashes)
current_hash = hashes['tree_inference']["global"]
previous_hash = hashes['modeltest']["global"]

def compare_njtree(wildcards):
	return [trigger("results/phylogeny-{bootstrap}/parameters.njtree.{inference}-{aligner}-{alitrim}.{hash}.yaml".format(bootstrap=wildcards.bootstrap, inference=wildcards.inference, aligner=wildcards.aligner, alitrim=wildcards.alitrim, hash=wildcards.hash), configfi)]

def get_concatenate_params(wildcards):
	return "results/phylogeny-"+wildcards.bootstrap+"/parameters.njtree.quicktree-"+wildcards.aligner+"-"+wildcards.alitrim+"."+wildcards.hash+".yaml"

include: "concatenate.smk"

rule read_params_per:
	input:
		trigger = compare_njtree,
		previous = previous_params_per
	output:
		"results/phylogeny-{bootstrap}/parameters.njtree.{inference}-{aligner}-{alitrim}.{hash}.yaml"
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
		"results/parameters.njtree."+current_hash+".yaml"
	shell:
		"""
		bin/read_write_yaml.py {input.trigger} {output} genetree_filtering,bootstrap_cutoff njtree,method njtree,options
		cat {input.previous} >> {output}
		"""
rule quicktree:
	input:
		rules.concatenate.output.stockholm_alignment
	output:
		tree = "results/phylogeny-{bootstrap}/quicktree/{aligner}-{alitrim}.{hash}/njtree.tre",
		checkpoint = "results/checkpoints/quicktree_{aligner}_{alitrim}_{bootstrap}.{hash}.done"
	params:
		additional_params = config["njtree"]["options"]["quicktree"],
	singularity: containers["quicktree"] 
	shell:
		"""
		quicktree -in a $(if [[ "{params.additional_params}" != "None" ]]; then echo "{params.additional_params}"; fi) {input} > {output.tree}
		touch {output.checkpoint}
		"""
def pull(wildcards):
	lis = []
	for i in config["njtree"]["method"]:
		for m in config['modeltest']['method']:
			for a in aligners:
				for t in trimmers:
					for b in bscuts:
						lis.append("results/checkpoints/"+i+"_"+a+"_"+t+"_"+str(b)+"."+tinference_hashes[str(b)][i][m][t][a]+".done")
						lis.append("results/phylogeny-"+str(b)+"/parameters.njtree."+i+"-"+a+"-"+t+"."+tinference_hashes[str(b)][i][m][t][a]+".yaml")
	return lis
		
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
