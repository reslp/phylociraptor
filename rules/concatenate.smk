#include: "functions.smk"
import yaml

def get_alignment_dir(wildcards):
	return "results/alignments/filtered/"+wildcards.aligner+"-"+wildcards.alitrim+"."+trimmer_hashes[wildcards.alitrim][wildcards.aligner]

def get_best_models(wildcards):
	return "results/modeltest/best_models_"+wildcards.aligner+"_"+wildcards.alitrim+"."+modeltest_hashes["iqtree"][wildcards.alitrim][wildcards.aligner]+".txt"

def get_genetree_stats(wildcards):
	return "results/modeltest/genetree_filter_"+wildcards.aligner+"_"+wildcards.alitrim+"."+modeltest_hashes["iqtree"][wildcards.alitrim][wildcards.aligner]+".txt"

def previous_params_per(wildcards):
	return "results/modeltest/parameters.modeltest."+wildcards.aligner+"-"+wildcards.alitrim+"."+modeltest_hashes["iqtree"][wildcards.alitrim][wildcards.aligner]+".yaml"

def previous_params_global(wildcards):
	return "results/modeltest/parameters.modeltest."+previous_hash+".yaml"

def get_concatenate_params(wildcards):
	if os.environ["MODE"] == "njtree":
		return "results/phylogeny/quicktree/bootstrap-cutoff-"+wildcards.bootstrap+"/parameters.njtree.quicktree-"+wildcards.aligner+"-"+wildcards.alitrim+"."+wildcards.hash+".yaml"
	elif os.environ["MODE"] == "mltree":
		return "results/phylogeny/raxml/bootstrap-cutoff-"+wildcards.bootstrap+"/parameters.mltree.raxml-"+wildcards.aligner+"-"+wildcards.alitrim+"."+wildcards.hash+".yaml"

rule concatenate:
	input:
		checkpoint = get_modeltest_checkpoint,
		params = get_concatenate_params
	output:
		alignment = "results/phylogeny/concatenate/bootstrap-cutoff-{bootstrap}/{aligner}-{alitrim}.{hash}/concat.fas",
		phylip_alignment = "results/phylogeny/concatenate/bootstrap-cutoff-{bootstrap}/{aligner}-{alitrim}.{hash}/concat.phy",
		stockholm_alignment = "results/phylogeny/concatenate/bootstrap-cutoff-{bootstrap}/{aligner}-{alitrim}.{hash}/concat.sto",
		statistics = "results/phylogeny/concatenate/bootstrap-cutoff-{bootstrap}/{aligner}-{alitrim}.{hash}/statistics.txt"
	benchmark:
		"results/statistics/benchmarks/tree-{bootstrap}/concatenate_{aligner}_{alitrim}.{hash}txt"
	params:
		wd = os.getcwd(),
		ids = config["species"],
		models = get_best_models,
		target_dir = "results/phylogeny/concatenate/bootstrap-cutoff-{bootstrap}/{aligner}-{alitrim}.{hash}",
		alidir = get_alignment_dir,
		genes = get_input_genes
	singularity:
		containers["concatnew"]	
	log:
		"log/mltree/concatenate-{aligner}-{alitrim}-{bootstrap}.{hash}.txt"
	shell:
		"""
		echo "$(date) - concatenate {wildcards.aligner}-{wildcards.alitrim}: Will use bootstrap cutoff {wildcards.bootstrap} before creating concatenated alignment" >> results/statistics/runlog.txt
		mkdir -p {params.target_dir}/algn
		for gene in $(echo "{params.genes}") 
		do
			cp {params.wd}/{params.alidir}/"$gene"_aligned_trimmed.fas {params.target_dir}/algn
		done
		# Now run concat:
		concat.py -d {params.target_dir}/algn --runmode concat -o {params.target_dir} --biopython --statistics 2>&1 | tee {log}
		
		# Now convert alignment to phylip:
		python -c "from Bio import AlignIO; alignment = AlignIO.read(open('{output.alignment}'), 'fasta'); print(' '+str(len(alignment))+' '+str(alignment.get_alignment_length())); print(''.join([str(seq.id)+'   '+str(seq.seq)+'\\n' for seq in alignment]));" > {output.phylip_alignment}
		# and stockholm (pfam) format:
		python -c "from Bio import AlignIO; alignment = AlignIO.read(open('{output.alignment}'), 'fasta'); print(alignment.format('stockholm'));" > {output.stockholm_alignment}	
		"""
	
