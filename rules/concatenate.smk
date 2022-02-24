def get_input_genes(wildcards):
	bs_cutoff = int(wildcards.bootstrap)
	list_of_genes = []
	with open("results/modeltest/genetree_filter_" + wildcards.aligner + "_" + wildcards.alitrim + ".txt") as file:
		for line in file:
			gene = line.split("\t")[0]
			bs_value = int(line.strip().split("\t")[-1])
			if bs_value >= bs_cutoff:
				list_of_genes.append(gene)
	return list_of_genes		

rule concatenate:
	input:
		checkpoint = "results/checkpoints/modeltest/aggregate_best_models_{aligner}_{alitrim}.done"
#		checkpoint = "results/checkpoints/filter_alignments_{alitrim}_{aligner}.done"
#		checkpoint = "results/statistics/filter-{aligner}-{alitrim}/alignment_filter_information_{alitrim}_{aligner}.txt"
#		checkpoint = "results/checkpoints/modes/modeltest.done"
	output:
		alignment = "results/phylogeny-{bootstrap}/concatenate/{aligner}-{alitrim}/concat.fas",
		phylip_alignment = "results/phylogeny-{bootstrap}/concatenate/{aligner}-{alitrim}/concat.phy",
		stockholm_alignment = "results/phylogeny-{bootstrap}/concatenate/{aligner}-{alitrim}/concat.sto",
		statistics = "results/phylogeny-{bootstrap}/concatenate/{aligner}-{alitrim}/statistics.txt"
	benchmark:
		"results/statistics/benchmarks/tree-{bootstrap}/concatenate_{aligner}_{alitrim}.txt"
	params:
		wd = os.getcwd(),
		ids = config["species"],
		datatype = config["filtering"]["seq_type"],
		models = "results/modeltest/best_models_{aligner}_{alitrim}.txt",
#		bootstrap_cutoff_file = "results/statistics/genetree_filter_{aligner}_{alitrim}.txt",
		genes = get_input_genes
	singularity:
		containers["concatnew"]	
	shell:
		"""
		echo "$(date) - concatenate {wildcards.aligner}-{wildcards.alitrim}: Will use bootstrap cutoff {wildcards.bootstrap} before creating concatenated alignment" >> results/statistics/runlog.txt
		mkdir -p results/phylogeny-{wildcards.bootstrap}/concatenate/{wildcards.aligner}-{wildcards.alitrim}/algn
		for gene in $(echo "{params.genes}") 
		do
			cp {params.wd}/results/alignments/filtered/{wildcards.aligner}-{wildcards.alitrim}/"$gene"_aligned_trimmed.fas results/phylogeny-{wildcards.bootstrap}/concatenate/{wildcards.aligner}-{wildcards.alitrim}/algn
		done
		# Now run concat:
		concat.py -d results/phylogeny-{wildcards.bootstrap}/concatenate/{wildcards.aligner}-{wildcards.alitrim}/algn --runmode concat -o results/phylogeny-{wildcards.bootstrap}/concatenate/{wildcards.aligner}-{wildcards.alitrim} --biopython --statistics
		
		# Now convert alignment to phylip:
		python -c "from Bio import AlignIO; alignment = AlignIO.read(open('{output.alignment}'), 'fasta'); print(' '+str(len(alignment))+' '+str(alignment.get_alignment_length())); print(''.join([str(seq.id)+'   '+str(seq.seq)+'\\n' for seq in alignment]));" > {output.phylip_alignment}
		# and stockholm (pfam) format:
		python -c "from Bio import AlignIO; alignment = AlignIO.read(open('{output.alignment}'), 'fasta'); print(alignment.format('stockholm'));" > {output.stockholm_alignment}	
		"""
	
