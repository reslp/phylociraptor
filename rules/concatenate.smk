rule concatenate:
	input:
		checkpoint = "results/checkpoints/filter_alignments_{alitrim}_{aligner}.done"
	output:
		alignment = "results/phylogeny/concatenate/{aligner}-{alitrim}/concat.fas",
		phylip_alignment = "results/phylogeny/concatenate/{aligner}-{alitrim}/concat.phy",
		stockholm_alignment = "results/phylogeny/concatenate/{aligner}-{alitrim}/concat.sto",
		statistics = "results/phylogeny/concatenate/{aligner}-{alitrim}/statistics.txt"
	benchmark:
		"results/statistics/benchmarks/tree/concatenate_{aligner}_{alitrim}.txt"
	params:
		wd = os.getcwd(),
		ids = config["species"],
		datatype = config["filtering"]["seq_type"],
		models = "results/modeltest/best_models_{aligner}_{alitrim}.txt"
	singularity:
		"docker://reslp/concat:0.3"
	shell:
		"""
		concat.py -d results/alignments/filtered/{wildcards.aligner}-{wildcards.alitrim} --runmode concat -o results/phylogeny/concatenate/{wildcards.aligner}-{wildcards.alitrim} --biopython --statistics
		
		# Now convert alignment to phylip:
		python -c "from Bio import AlignIO; alignment = AlignIO.read(open('{output.alignment}'), 'fasta'); print(' '+str(len(alignment))+' '+str(alignment.get_alignment_length())); print(''.join([str(seq.id)+'   '+str(seq.seq)+'\\n' for seq in alignment]));" > {output.phylip_alignment}
		# and stockholm (pfam) format:
		python -c "from Bio import AlignIO; alignment = AlignIO.read(open('{output.alignment}'), 'fasta'); print(alignment.format('stockholm'));" > {output.stockholm_alignment}	
		"""
	
