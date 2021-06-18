rule concatenate:
	input:
		checkpoint = "checkpoints/filter_align.done"
	output:
		alignment = "results/phylogeny/concat.fas",
		phylip_alignment = "results/phylogeny/concat.phy",
		stockholm_alignment = "results/phylogeny/concat.sto",
		statistics = "results/phylogeny/statistics.txt"
	benchmark:
		"results/statistics/benchmarks/tree/concatenate.txt"
	params:
		wd = os.getcwd(),
		ids = config["species"],
		datatype = config["filtering"]["seq_type"],
		models = "results/modeltest/best_models.txt"
	singularity:
		"docker://reslp/concat:0.3"
	shell:
		"""
		concat.py -d results/filtered_alignments/ --runmode concat -o results/phylogeny/ --biopython --statistics
		
		# Now convert alignment to phylip:
		python -c "from Bio import AlignIO; alignment = AlignIO.read(open('{output.alignment}'), 'fasta'); print(' '+str(len(alignment))+' '+str(alignment.get_alignment_length())); print(''.join([str(seq.id)+'   '+str(seq.seq)+'\\n' for seq in alignment]));" > {output.phylip_alignment}
		# and stockholm (pfam) format:
		python -c "from Bio import AlignIO; alignment = AlignIO.read(open('{output.alignment}'), 'fasta'); print(alignment.format('stockholm'));" > {output.stockholm_alignment}	
		"""
	
