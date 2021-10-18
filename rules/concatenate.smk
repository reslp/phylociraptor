rule concatenate:
	input:
#		checkpoint = "results/checkpoints/filter_alignments_{alitrim}_{aligner}.done"
		checkpoint = "results/statistics/filter-{aligner}-{alitrim}/alignment_filter_information_{alitrim}_{aligner}.txt"
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
		models = "results/modeltest/best_models_{aligner}_{alitrim}.txt",
		bootstrap_cutoff_file = "results/statistics/genetree_filter_{aligner}_{alitrim}.txt"
	singularity:
		"docker://reslp/concat:0.3"
	shell:
		"""
		if [[ -f {params.wd}/{params.bootstrap_cutoff_file} ]]; then #in case gene trees were filtered with bootstrap cutoff before
			echo "$(date) - concatenate {wildcards.aligner}-{wildcards.alitrim}: Will use bootstrap cutoff before creating concatenated alignment" >> results/statistics/runlog.txt
			mkdir -p results/phylogeny/concatenate/{wildcards.aligner}-{wildcards.alitrim}/algn
			for gene in $(cat {params.wd}/{params.bootstrap_cutoff_file} | awk -F"\t" '{{if ($8 == "OK") {{print $1;}}}}');
			do
				cp {params.wd}/results/alignments/filtered/{wildcards.aligner}-{wildcards.alitrim}/"$gene"_aligned_trimmed.fas results/phylogeny/concatenate/{wildcards.aligner}-{wildcards.alitrim}/algn
			done
			concat.py -d results/phylogeny/concatenate/{wildcards.aligner}-{wildcards.alitrim}/algn --runmode concat -o results/phylogeny/concatenate/{wildcards.aligner}-{wildcards.alitrim} --biopython --statistics
		else	#take all alignments
			concat.py -d results/alignments/filtered/{wildcards.aligner}-{wildcards.alitrim} --runmode concat -o results/phylogeny/concatenate/{wildcards.aligner}-{wildcards.alitrim} --biopython --statistics
		fi
		
		# Now convert alignment to phylip:
		python -c "from Bio import AlignIO; alignment = AlignIO.read(open('{output.alignment}'), 'fasta'); print(' '+str(len(alignment))+' '+str(alignment.get_alignment_length())); print(''.join([str(seq.id)+'   '+str(seq.seq)+'\\n' for seq in alignment]));" > {output.phylip_alignment}
		# and stockholm (pfam) format:
		python -c "from Bio import AlignIO; alignment = AlignIO.read(open('{output.alignment}'), 'fasta'); print(alignment.format('stockholm'));" > {output.stockholm_alignment}	
		"""
	
