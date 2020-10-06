
rule modeltest:
	input:
		rules.part2.output
	output:
		models = "results/modeltest/best_models.txt",
		checkpoint = "results/checkpoints/modeltest.done"
	params:
		wd = os.getcwd()
	singularity:
		"docker://reslp/iqtree:2.0rc2"
	threads: 16
	shell:
		"""
		mkdir -p results/modeltest
		echo "name\tmodel\tnseq\tnsites\tconstsites\tinvsites\tparssites\tdistpatt" > {params.wd}/results/modeltest/best_models.txt
		for file in $(ls results/filtered_alignments/*.fas);
		do
			outname=$(basename $file)
			iqtree -m TESTONLY -s $file -msub nuclear -redo --prefix {params.wd}/results/modeltest/$outname -T 16
			printf "$outname\t" >> {params.wd}/results/modeltest/best_models.txt
			cat {params.wd}/results/modeltest/$outname.log | grep "Best-fit model:" | awk -F ":" '{{print $2}}' | awk -F " " '{{printf ("%s\t", $1)}}' >> {params.wd}/results/modeltest/best_models.txt
			cat {params.wd}/results/modeltest/$outname.iqtree | grep "Input data" | awk -F ":" '{{print $2}}' | awk -F " " '{{printf ("%s\t", $1)}}' >> {params.wd}/results/modeltest/best_models.txt
			cat {params.wd}/results/modeltest/$outname.iqtree | grep "Input data" | awk -F ":" '{{print $2}}' | awk -F " " '{{printf ("%s\t", $4)}}' >> {params.wd}/results/modeltest/best_models.txt
			cat {params.wd}/results/modeltest/$outname.iqtree | grep "Number of constant sites" | awk -F ":" '{{print $2}}' | awk -F " " '{{printf ("%s\t", $1)}}' >> {params.wd}/results/modeltest/best_models.txt
			cat {params.wd}/results/modeltest/$outname.iqtree | grep "Number of invariant" | awk -F ":" '{{print $2}}' | awk -F " " '{{printf ("%s\t", $1)}}' >> {params.wd}/results/modeltest/best_models.txt
			cat {params.wd}/results/modeltest/$outname.iqtree | grep "Number of parsimony" | awk -F ":" '{{print $2}}' | awk -F " " '{{printf ("%s\t", $1)}}' >> {params.wd}/results/modeltest/best_models.txt
			cat {params.wd}/results/modeltest/$outname.iqtree | grep "Number of distinct" | awk -F ":" '{{print $2}}' | awk -F " " '{{print $1}}' >> {params.wd}/results/modeltest/best_models.txt
			# currently iqtree output will be removed and only the best model is saved.
			rm {params.wd}/results/modeltest/$outname.*
		done
		touch {output.checkpoint}
		"""

if config["phylogeny"]["concat"] == "yes":
	rule concatenate:
		input:
			checkpoint = rules.part2.output,
			models = rules.modeltest.output.models
		output:
			checkpoint = "results/checkpoints/concatenate.done",
			alignment = "results/phylogeny/concat.fas",
			partitions = "results/phylogeny/partitions.txt",
			phylip_alignment = "results/phylogeny/concat.phy"
		params:
			wd = os.getcwd(),
			ids = config["species"]
		singularity:
			"docker://reslp/concat:0.2"
		shell:
			"""
			tail -n +2 {params.ids} | awk -F "," '{{print $1;}}' | sed 's/ /_/g' > results/phylogeny/ids.txt	
			concat.py -d results/filtered_alignments/ -t results/phylogeny/ids.txt --runmode concat -o results/phylogeny/ --biopython --statistics
			awk 'FNR==NR{{a[$1]=$2;next}}{{print $0"\\t"a[$1]}}' {input.models} results/phylogeny/statistics.txt | awk -F"\\t" 'NR>1{{split($1,b,"_"); print $5", " b[1]"="$2"-"$3}}' > results/phylogeny/partitions_unformated.txt
			# correct some model names to make them raxml compatible:
			# it is not quite clear which models are compatible. During more extensive tests additional problems should show up
			cat results/phylogeny/partitions_unformated.txt | sed 's/JTTDCMut/JTT-DCMut/' > {output.partitions}
			touch {output.checkpoint}
			python -c "from Bio import AlignIO; alignment = AlignIO.read(open('{output.alignment}'), 'fasta'); print(' '+str(len(alignment))+' '+str(alignment.get_alignment_length())); print(''.join([str(seq.id)+'   '+str(seq.seq)+'\\n' for seq in alignment]));" > {output.phylip_alignment}
			"""	
		
	if "iqtree" in config["phylogeny"]["method"]:
		rule iqtree:
			input:
				rules.part2.output
			output:
				checkpoint = "results/checkpoints/iqtree.done"
			singularity:
				"docker://reslp/iqtree:2.0rc2"
			params:
				wd = os.getcwd(),
				nt = "AUTO",
				bb = config["iqtree"]["bootstrap"],
				m = config["iqtree"]["model"]
			threads:
				config["iqtree"]["threads"]
			shell:
				"""
				rm -rf results/phylogeny/concatenated/algn
				mkdir -p results/phylogeny/concatenated/
				cd results/phylogeny/concatenated/
				mkdir algn
				cp {params.wd}/results/filtered_alignments/*.fas algn
				iqtree -p algn/ --prefix concat -bb {params.bb} -nt {params.nt} -m {params.m} -redo -T {threads}
				rm -r algn
				cd {params.wd}
				touch {output.checkpoint}
				"""
	else:  #checkpoint files need to be created anyway
        	rule iqtree:
                	output:
                        	checkpoint = "results/checkpoints/iqtree.done"
                	shell:
                        	"""
                        	touch {output.checkpoint}
				"""
	
	if "raxml" in config["phylogeny"]["method"]:
		rule raxmlng:
			input:
				alignment = rules.concatenate.output.alignment,
				partitions = rules.concatenate.output.partitions
			output:
				checkpoint = "results/checkpoints/raxmlng.done",
				alignment = "results/phylogeny/raxmlng/concat.fas",
				partitions = "results/phylogeny/raxmlng/partitions.txt"
			singularity:
				"docker://reslp/raxml-ng:1.0.0"
			params:
				threads = config["raxmlng"]["threads"],
				bs = config["raxmlng"]["bootstrap"],
				wd = os.getcwd(),
			shell:
				"""
				cp {input.alignment} {output.alignment}
				cp {input.partitions} {output.partitions}
				cd results/phylogeny/raxmlng
				raxml-ng --msa concat.fas --prefix raxmlng -threads {params.threads} --bs-trees {params.bs} --model partitions.txt --all
				cd {params.wd}
				touch {output.checkpoint}
				"""
	else:
		rule raxmlng:
			output: 
				checkpoint = "results/checkpoints/raxmlng.done"
			shell:
				"""
				touch {output.checkpoint}
				"""
	
	if "phylobayes" in config["phylogeny"]["method"]:
		rule phylobayes:
			input:
				alignment = rules.concatenate.output.phylip_alignment,
			output:
				checkpoint = "results/checkpoints/phylobayes.done",
				alignment = "results/phylogeny/phylobayes/concat.phy",
			singularity:
				"docker://reslp/phylobayes-mpi:1.8b"
			params:
				model = config["phylobayes"]["model"],
				threads = config["phylobayes"]["threads"],
				wd = os.getcwd()
			shell:
				"""
				cp {input.alignment} {output.alignment}
				cd results/phylogeny/phylobayes
				unset PE_HOSTFILE
				mpirun -np {params.threads} pb_mpi -d concat.phy {params.model} phylobayes
				cd {params.wd}
				touch {output.checkpoint}
				"""
	else:
		rule phylobayes:
			output:
				checkpoint = "results/checkpoints/phylobayes.done"
			shell:
				"""
				touch {output.checkpoint}
				"""

if config["phylogeny"]["species_tree"] == "yes":
	rule iqtree_gene_trees:
		input:
			rules.part2.output,
		output:
			checkpoint = "results/checkpoints/iqtree_gene_trees.done",
			trees = "results/phylogeny/gene_trees/loci.treefile"
		params:
			wd = os.getcwd(),
		threads:
			config["iqtree"]["threads"]
		singularity:
			"docker://reslp/iqtree:2.0rc2"
		shell:
			"""
			rm -rf results/phylogeny/gene_trees/algn
			mkdir -p results/phylogeny/gene_trees/algn
			cd results/phylogeny/gene_trees
			cp {params.wd}/results/filtered_alignments/*.fas algn
			iqtree -S algn/ --prefix loci -nt AUTO -m MFP -redo -T {threads}
			rm -r algn
			cd {params.wd}
			touch {output.checkpoint}
			"""
	
	rule astral_species_tree:
		input:
			trees = rules.iqtree_gene_trees.output.trees,
			checkpoint = rules.iqtree_gene_trees.output.checkpoint
		output:
			species_tree = "results/phylogeny/astral/species_tree.tre",
			checkpoint = "results/checkpoints/astral_species_tree.done"
		params:
			wd = os.getcwd()
		singularity:
			"docker://reslp/astral:5.7.1"
		shell:
			"""
			java -jar /ASTRAL-5.7.1/Astral/astral.5.7.1.jar -i {input.trees} -o {output.species_tree}
			touch {output.checkpoint}
			"""
			
else: #checkpoint files need to be created anyway
	rule iqtree_gene_trees:
		output:
			checkpoint = "results/checkpoints/iqtree_gene_trees.done"
		shell:
			"""
			touch {output.checkpoint}
			"""
	rule astral_species_tree:
		output:
			checkpoint = "results/checkpoints/astral_species_tree.done"
		shell:
			"""
			touch {output.checkpoint}
			"""

