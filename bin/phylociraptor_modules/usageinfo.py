
def hi():
	print("hi")

#some basic variables needed:
try:
	with open(".version", "r") as file:
		version = file.readline().strip("\n")
except FileNotFoundError:
	version = "unkown"

phylociraptor = """
			     Welcome to
           __          __           _                  __            
    ____  / /_  __  __/ /___  _____(_)________ _____  / /_____  _____
   / __ \/ __ \/ / / / / __ \/ ___/ / ___/ __ `/ __ \/ __/ __ \/ ___/
  / /_/ / / / / /_/ / / /_/ / /__/ / /  / /_/ / /_/ / /_/ /_/ / /    
 / .___/_/ /_/\__, /_/\____/\___/_/_/   \__,_/ .___/\__/\____/_/     
/_/          /____/                         /_/                      

	  the rapid phylogenomic tree calculator, v%s
""" % version

default_help = phylociraptor + """

Usage: phylociraptor <command> <arguments>

Commands:
	setup			Setup pipeline
	orthology		Infer orthologs in a set of genomes
	filter-orthology	Filter orthology results
	align			Create alignments for orthologous genes
	filter-align		Trim and filter alignments
	modeltest		Calculate gene-trees and perform modeltesting
	mltree			Calculate Maximum-Likelihood phylogenomic trees
	speciestree		Calculate species tree
	njtree			Calculate Neighbor-Joining tree

	report			Create a HTML report of the run
	check			Quickly check status of the run
	util			Utilities for a posteriori analyses of trees

	-v, --version 		Print version
	-h, --help		Display help

Examples:
	To see options for the setup step:
	./phylociraptor setup -h

	To run orthology inferrence for a set of genomes on a SLURM cluster:
	./phylociraptor orthology -t slurm -c data/cluster-config-SLURM.yaml

	To filter alignments overwriting the number of parsimony informative sites set in the config file:
	./phylciraptor filter-align --npars_cutoff 12
        
"""
standard_arguments= """
Argumemts:
	-t, --cluster		Specify cluster type. Options: slurm, sge, torque, serial.
	-c, --cluster-config	Specify Cluster config file path. Default: data/cluster-config-CLUSTERTYPE.yaml.template
	-f, --force		Soft force runmode which has already been run.
	-F, --FORCE		Hard force runmode recreating all output.
	
	--dry			Make a dry run.
	--verbose		Display more output.
	-h, --help		Display help.
"""

additional_arguments = """

Additonal customization (optional):
	--singularity=		Pass additional arguments to the singularity containers (eg. additional bindpoints).
				Have to be put under quotes " or '
	--snakemake=		Pass additional arguments to snakemake.
				Have to be put under quotes " or ' 
	--rerun-incomplete	This can be used if the analysis fails with an error indicating incomplete files.
				Will be passed on to snakemake. Equivalent to --snakemake="--rerun-incomplete". """

setup_help = """
phylociraptor setup - will prepare your analysis

Usage: phylociraptor setup <arguments>
""" + standard_arguments + additional_arguments + """
	--config-file		Custom config-file path. (Default: data/config.yaml)	
	--busco_set		BUSCO set to download. (Default: value from config.yaml)
	--samples_csv		Samples CSV file path. (Default: value from config.yaml)
	--add_genomes		Will only add additional genomes specified in th CSV file.
"""

orthology_help = """
phylociraptor orthology - Will infer orthologous genes in a set of genomes.

Usage: phylociraptor orthology <arguments>
""" + standard_arguments + additional_arguments + """
	--config-file           Custom config-file path. (Default: data/config.yaml)
	--busco_threads		Number of threads for each BUSCO run. (Default: value from config.yaml)
	--augustus_species	Pretrained species for Augustus. (Deafult: value from config.yaml)
	--additional_params	Additional parameter passed on to BUSCO. Must be placed inside quotes. (Default: value from config.yaml)
	        
"""

forthology_help = """
phylociraptor filter-orthology - Will filter orthology results produced by phylociraptor orthology.

Usage: phylociraptor filter-orthology <arguments>

""" + standard_arguments + additional_arguments + """
	--config-file           Custom config-file path. (Default: data/config.yaml)
	--dupseq		Set how occasionally found duplicated sequences should be handled.
				Options: persample; remove only samples with duplicated sequences
					 perfiler; remove complete file
				(Default: value from config.yaml)
	--cutoff		Minimum BUSCO completeness for a sample to be kept.
				(Default: value from config.yaml)
	--minsp			Mimimum number of species that need to have a BUSCO gene for it to be kept.
				(Default: value from config.yaml)
	--seq_type		Type of sequence data to use. Options (aa, nu).
				(Default: value from config.yaml)
        
"""

align_help = """
phylociraptor align - Will create alignments for a set of single-copy orthologous genes.

Usage: phylociraptor align <arguments>

""" + standard_arguments + additional_arguments + """
	--config-file           Custom config-file path. (Default: data/config.yaml)
	--method		Alignment method. Options: mafft (Default: mafft; read from config.yaml)
	--parameters        	Commandline arguments for alignment method. (Default: read from config.yaml)
	--threads		Number of threads for alignment step. (Default: read from config.yaml)	

"""

falign_help = """
phylociraptor filter-align - Will filter alignments.

Usage: phylociraptor filter-align <arguments>

""" + standard_arguments + additional_arguments + """
	--config-file           Custom config-file path. (Default: data/config.yaml)
	--min_parsimony_sites	Minimum number of parsimony informative sites in each alignments.
	--method		Trimming method. Options: trimal, aliscore (Default: read from config.yaml)
	--parameters        	Commandline arguments for trimming method. (Default: read from config.yaml)

"""

tree_help = """
phylociraptor tree - Will calculate phylogenomic trees based on supermatrix (concatenated alignment).

Usage: phylociraptor tree <arguments>

""" + standard_arguments + additional_arguments + """
	--config-file           Custom config-file path. (Default: data/config.yaml)
"""

sptree_help = """
phylociraptor speciestree - Will calculate single-gene trees and a species tree.

Usage: phylociraptor speciestree <arguments>

""" + standard_arguments + additional_arguments + """
	--config-file           Custom config-file path. (Default: data/config.yaml)
"""

njtree_help = """
phylociraptor njtree - Will calculate a NJ tree.

Usage: phylociraptor njtree <arguments>

""" + standard_arguments + additional_arguments + """
	--config-file           Custom config-file path. (Default: data/config.yaml)
"""


model_help = """
phylociraptor modeltest - Will perform substitution model tests and calculate a gene tree for each alignment.

Usage: phylociraptor modeltest <arguments>

""" + standard_arguments + additional_arguments + """
	--config-file           Custom config-file path. (Default: data/config.yaml)
"""

report_help = """
phylociraptor report - Will create a HTML report of the run

Usage: phylociraptor report <arguments>

Argumemts:
	--verbose               Display more output.
	--figure		Create a single figure report.
				This only works after modeltest has finished.
	--config-file           Relative custom config-file path. Only required for --figure (Default: data/config.yaml)
	-h, --help              Display help.

"""

check_help = """
phylociraptor check - Quickly check the status of the run.

Usage: phylociraptor check <arguments>

Argumemts:
	--verbose               Display more output.
	-h, --help              Display help.

"""

util_help = """
phylociraptor util - Utilities for a-posteriori analyses of phylociraptor results.

Usage: phylociraptor util <arguments>

Argumemts:
	get-lineage		retrieve full lineage information for all included samples.
	estimate-conflict	estimate conflict between trees based on the occurence of quartets.
	plot-tree		plot one or more trees.
	plot-conflict		plot conflicts between two trees based on quartet comparison.
	plot-similarity		plot heatmap of quartet similarity.
	-h, --help              Display help.

"""

util_lineage_help = """
phylociraptor util get-lineage - Download lineage information from NCBI for included samples.

This is used in other phylociraptor utils to enhance functionality.

Usage: phylociraptor util get-lineage <arguments>

Arguments:
	-g, --genomefile	Relative path to download_genomes_statistics.txt file
				Default: results/statistics/downloaded_genomes_statistics.txt	
	-o, --outfile		Output file name.
	--quiet			Suppress on-screen output.

"""

util_estimate_conflict_help = """
phylociraptor util estimate-conflict - Estimates conflicts between trees based on quartets.

Usage: phylociraptor util estimate-conflict <arguments>

Required Arguments:
	-i, --intrees		Relative paths to input trees in results folder, seperated by commas. (Default: all)
	-o, --outprefix		Output file name prefix.
	-n, --nquartets		Total number of quartets to be calculated. Use either this or --stopby.
	-b, --stopby		Stoping criterion. There are two options:
					-b conflicts=100 stops as soon as 100 conflicts have been found.
					--stopby tipcoverage=100 stops as soon as every tip is in at least 100 quartets.
Optional Arguments:
	-s, --seed		Random seed number for reproducibility. (Default: random)
	-l, --lineagefile	Lineagefile created with phylociraptor util get-lineage.
				Mandatory when using-a/--selecttaxa.
	-a, --selecttaxa	Sample quartets only from tips belonging to specific taxa.
				Requires output from phylociraptor util get-lineage.
					Examples: --selecttaxa order=Diperta,Hymenoptera
	-t, --threads		Number of threads. (Default: 1)
	
	--quiet			Suppress on-screen output.

"""

util_plot_tree_help = """
phylociraptor util plot-tree - Plot one or more phylogenomic trees.

Usage: phylociraptor util plot-tree <arguments>

Required Arguments:
	-i, --intrees		Relative paths to input trees in results folder, separated by commas. (Default: all)
	
Optional Arguments:
	-o, --outprefix		Output file name prefix.
	-l, --lineagefile	Lineagefile created with phylociraptor util get-lineage.
	-e, --level		Taxonomic level in lineage file which should be plotted.
	-s, --seed		Random seed number for reproducibility. (Default: random)
	-g, --outgroup		Comma seperated list of tips which should be used as Outgroup.
				Trees will be rerooted accordingly.
				(Default: none)
	--quiet			Suppress on-screen output.

"""

util_plot_conflict_help = """
phylociraptor util plot-conflict - Visualizes conflicts between two trees.

Usage: phylociraptor util plot-conflict <arguments>

Required Arguments:
	-i, --intrees		Two trees for which conflicts should be visualized, separated by a dash.
				Naming follows the first column in the *.treelist.csv file.
				Example: -i T1,T5 will plot conflicts between T1 and T5 in the
					 treelist file.
	-q, --quartetfile	*.quartets.csv file from phylociraptor util estimate-conflict.
	-r, --treelist		*.treelist.csv file from phylociraptor util estimate-conflict.
	
Optional Arguments:
	-s, --seed		Random seed number for reproducibility. (Default: random)
	-l, --lineagefile	Lineagefile created with phylociraptor util get-lineage.
				Is required when -e/--level should be used.
				(Default: none)
	-e, --level		Taxonomic level in lineage file which should be plotted.
				Has to be used together with -l/--lineagefile.
				(Default: none)
	-g, --outgroup		Comma seperated list of tips which should be used as Outgroup.
				Trees will be rerooted accordingly.
				(Default: none)
	--quiet			Suppress on-screen output.

"""
util_plot_similarity_help = """
phylociraptor util plot-similarity - Create similarity heatmap for trees and calculate tip2tip distances.

Usage: phylociraptor util plot-similarity <arguments>

Required Arguments:
	-m, --simmatrix		*.similarity_matrix.csv file from phylociraptor util estimate-conflict.
	
Optional Arguments:
	-r, --treelist		*.treelist.tsv file from phylociraptor util estimate-conflict.
				When this file is provided, a tip2tip distance analysis will be performed.
				The results will be plotted as PCA and provided as CSV file.
				(Default: none)
	-s, --seed		Random seed for reproducibility. (Default: random)
	-n, --ndistances	Number of tip2tip distances to be calculated. (Default: all)
				Only meaningful when combined with -r/--treelist
	-t, --threads		Number of threads for tip2tip distance analysis. (Default: 1)
				Only meaningful when combined with -r/--treelist
	--quiet			Suppress on-screen output.

"""
