
outfile_dict = {
	"setup": [],
	"orthology": ["results/checkpoints/modes/phylogenomics_setup.done"],
	"filter-orthology": ["results/checkpoints/modes/orthology.HASH.done"],
	"align": ["results/checkpoints/modes/filter_orthology.HASH.done"],
	"filter-align": ["results/checkpoints/modes/align.HASH.done"],
	"speciestree": ["results/checkpoints/modes/modeltest.HASH.done"],
	"njtree": ["results/checkpoints/modes/modeltest.HASH.done"],
	"mltree": ["results/checkpoints/modes/modeltest.HASH.done"],
	"modeltest": ["results/checkpoints/modes/filter_align.HASH.done"],
	}

#steps_to_check = ["setup", "orthology", "filter-orthology", "align", "filter-align", "njtree", "modeltest", "mltree", "speciestree"]
steps_to_check = ["setup", "orthology", "filter-orthology", "align", "filter-align", "njtree", "modeltest", "speciestree", "mltree"]
checkpoint_file_dict = {
	"setup": "results/checkpoints/modes/phylogenomics_setup.done",
	"orthology": "results/checkpoints/modes/orthology.HASH.done",
	"filter-orthology": "results/checkpoints/modes/filter_orthology.HASH.done",
	"align": "results/checkpoints/modes/align.done",
	"filter-align": "results/checkpoints/modes/filter_align.done",
	"speciestree": "results/checkpoints/modes/speciestree.done",
	"njtree": "results/checkpoints/modes/njtree.done",
	"mltree": "results/checkpoints/modes/trees.done",
	"modeltest": "results/checkpoints/modes/modeltest.done"
}
	

outdir_dict = {
	"setup": ["results/orthology/busco/busco_set", "results/assemblies", "results/downloaded_genomes"],
	"orthology": ["results/orthology/busco"],
	"align": ["results/alignments"],
	"filter-orthology": ["results"],
	"align": ["results"],
	"filter-align": ["results/alignments/trimmed", "results/alignments/filtered"],
	"speciestree": [""], #donefile will have to do as check, because there are several possible output folder combinations for this step
	"njtree": [""], #donefile will have to do as check, because there are several possible output folder combinations for this step
	"mltree": [""], #donefile will have to do as check, because there are several possible output folder combinations for this step
	"modeltest": ["results/modeltest"]
}
