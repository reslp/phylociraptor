
outfile_dict = {
	"setup": ["results/checkpoints/modes/phylogenomics_setup.done"],
	"orthology": ["results/checkpoints/modes/phylogenomics_setup.done"],
	"filter-orthology": ["results/checkpoints/modes/phylogenomics_setup.done", "results/checkpoints/modes/orthology.done"],
	"align": ["results/checkpoints/modes/phylogenomics_setup.done", "results/checkpoints/modes/orthology.done", "results/checkpoints/modes/filter_orthology.done"],
	"filter-align": ["results/checkpoints/modes/phylogenomics_setup.done", "results/checkpoints/modes/orthology.done", "results/checkpoints/modes/filter_orthology.done", "results/checkpoints/modes/align.done"],
	"speciestree": ["results/checkpoints/modes/phylogenomics_setup.done", "results/checkpoints/modes/orthology.done", "results/checkpoints/modes/filter_orthology.done", "results/checkpoints/modes/align.done", "results/checkpoints/modes/filter_align.done"],
	"njtree": ["results/checkpoints/modes/phylogenomics_setup.done", "results/checkpoints/modes/orthology.done", "results/checkpoints/modes/filter_orthology.done", "results/checkpoints/modes/align.done", "results/checkpoints/modes/filter_align.done"],
	"mltree": ["results/checkpoints/modes/phylogenomics_setup.done", "results/checkpoints/modes/orthology.done", "results/checkpoints/modes/filter_orthology.done", "results/checkpoints/modes/align.done", "results/checkpoints/modes/filter_align.done"],
	"modeltest": ["results/checkpoints/modes/phylogenomics_setup.done", "results/checkpoints/modes/orthology.done", "results/checkpoints/modes/filter_orthology.done", "results/checkpoints/modes/align.done", "results/checkpoints/modes/filter_align.done"],
	"report": ["results/checkpoints/modes/phylogenomics_setup.done"]
	}

steps_to_check = ["setup", "orthology", "filter-orthology", "align", "filter-align", "njtree", "modeltest", "mltree", "speciestree"]
checkpoint_file_dict = {
	"setup": "results/checkpoints/modes/phylogenomics_setup.done",
	"orthology": "results/checkpoints/modes/orthology.done",
	"filter-orthology": "results/checkpoints/modes/filter_orthology.done",
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
