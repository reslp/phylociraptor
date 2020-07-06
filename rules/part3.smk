if config["phylogeny"]["concat"] == "yes":
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
                       bb = "1000",
                       m = "WAG"
               threads:
                       config["iqtree"]["threads"]
               shell:
                       """
                       rm -rf results/phylogeny/concatenated/algn
                       mkdir -p results/phylogeny/concatenated/
                       cd results/phylogeny/concatenated/
                       mkdir algn
                       cp {params.wd}/results/trimmed_alignments/*.fas algn
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


