.. role:: bash(code)
   :language: bash


.. _BUSCO: https://busco-archive.ezlab.org/
.. _YAML: https://en.wikipedia.org/wiki/YAML
.. _Augustus: http://bioinf.uni-greifswald.de/augustus/
.. _mafft: https://mafft.cbrc.jp/alignment/server/
.. _trimal: http://trimal.cgenomics.org/
.. _aliscore: https://bonn.leibniz-lib.de/en/research/research-centres-and-groups/aliscore
.. _raxml-ng: https://github.com/amkozlov/raxml-ng
.. _iqtree: http://www.iqtree.org/
.. _astral: https://github.com/smirarab/ASTRAL
.. _NCBI Genome Browser: https://www.ncbi.nlm.nih.gov/genome/browse#!/overview/
.. _biomartr: https://github.com/ropensci/biomartr
.. _quicktree: https://github.com/khowe/quicktree
 
======================
phylociraptor runmodes
======================

A phylociraptor analysis is split into different parts, which also correspond to typical steps in a phylogenomic analysis. Each part is represented by it's own runmode and can be run independently. However runmodes can depend on each other (eg. for calculating a tree you will first have to create alignments) and are therefore typically executed in this order:

.. code-block:: bash

	$ phylociraptor setup
	$ phylociraptor orthology
	$ phylociraptor filter-orthology
	$ phylociraptor align
	$ phylociraptor filter-align
	$ phylociraptor model
	$ phylociraptor njtree
        $ phylocirpator mltree
        $ phylocirpator speciestree
	$ phylociraptor check
	$ phylociraptor report

We will now talk about the different runmodes individually.


------------------------------------
orthology (Orthology inference)
------------------------------------

A prerequisite of gene-based phylogenomics is to establish orthologous relationships and identify a set of suitable genes (usually single-copy orthologs). phylociraptor currently uses the BUSCO approach to search for complete single-copy orthologs (genes which are present with only a single copy in eachh genome) in the specified set of assemblies. BUSCO is well established and widely used in contemporary genomic analyses. Phylociraptor will automatically run BUSCO on the set of specified genomes, extract single-copy genes and combine them into individual files for each BUSCO gene. This step is invoked by running :bash:`phylociraptor orthology`.  
After this step you will see a new directory :bash:`results/busco_sequences` which contains a file for each BUSCO gene. Each gene file contains all the sequences of that gene which were found in the set of specified genomes. This selection is based on a table produced by phylocirpator which is found in :bash:`results/busco_table`. 
There are many additional approaches to infer orthology and in the future we plan to add some of them to phylociraptor.

For additional information on the inner workings of BUSCO please go `here <https://busco-archive.ezlab.org/>`_.

--------------------------------------
filter-orthology Filter orthologs)
--------------------------------------

When orthology inferrence has successfully finished (after running :bash:`phylocirpator orthology`) it is necessary to filter the results. Due to a number of reasons (eg. low assembly quality, poor representation of taxonomic group in BUSCO set, etc.)  it can happen that BUSCO performance is low and we therefore want to remove samples with too-low performance from downstream analyses. To do this, phylociraptor offers the :bash:`phylociraptor filter-orthology` runmode. Here phylociraptor will perform two filtering steps, which can be set in the :bash:`config.yaml` file under the section :bash:`filtering` (see also the corresponding section above).

1. Remove samples with low overall BUSCO performance based on the percentage of found complete single-copy orthologs. (sample-based filtering)
2. Remove BUSCO genes which have only been found in a small number of samples. (gene-based filtering)

For the first filtering step you need to specify a number between 0 and 1. 0 means include all species and 1 means include only the samples for which all BUSCO genes were found in a single copy. Both scenarios are unrealistic. On the one hand, including also species for which only a very small number of single copy orthologous genes have been found could influence phylogenetic placement and quality of the tree. On the other hand, by using 1 you assume that BUSCO was able to identify each BUSCO gene in each sample, which is very unlikely as well.  
As you can see, filtering is a trade-off. Increasing this value will lower the number of samples included in the analysis, while keeping it too low could impact phylogenetic placement.  

The second filtering step is also important. For each gene used in a phylogenomic analysis you will want a reasonably high number of sequences from different samples and a small number of missing data. The number you use here should be between 3 (the minimum number of sequences you need to calculate a phylogeny) and the number of samples you have in your dataset. 

-------------------------------------
align (Create and trim alignments)
-------------------------------------

During this step phylociraptor creates individual alignments for each recovered single-copy orthologous gene. Alignment is currently done using `mafft`_ but we plan to add additional aligners in the future. According to the setting specified in the :bash:`config.yaml` file (see above) mafft will be run for each gene. Each alignment will be placed in the directory `results/alignments/full`. Individual alignments are in FASTA format and can be downloaded and inspected.

The corresponding runmode of phylociraptor is :bash:`phylociraptor align`

.. note::

   Alignment and trimming are executed together in the runmode :bash:`-m align` . 

After alignments have been generated, each alignment is trimmed to filter out positions and sequences (depending on the selected trimming strategy). Phylociraptor supports `trimal`_ and AliScore/Alicut for alignment trimming.

-----------------------------------
filter-align (Filter alignments)
-----------------------------------

When alignment is finished, phylociraptor provides an additional step to filter alignments by running :bash:`phylociraptor filter-align`. This runmode performs two steps. First it will trim alignments using `trimal`_ or `aliscore`_ or both depending on what was specified in the `config.yaml` file. Trimal and aliscore will remove sites and/or sequences from the alignments based on the specified settings. Thus, as a second step after trimming, the alignments have to be reevalueated if they should be kept for the subsequent phylogenetic reconstructions. Trimmed alignments are filtered based on two criteria:

1. First, alignments will be filtered based on the number of parsimony informative sites in the alignment. This value can be set in the :bash:`config.yaml` file.
2. Second, alignments will be filtered again for the number of sequences they contain. This step is similar to the filtering down in :bash:`phylociraptor filter-orthology`. It is necessary to do this twice, since the number of sequences in each alignment could have changed after trimming.

phylociraptor will output trimmed alignments to :bash:`results/alignments/trimmed` and filtered alignments to :bash:`results/alignments/filtered`. The files in the later folder will be used for subsequent steps.

-------------------------------------
model (Substitution model testing)
-------------------------------------

During this step phylociraptor will determine the best substitution model for each gene. It uses the `iqtree -m TESTONLY` mode from IQ-Tree. Look `here <http://www.iqtree.org/doc/Tutorial#choosing-the-right-substitution-model>`_ for additional information on how this works.

The information on the best substitution model is available in the `results/modeltest` directory. Due to the reason that iqtree and raxml support different numbers of substitution models and because they are named differently, some model names infered by iqtree may be incompatible with raxml.
phylociraptor tries to resolve these discrepancies automatically to make sure that the models inferred with iqtree also work with raxml. This does not work in every case and it is hard to anticipate which models work and which don't. If you encounter a problematic model with raxml please let us know.

-------------------------------------
mltree (Calculate ML phylogenies)
-------------------------------------

This runmode allows to calculate maximum-likelihood trees from concatenated (supermatrix) alignments of all genes which pass the filtering step.
The trees can be calculated using iqtree or raxml. phylociraptor will create the partition file necessary for raxml (iqtree does this automatically).
If `phylocripator model` has been run before, phylociraptor will pass the best models on to the tree inference software to save time.
Otherwise a model (or modeltest) can be specified in the `config.yaml` file.

-----------------------------------------
speciestree (Calculate species trees)
-----------------------------------------

phylociraptor calculates a species tree using `astral`. Astral takes pre-calculated gene trees as input. Phylociraptor checks if gene-trees have been already calculated and creates them in case they are not yet available.

To calculate individal gene trees phylociraptor uses iqtree.

------------------------------------------
njtree (Calculate NJ tree)
------------------------------------------

To get a fast first tree you can run `phylociraptor njtree`. This will calculate a Neighbor-Joining tree using `quicktree`. This usually only takes seconds.

------------------------------------------
check (Check status of the run)
------------------------------------------

`phylciraptor check` will give a quick (and dirty) overview about which steps have already been run. This can be helpful to keep an overview of how many steps have already finished in cases where there are hundreds or thousands of jobs submitted to a cluster. `phylociraptor check` is however superficial and can only help to quickly assess the status of the pipeline. It shows DONE for each step that has finsihed, INCOMPLETE for steps which have either not finished or not run at all and NOT EVALUATED for steps which require other steps to have been run before.


------------------------------------------
report (Create an HTML report summarizing the results)
------------------------------------------

`phylociraptor report` will create an HTML report of the run. It includes statistics calculated during each step. It can be run after each step of phylociraptor and is intended to help to decide on meaningful setting for the next analysis steps.



