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
.. _clustalo: http://www.clustal.org/omega/
 
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
	$ phylociraptor modeltest
	$ phylociraptor njtree
        $ phylocirpator mltree
        $ phylocirpator bitree
        $ phylocirpator speciestree
	$ phylociraptor check
	$ phylociraptor report

We will now talk about the different runmodes individually.


---------
orthology
---------

.. note::

   This mode will infer orthology relationships and identify single-copy orthologs

A prerequisite of gene-based phylogenomics is to establish orthologous relationships and identify a set of suitable genes (usually single-copy orthologs). phylociraptor currently uses the BUSCO approach to search for complete single-copy orthologs (genes which are present with only a single copy in each genome) in the specified set of assemblies. BUSCO is well established and widely used in contemporary genomic analyses. Phylociraptor will automatically run BUSCO on the set of specified genomes, extract single-copy genes and combine them into individual files. This step is invoked by running :bash:`phylociraptor orthology`.  
After this step you will see a new directory :bash:`results/busco_sequences` which contains a file for each BUSCO gene. Each gene file contains all the sequences of that gene which were found in the set of specified genomes. This selection is based on a table produced by phylocirpator which is found in :bash:`results/orthology/busco/busco_table.txt`. 
There are many additional approaches to infer orthology and in the future we plan to add some of them to phylociraptor as alternatives to BUSCO.

For additional information on the inner workings of BUSCO please go `here <https://busco-archive.ezlab.org/>`_.

.. hint::

   Phylociraptor supports genome and transcriptome assmeblies as well as sets of protein sequences  as input.

----------------
filter-orthology
----------------

.. note::

   This mode will filter the identfied orthologs

When orthology inferrence has successfully finished (after running :bash:`phylocirpator orthology`) it is necessary to filter the results. Due to a number of reasons (eg. low assembly quality, poor representation of taxonomic group in BUSCO set, etc.) it can happen that BUSCO performance is low and we therefore want to remove samples with low BUSCO scores from downstream analyses. To do this, phylociraptor offers :bash:`phylociraptor filter-orthology`. When this is run, phylociraptor will perform two filtering steps, which can be set in the :bash:`config.yaml` file under the section :bash:`filtering` (see also the corresponding section above):

1. Sample-based filtering: Remove samples with low overall BUSCO performance based on the percentage of found complete single-copy orthologs. (cutoff)
2. Gene-based filtering: Remove BUSCO genes which have only been found in a specified number of samples. (minsp) 

For the first filtering step you need to specify a number between 0 and 1. 0 includes all species and 1 includes only the samples for which ALL BUSCO genes were found in a single copy (this species has 100% BUSCO completeness). Both scenarios are unrealistic. On the one hand, including also species for which only a very small number of single copy orthologous genes have been found could influence phylogenetic placement and the quality of the tree. On the other hand, by using 1 you assume that BUSCO was able to identify each BUSCO gene in each sample, which is very unlikely as well.  
As you can see, filtering is a trade-off. Increasing this value will lower the number of samples included in the analysis, while keeping it too low could impact phylogenetic placement.  

The second filtering step is also important. For each gene used in a phylogenomic analysis you will want a reasonably high number of sequences from different samples and a small number of missing data. The number you use here should be between 3 (the minimum number of sequences you need to calculate a phylogeny) and the number of samples you have in your dataset. Mind you that setting this to a very high value will potentially exclude many genes from subsequent analysis.

Results from this step can be found here: :bash:`results/orthology/busco/busco_sequences_deduplicated`.

-----
align
-----


.. note::

   This mode will create multiple sequence alignments of all orthologs

During this step phylociraptor creates individual alignments for each recovered single-copy orthologous gene. Alignments can be created using `mafft`_  or `clustalo`_ . According to the setting specified in the :bash:`config.yaml` file (see above) the aligner will be run for each gene. Each alignment will be placed in the directory `results/alignments/full`. Individual alignments are in FASTA format and can be downloaded and inspected.

The corresponding runmode of phylociraptor is :bash:`phylociraptor align`

------------
filter-align
------------

.. note::

   This mode will trim and filter multiple sequence

When alignment is finished, phylociraptor provides an additional step to filter alignments by running :bash:`phylociraptor filter-align`. This runmode performs two steps. First it will trim alignments using `trimal`_ or `aliscore`_ or both depending on what was specified in the `config.yaml` file. Trimal and aliscore will remove sites and/or sequences from the alignments based on the specified settings. Thus, as a second step after trimming, the alignments have to be reevalueated if they should be kept for the subsequent phylogenetic reconstructions. Trimmed alignments are filtered based on two criteria:

1. First, alignments will be filtered based on the number of parsimony informative sites in the alignment. This value can be set in the :bash:`config.yaml` file.
2. Second, alignments will be filtered again for the number of sequences they contain. This step is similar to the filtering down in :bash:`phylociraptor filter-orthology`. It is necessary to do this twice, since the number of sequences in each alignment could have changed after trimming.

phylociraptor will output trimmed alignments to :bash:`results/alignments/trimmed` and filtered alignments to :bash:`results/alignments/filtered`. The files in the later folder will be used for subsequent steps.

---------
modeltest
---------

.. note::

   This mode will perform substitution model testing

Phylociraptor can determine the best substitution model for each gene. It uses the :bash:`iqtree` to infer the best substitution model and it will use this model to calculate a maximum-likelihood gene tree. Look `here <http://www.iqtree.org/doc/Tutorial#choosing-the-right-substitution-model>`_ for additional information on how this works. 

The information on the best substitution model is available in the `results/modeltest` directory. Due to the reason that iqtree and raxml support different numbers of substitution models and because they are named differently, some model names infered by iqtree may be incompatible with raxml.
phylociraptor tries to resolve these discrepancies automatically to make sure that the models inferred with iqtree also work with raxml. This does not work in every case and it is hard to anticipate which models work and which don't. If you encounter a problematic model with raxml please let us know by raising an issue on GitHub.

------
mltree
------

.. note::

   Calculate maximum-likelihood trees


This runmode allows to calculate maximum-likelihood trees from concatenated (supermatrix) alignments of all genes which pass the filtering step.
The trees can be calculated using iqtree or raxml. phylociraptor will create the partition file necessary for raxml (iqtree does this automatically) and create a concatenated alignment of all single-gene alignments which survived the filtering step. Results from this analysis step can be found in :bash:`results/phylogeny/concatenate`. 

If `phylocripator model` has been run before, phylociraptor will pass the best models estimated in this step on to the tree inference software to save time.
Otherwise a model (or modeltest) can be specified in the `config.yaml` file.

-----------
speciestree
-----------

.. note::

   Calculate species trees

phylociraptor calculates species trees using `astral`. Astral takes pre-calculated gene trees as input. Phylociraptor checks if gene-trees have been already calculated and creates them in case they are not yet available.

Individal gene trees are calculated with iqtree.

------
bitree
------

.. note::

   Calculate bayesian trees


phylociraptor calculates bayesian trees using phylobayes. Phylobayes uses concatenated (supermatrix) alignments all genes which pass the filtering steps as input.
Phylociraptor can start multiple chains as specified in the config file and it allows to monitor running chains and calculating consensus trees using `phylocripator util`.

------
njtree
------

.. note::

   Calculate neighbor-joining trees

To get a fast first tree you can run `phylociraptor njtree`. This will calculate a Neighbor-Joining tree using `quicktree`. This usually only takes seconds and even on a cluster it is typically not necessary to use batch job submission.

-----
check
-----

.. note::

   Get an overview of your phylociraptor analysis

`phylciraptor check` will give a quick (and dirty) overview about which steps have already been run. This can be helpful to keep an overview of how many steps have already finished in cases where there are hundreds or thousands of jobs submitted to a cluster. `phylociraptor check` is however superficial and can only help to quickly assess the status of the pipeline. It shows DONE for each step that has finsihed, INCOMPLETE for steps which have either not finished or not run at all and NOT EVALUATED for steps which require other steps to have been run before.

------
report
------

.. note::

   Create a report for your phylociraptor analysis

`phylociraptor report` will create an HTML report of the run. It includes statistics calculated during each step. It can be run after each step of phylociraptor and is intended to help to decide on meaningful setting for the next analysis steps.



