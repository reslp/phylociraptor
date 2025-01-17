
.. role:: bash(code)
   :language: bash


.. _BUSCO: https://busco-archive.ezlab.org/
.. _YAML: https://en.wikipedia.org/wiki/YAML
.. _Augustus: http://bioinf.uni-greifswald.de/augustus/
.. _mafft: https://mafft.cbrc.jp/alignment/server/
.. _trimal: http://trimal.cgenomics.org/
.. _raxml-ng: https://github.com/amkozlov/raxml-ng
.. _iqtree: http://www.iqtree.org/
.. _astral: https://github.com/smirarab/ASTRAL
.. _NCBI Genome Browser: https://www.ncbi.nlm.nih.gov/genome/browse#!/overview/
.. _biomartr: https://github.com/ropensci/biomartr
 
====================================
Understanding phylociraptor's output
====================================

Phylociraptor produces a lot of different output. The philsophy behind this is to make each step transpartent so that every file can be traced back to the original input files. Phylociraptor uses a pre-defined directory structure which is identical for every run. All results are contained in the :bash:`results` directory in the phylociraptor folder.

This is the content of the results folder when every step of phylociraptor has finished successfully:

.. code-block:: bash

    $ ls results/
    alignments  assemblies  checkpoints  downloaded_genomes  modeltest  orthology  phylogeny-0  phylogeny-60  phylogeny-70  report.html  report_figure.pdy  statistics

We will go through every folder now in alphabethical order.


.. caution::

   If you compare what is written here to your own phylociraptor results you will notice many directory- and filenames additionally contain combinations of numbers and letters.
   These so-called `hashes <https://en.wikipedia.org/wiki/Hash_function>`_ are used to track files and folder of analyses created with different parameter combinations from the config file. This is necessary to make analyses reproducible down to every set parameter in your analysis. For simplicity we have omitted hashes here.

 

----------
alignments
----------

This folder will only show up in results after runnning :bash:`phylociraptor align`.

.. code-block:: bash

    $ ls alignments/
    filtered  full  trimmed 

At the end of a phylociraptor run the  folder alignments should include three folders. Each of them include another set of folders depending on which aligners have been used.Let's have a look into the full folder:

.. code-block:: bash

    $ ls alignments/full
    clustalo  mafft  muscle 

The folder :bash:`full` includes a folder for each aligner. These folders contain the actual full (untrimmed) alignments for each orthologous gene in FASTA format.

The same logic applies to the :bash:`trimmed` and :bash:`filtered` folders, however the subfolders there are named slightly different to also incorporate the used trimmer:

.. code-block:: bash

   $ ls alignments/trimmed
   clustalo-aliscore  clustalo-trimal  mafft-aliscore  mafft-trimal  muscle-aliscore  muscle-trimal

As you can see we have six folders now. This is because phylocirpator will perform trimming for every aligner and trimmer combination and in this example we used three aligners (clustalo, mafft and muscle) and two trimmers (trimal, aliscore).

----------
assemblies
----------

This folder is created during the :bash:`phylociraptor setup` step.

It contains symlinks to all input data used for a phylociraptor run. Regardless of its name, it may contain (a mixture of) genomes, transcriptomes and proteoms.

.. code-block:: bash

   $ ls -l assemblies
   total 2
   lrwxrwxrwx 1 reslp p71312 101 30. Mär 10:14 Amanita_muscaria.fasta.gz -> /binfl/lv71312/reslp/phylociraptor/results/downloaded_genomes/Amanita_muscaria_genomic_genbank.fna.gz
   lrwxrwxrwx 1 reslp p71312 102 30. Mär 10:14 Neurospora_crassa.fasta.gz -> /binfl/lv71312/reslp/phylociraptor/results/downloaded_genomes/Neurospora_crassa_genomic_genbank.fna.gz
   lrwxrwxrwx 1 reslp p71312 103 30. Mär 10:14 Tuber_melanosporum.fasta.gz -> /binfl/lv71312/reslp/phylociraptor/results/downloaded_genomes/Tuber_melanosporum_genomic_genbank.fna.gz
   lrwxrwxrwx 1 reslp p71312 101 30. Mär 10:14 Usnea_hakonensis.fasta.gz -> /binfl/lv71312/reslp/phylociraptor/results/downloaded_genomes/Usnea_hakonensis_genomic_genbank.fna.gz

As you can see from the listing above, in this case the directory contains symlinks to four genomes. You can also see the actual location of the files is in :bash:`results/downloaded_genomes` we will get to this folder later. Phylociraptor uses symlinks here to save disk-space will maintaining a consitent location and naming scheme for all files in the directory.


-----------
checkpoints
-----------

This folder contains a set of files which help phylociraptor to keep track about which parts have been finished. Usually these files are irrelevant for the user. They are only for internal use of the pipeline.

------------------
downloaded_genomes
------------------

This folder will only appear when :bash:`phylociraptor setup` has finished. As mentioned above, it contains downloaded genomes from the NCBI database and additional information for each genome which will be displayed in the run report.

.. code-block:: bash

   $ ls results/downloaded_genomes
   Amanita_muscaria_db_genbank.tsv          download_overview.txt                     not_downloaded.txt                 Tuber_melanosporum_genomic_genbank.fna.gz
   Amanita_muscaria_genomic_genbank.fna.gz  Neurospora_crassa_db_genbank.tsv          successfully_downloaded_1.txt      Usnea_hakonensis_db_genbank.tsv
   assembly_summary_genbank.txt             Neurospora_crassa_genomic_genbank.fna.gz  successfully_downloaded.txt        Usnea_hakonensis_genomic_genbank.fna.gz
   download_overview_1.txt                  not_downloaded_1.txt                      Tuber_melanosporum_db_genbank.tsv

As you can see there are several files in this folder. The fna.gz file is the gzipped genome assembly downloaded from NCBI. It is accompanied by a .tsv file which contains information about the assembly such was the accession number, the download URL, the Institution which uploaded the genome and more. The information in all .tsv files will be summarized in the HTML report. 

Additionally, the folder contains files which are for internal use in phylociraptor. They may only be relevant when debugging problems and are not directly relevant for the user. All these files end with .txt. 

---------
modeltest
---------

This folder contains all results from modeltesting and gene-tree calculations. It will only appear after running :bash:`phylociraptor modeltest`.

.. code-block:: bash

   $ ls results/modeltest
   best_models_clustalo_aliscore.txt  best_models_muscle_aliscore.txt  genetree_filter_clustalo_aliscore.txt  genetree_filter_muscle_aliscore.txt  muscle-aliscore
   best_models_clustalo_trimal.txt    best_models_muscle_trimal.txt    genetree_filter_clustalo_trimal.txt    genetree_filter_muscle_trimal.txt    muscle-trimal
   best_models_mafft_aliscore.txt     clustalo-aliscore                genetree_filter_mafft_aliscore.txt     mafft-aliscore
   best_models_mafft_trimal.txt       clustalo-trimal                  genetree_filter_mafft_trimal.txt       mafft-trimal

This folder contains several text files which summarize the best substitution models estimated by the iqtree MFP algorithm according to BIC in .txt format. These files are used for subsequent steps when creating a partitioned concatenated alignment for Maximum-Likelihood tree calculation.

Additionally the subfolders contain all gene-tree calculation results. They are consistantly named with a combination of the name of the aligner and trimmer. For example the folder clustalo-trimal contains all gene tree results for each orthologous gene aligned with clustalo and trimmed with trimal:

.. code-block:: bash

   $ ls results/modeltest/clustalo-trimal
   EOG092C5OAL  EOG092C5OPO  EOG092C5Q82  EOG092C5S4U  EOG092C5U93  EOG092C5UXM  EOG092C5V3D  EOG092C608
   $ ls results/modeltest/clustalo-trimal/EOG092C5OAL
   EOG092C5OAL.bionj  EOG092C5OAL.ckp.gz  EOG092C5OAL.contree  EOG092C5OAL.iqtree  EOG092C5OAL.log  EOG092C5OAL.mldist  EOG092C5OAL.model.gz  EOG092C5OAL.splits.nex  EOG092C5OAL.treefil

The .treefile from each individual analysis will be used later when a speciestree is calculated and when the bootstrap cut-off is applied.

---------
orthology
---------

In this folder all results from the :bash:`phylociraptor orthology` step is stored. 

.. code-block:: bash

   $ ls results/orthology
   busco
   $ ls results/orthology/busco
   busco_runs  busco_sequences  busco_sequences_deduplicated  busco_set  busco_table.txt


On the first layer there is only a single folder :bash:`busco`, which is currently the only method to infer orthologous genes for the phylogeny.
The :bash:`busco` folder contains several additional folders:

:bash:`busco_runs` contains the results of the individual busco runs for each sample in separate folders.

:bash:`busco_sequences` contains extracted BUSCO sequences from each sample combined into a single file for each BUSCO.

:bash:`busco_sequences_deduplicated` BUSCO sometimes reports more than one sequences even for single-copy BUSCOs. This folder contains FASTA files of each BUSCO gene from :bash:`busco_sequences` but with duplicates removed.

:bash:`busco_set` contains the used busco set as spcified in the config.yaml file.

Important for the pipline is the file :bash:`busco_table.txt`. It is a condensed summary of all individual BUSCO runs and is used to determine which samples will be used in subsequent steps based on BUSCO completeness.

-----------
phylogeny-* 
-----------

This is probably the most important output folder. It contains the phylogenomic trees. When phylociraptor is finish calculating trees you should see one or more folders starting with :bash:`phylogeny` in your results folder. These folder contain results from phylogenomic tree calculations. The reason there are more than one is because you can subset which genes should be used based on the mean bootstrap support value of the gene trees. So, the folder :bash:`phylogeny-80` contains phylogenetic results only based on genes for which the mean bootstrap support was >80%.

Let's have a look at a phylogeny folder more closely:

.. code-block:: bash

   $ ls results/phylogeny-0
   astral  concatenate  iqtree  njtree  raxml

As you can see there are three folders. Let us investigate the concatenate folder first, which contains information on concatenated alignments:

.. code-block:: bash

   $ ls results/phylogeny-0/concatenate
   clustalo-aliscore  clustalo-trimal  mafft-aliscore  mafft-trimal  muscle-aliscore  muscle-trimal
   $ ls results/phylogeny-0/concatenate/clustalo-aliscore
   algn  concat.fas  concat.phy  concat.sto  partitions.txt  partitions_unformated.txt  statistics.txt

The concatenate folder contains several subfolders for each aligner and trimmer combination. Each of these folders contains the following items:

The folder :bash:`align`, which contains individual trimmed and filtered single gene alignments which are used to create the concatenated alignment.

The files :bash:`concat.fas`, :bash:`concat.phy` and :bash:`concat.sto` contain the resulting concatenated alignment in FASTA, PHYLIP and STOCKHOLM format.

The files :bash:`partitions.txt` and :bash:`partitions_unformated.txt` contain the partitioning scheme in a format raxml understands. These two files could differ slightly. This is because the names used to specify substitution models differs slightly between raxml and iqtree (which we use to estimate the model). Phylociraptor converts the modelnames from iqtree notation (:bash:`partitions_unformated.txt`) to raxml notation (:bash:`partitions.txt`).

The :bash:`statistics.txt` file includes statistics of the genes used in the concatenated alignment, such as parsimony informative sites, length, etc.



The next folder we can have a look at is the :bash:`iqtree` folder. It includes results from iqtree maximum-likelihood calculations based on the concatenated datasets:

.. code-block:: bash

   $ ls results/phylogeny-0/iqtree
   clustalo-aliscore  clustalo-trimal  mafft-aliscore  mafft-trimal  muscle-aliscore  muscle-trimal
   $ ls results/phylogeny-0/iqtree/mafft-aliscore
   algn  concat.best_model.nex  concat.bionj  concat.ckp.gz  concat.contree  concat.iqtree  concat.log  concat.mldist  concat.nex  concat.splits.nex  concat.treefil

Again this folder contains subfolder for each aligner and trimmer combination you specified in the config file. Each subfolder contains all output from iqtree as well as a folder (:bash:`algn`) with all alignments used to calculate the tree.

The folder :bash:`raxml` contains all output from raxml for the different aligner and trimmer combinations.

.. code-block:: bash

   $ ls results/phylogeny-0/raxml
   clustalo-aliscore  clustalo-trimal  mafft-aliscore  mafft-trimal  muscle-aliscore  muscle-trimal
   $ ls results/phylogeny-0/raxml/mafft-trimal
   concat.fas      raxmlng.raxml.bestModel  raxmlng.raxml.bootstraps  raxmlng.raxml.mlTrees  raxmlng.raxml.startTree
   partitions.txt  raxmlng.raxml.bestTree   raxmlng.raxml.log         raxmlng.raxml.rba      raxmlng.raxml.support

It includes the :bash:`concat.fas` file which is the concatenated alignment along with the :bash:`partitions.txt` file with partition specifications as well as the output from raxml.

In the folder :bash:`astral` you will find species tree results from Astral:

.. code-block:: bash

   $ ls results/phylogeny-0/astral
   clustalo-aliscore  clustalo-trimal  mafft-aliscore  mafft-trimal  muscle-aliscore  muscle-trimal
   $ ls results/phylogeny-0/astral/muscle-trimal
   species_tree.tre  trees_muscle_trimal.tre

Again there is a folder for each aligner and trimmer combination. Each folder contains two .tre files in Newick format. The :bash:`species_tree.tre` file is the species tree calculated by astral and the file :bash:`trees_muscle_trimal.tre` contains all individual gene-trees used to calculate the species tree.

The last phylogeny results folder is the :bash:`njtree` folder which contains Neighbor-Joining trees calculated with quicktree.

.. code-block:: bash

   $ ls results/phylogeny-0/njtree
   clustalo-aliscore  clustalo-trimal  mafft-aliscore  mafft-trimal  muscle-aliscore  muscle-trimal
   $ ls results/phylogeny-0/njtree/muscle-aliscore
   njtree.tre

For each aligner and trimmer combinations it contains a single Newick tree file with the calculated NJ-tree.





 

