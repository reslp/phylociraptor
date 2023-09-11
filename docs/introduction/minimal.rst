.. role:: bash(code)
    :language: bash

============================
Tutorial - A minimal example
============================

This document guides you through a first example on how to run phylociraptor. We'll analyse a minimal dataset consists of 6 fungal genomes, three belonging to the Ascomycota and three to
Basidiomycota which will be downloaded during the first step of the run. Downstream analyses are carried out on 10 single-copy orthologous genes which will be identified in each genome.
Data for all genes will be aligned and trimmed and we will calculate single-gene trees for each gene. After this, we will perform Maximum-Likelihood analyses on the full dataset and infer species trees.

Running through the main parts of the tutorial should take you about 1 hour to complete. Then there are further things to play around with if you don't want to stop at this point.

------------------------
Dataset and preparations
------------------------

.. note::
        
	The tutorial below assumes you have three things set up:
	 - Snakemake (tested with version 6.0.2; perhaps loaded as a conda environment)
	 - Singularity (tested with version 3.11.4-focal)
	 - You have a local copy of the phylociraptor repo on your system and your present working directory (command prompt) is at the base level of this repo



Input files for this dataset are located in :bash:`data/test_cases/minimal`.
Two config files are mandatory:
 - The *settings file* contains settings applied in the following analyses.
 - The *data file* specifies the data to be processed. 

Check out the content of these files, e.g. with the following commands.

.. code-block:: bash
        
        $ cat data/test_cases/minimal/config.yaml # settings file
        $ cat data/test_cases/minimal/minimal.csv # data file


-----------------------
Analysis setup
-----------------------

At first it is necessary to download the genomes which should be analyzed. Specific accession numbers to be downloaded from NCBI are specified in the *data file*.

You can start the process with the following command (runtime: about 3 mintues, depending on the speed of your connection):

.. code-block:: bash

        $ ./phylociraptor setup --config-file data/test_cases/minimal/config.yaml -t local=4 --verbose

This command runs in local mode (without job submission to an HPC job scheduling system) and uses up to 4 CPU cores (``-t local=4``). Since in the *settings file* we have specified :bash:`concurrency: 4`, the data download will be performed in four parallel batches. So, with four cores available, four downloads should run in parallel. If you had limited the use of resources to, say ``-t local=2``, then the download would start with two processes in parallel, and then phylociraptor starts the next two whenever cores become available. 

.. note::
        
        The flag :bash:`--verbose` in the above command and all the commands below can also be omitted. It is used here to viualize pipeline output on screen.

The setup step above has downloaded a bunch of genomes from NCBI, according to accession numbers we've specified in the *data file*. The downloaded genome assemblies are in :bash:`results/assemblies`. Metadata of downloaded genomes is in :bash:`results/downloaded_genomes`.
Also, it downloaded a BUSCO reference gene set. We had specified ``fungi_odb9`` in the *settings file*. This set contains 290 genes. For this tutorial we want to limit our analyses to 10 random genes from this set. Phylociraptor has a utility to subsample an existing BUSCO set. The following command will take the original BUSCO set, pick 10 random genes and create a new BUSCO set. To ensure that we'll all get the same set of genes we specify a seed for the random number generator. 

Let's run (runtime: 2 seconds):

.. code-block:: bash

        $ ./phylociraptor util modify-busco -b fungi_odb9 -n 10 --seed 42

You will be informed that the new BUSCO set is called ``fungi_odb9-seed-42-genes-10``. We'll have to adjust our *settings file* accordingly to use this set for subsequent analyses. Let's make a copy of the file.

.. code-block:: bash

        $ cp data/test_cases/minimal/config.yaml data/config.yaml 

Then, make the change in the file either by editing it with your favourite text editor or do it programmatically - we love a good `sed`..

.. code-block:: bash

        $ sed -i 's/fungi_odb9/fungi_odb9-seed-42-genes-10/' data/config.yaml

.. note::
        
        If you wanted a specific set of genes from the BUSCO set you could also specify this explicitly.

For future reference:

.. code-block:: bash

	$ ./phylociraptor util modify-busco -b fungi_odb9 -g EOG092C5OAL,EOG092C5OPO,EOG092C5Q82,EOG092C5S4U,EOG092C5T6H


--------------------------------------------------
Identify single-copy orthologs
--------------------------------------------------

Next, we'll infer single-copy orthologs in the six selected genomes via BUSCO. 

Run this command (runtime: 6 minutes):

.. code-block:: bash
        
        $ ./phylociraptor orthology --verbose -t local=4

Notice, that we don't specify the *settings file* explicitly in the above command (unlike when we called ``phylociraptor setup``). Now, phylociraptor uses the default, which is `data/config.yaml`. Remember that we put a copy of the *settings file* there.

.. note::
        
        We will run this command using four threads as indicated with ``-t local=4``. It is also possible to omit the the number of threads and use just ``-t local``. In this case phylociraptor will use as many threads as are available on your machine/server/node. If you inspect the *settings file* you'll see that we had specified two threads for BUSCO. With the resources specified in the above command (``-t local=4``) phylociraptor will thus run two BUSCO jobs, each using two threads in parallel, whenever possible. Parallelization could be achieved at two levels, 1) the software and 2) the number of instances of a given step that are being run in parallel.


Results of this step can be found in :bash:`results/orthology`.


--------------------------------
Filter orthology results
--------------------------------

This step is used to filter orthology results based on how many orthologs where found in each genome. The process will again be split into four batches. Consult the *settings file* to see the filters we have chosen. 

To run this command (runtime: 30 seconds):

.. code-block:: bash
        
        $ ./phylociraptor filter-orthology --verbose -t local=4


.. note::
       
        In this example nothing will be filtered since the ten selected genes should be present in all genomes.


------------------------------
Align single copy orthologs
------------------------------

In this step we will create alignments of the single-copy orthologs recovered in each genome using *MAFFT* and *Clustal omega*. Note, that, since we have specified a single thread for alignment jobs in the *settings file*, this (`-t local=4`) will run up to 4 alignments in parallel at a given moment.

Let's go (runtime: 1 minute): 


.. code-block:: bash

        $ ./phylociraptor align --verbose -t local=4

Results will be in :bash:`results/alignments/full`. 


------------------------------
Filter alignments
------------------------------

Now, each alignment will be trimmed using *Aliscore/Alicut* and *trimAl* to remove poorly aligned regions. In the same step we'll filter alignments based on their information content (number of parsimony informative sites) and compositional heterogeneity. This is done for each aligner and trimmer combination, so after this step you will have four sets of alignments. 

Let's trim/filter our alignments (runtime: 2 minutes): 

 
.. code-block:: bash

        $ ./phylociraptor filter-align --verbose -t local=4

Results of this step will be in:
 - :bash:`results/alignments/trimmed` (after trimming), and
 - :bash:`results/alignments/filtered` (alignments passing filters). 

The previous step will have filtered out a few alignments. Let's see what happened, for example, with the set processed through *Clustal omega* and *trimAl*:

.. code-block:: bash

        $ cat results/statistics/filter-clustalo-trimal.d0171a17c0/alignment_filter_information_trimal_clustalo.txt



----------------------------------
Infer the best substitution models
----------------------------------

Using filtered alignments we will now infer the best substitution model for each alignment and also calculate single gene-trees with IQ-Tree. This is done for each alignment in every aligner and trimmer combination.

Let's go (runtime: 4 minutes):


.. code-block:: bash

        $ ./phylociraptor modeltest --verbose -t local=4

The results will be located in :bash:`results/modeltest`

----------------------------------------
Calculate a full Maximum-Likelihood tree
----------------------------------------

Now it is time to calculate Maximum-Likelihood trees based on concatenated supermatrices of the individual gene alignments. *IQ-Tree* will infer trees for every aligner and trimmer combination. The analyses will be partitioned using the best substitution models inferred during the step above. Additionally, phylociraptor will only include genes for which the gene tree showed average bootstrap support above specified bootstrap values. In this example the average bootstrap filters used are 50, 60 and 70, as specified in the *settings file*.

The command is as follows (runtime: 4 minutes):

.. code-block:: bash

        $ ./phylociraptor mltree --verbose -t local=4


Results are in:
 - :bash:`results/phylogeny/iqtree/bootstrap-cutoff-50/`
 - :bash:`results/phylogeny/iqtree/bootstrap-cutoff-60/`
 - :bash:`results/phylogeny/iqtree/bootstrap-cutoff-70/`.


----------------------------------------
Calculate a species tree
----------------------------------------

Finally, we will calculate species trees for every aligner and trimmer combination and every bootstrap cutoff value.

This should be fast (runtime: 1 minute).

.. code-block:: bash

        $ ./phylociraptor speciestree --verbose -t local=4


Results are in:
 - :bash:`results/phylogeny/astral/bootstrap-cutoff-50/`
 - :bash:`results/phylogeny/astral/bootstrap-cutoff-60/`
 - :bash:`results/phylogeny/astral/bootstrap-cutoff-70/`.


----------------------------------------
Generate a report
----------------------------------------

.. note::

	The command below could be run at any point during the analyses to see the current status/stats.

Let's create a HTML report with the results of your analyses (runtime: 1 minute). 

.. code-block:: bash

        $ ./phylociraptor report

You can find the report here: :bash:`results/report.html`. Do have a look.

Let's recap what the nine phylociraptor commands that we've executed above have generated for us:

- 6 fungal genome assemblies downloaded
- 10 single-copy orthologs identified in each genome
- 20 multiple sequence alignments using mafft (n=10) and clustal omega (n=10)
- 24 trimmed alignments using trimal and aliscore/alicut (could have been 40, but 16 were eliminated by our filters)
- 24 maximum-likelihood gene trees calculated with iqtree
- 12 concatenated maximum-likelihood phylogenies using 3 different bootstrap-cutoff values
- 12 species trees using 3 different bootstrap-cutoff values
- 1 comprehensive report of our analyses in HTML format


**Kudos!**

.. note::
        
        Thanks for going through this tutorial with us. If you don't have enough yet, there's a few more things we'd like to show you below.


----------------------------------------
Further exploration of software- and parameter space
----------------------------------------

In the above tutorial we had not yet enabled all software implemented in phylociraptor. Let's also trim our alignments with the third piece of software implemented.

Enable trimming with *BMGE*, by adjusting the *settings file*. Open ``data/config.yaml`` in your favourite text editor search for the section ``trimming:``, and change:

.. code-block:: bash

        trimming:
	    method: ["trimal", "aliscore"]

to:

.. code-block:: bash

        trimming:
	    method: ["trimal", "aliscore", "bmge"]

Now, let's see what would happen if you reran the alignment trimming step after this change. Note, that we add ``--dry`` to the command, which will result in a so-called *dry-run*, i.e. don't actually execute, but just show me what would happen.

.. code-block:: bash

        $ ./phylociraptor filter-align --verbose -t local=4 --dry

You'll see that phylociraptor proposes to run the necessary steps to add *BMGE* to our analyses.

Let's do it for real (runtime: 1 minute):

.. code-block:: bash

        $ ./phylociraptor filter-align --verbose -t local=4

Continue with modeltesting for the new alignments (runtime: 2 minutes), 

.. code-block:: bash

        $ ./phylociraptor modeltest --verbose -t local=4

and the inference of maximum-likelihood trees using *IQ-Tree* for the new datasets (runtime: 3 minutes).

.. code-block:: bash

        $ ./phylociraptor mltree --verbose -t local=4

Neat, no? Now, let's say we want to also include the aligner *MUSCLE*.

Open ``data/config.yaml`` in your favourite text editor search for the section ``alignment:``, and change:

.. code-block:: bash

        alignment:
            method: ["mafft", "clustalo"]

to:

.. code-block:: bash

        alignment:
            method: ["mafft", "clustalo", "muscle"]

Rerun the alignment (add *MUSCLE*), filter-align (trim/filter for all new combinations), modeltest (for all new combinations) and mltree (for all new combinations) steps (runtime: 5 minutes).

.. code-block:: bash

        $ ./phylociraptor align --verbose -t local=4
        $ ./phylociraptor filter-align --verbose -t local=4
        $ ./phylociraptor modeltest --verbose -t local=4
        $ ./phylociraptor mltree --verbose -t local=4
        $ ./phylociraptor speciestree --verbose -t local=4

.. note::
        
        While this is running, or at any time, really, you can get a quick overview of the current progress of you analyses by running the following command (open a new terminal window and navigate to your working directory first, in case you want to check while a step is in progress):


.. code-block:: bash

        	$ ./phylociraptor check --verbose


We could also infer phylogenomic trees with *RAxML-NG* - we'd just need to enable it in the *settings file* and (re-)run ``phylociraptor mltree``.

.. note::

	Running *RAxML-NG* for all datasets would actually take about 3 hours (if you stick to the 4 computational threads) and this is entirely optional. You can skip this step, and move straight to the next section, unless you have the time.


Open ``data/config.yaml`` in your favourite text editor and search the section ``mltree:``, and change:

.. code-block:: bash

        mltree:
	    method: ["iqtree"]

to:

.. code-block:: bash

        mltree:
	    method: ["iqtree", "raxml"]

Now, let's see what would happen if you reran the mltree inference of phylociraptor after this change. Note, that we add ``--dry`` to the command, which will result in a so-called *dry-run*, i.e. don't actually execute, but just show me what would happen.

.. code-block:: bash

        $ ./phylociraptor mltree --verbose -t local=4 --dry

As expected, phylociraptor would prepare the data and run raxml for the 27 aligner/trimmer/bootstrap cutoff combinations that we have already completed with *IQ-Tree*. 

If you are happy to get into this, let's do this (runtime: 180 minutes):

.. code-block:: bash

        $ ./phylociraptor mltree --verbose -t local=4



----------------------------------------
Exploration of results and post-processing
----------------------------------------

Now, let's do some post-processing. Plot a few trees, evaluate conflicts between trees, etc.

First, let's generate an updated report.

.. code-block:: bash

        $ ./phylociraptor report

We can also plot the gist of the analyses as PDF.

.. code-block:: bash

        $ ./phylociraptor report --figure

If you're interested in inspecting the phylogenomic trees that have been inferred you can copy the textual representation of a given tree from the report and visualise it e.g. with [IToL](https://itol.embl.de/upload.cgi).

Phylociraptor also has a utility to plot trees to PDFs. Let's try. The random number seed in this case controls the colors of the tips in the tree, so if you keep using the same seed a parituclar taxon should always be displayed in the same color, even if the topology is different.

.. code-block:: bash

	$ ./phylociraptor util plot-tree -i results/phylogeny/iqtree/bootstrap-cutoff-70/clustalo-trimal.727b788ba6/concat.treefile -g Neurospora_crassa,Usnea_hakonensis --seed 42

This will produce a PDF with the name ``iqtree-clustalo-trimal-70-none-tree.pdf``. 

If the sample names in the *data file* are actually valid species binomials you can annotate the tree with taxonomic information. First, let's query Genbank for the taxonomic information for the taxa included in our analyses.

.. code-block:: bash

	$ ./phylociraptor util get-lineage -d data/test_cases/minimal/minimal.csv --outfile lineage_information.txt --force

Then annotate the tree with the taxonomic information at the level *class*.

.. code-block:: bash

	$ ./phylociraptor util plot-tree -i results/phylogeny/iqtree/bootstrap-cutoff-70/clustalo-trimal.727b788ba6/concat.treefile -g Neurospora_crassa,Usnea_hakonensis --seed 42 -l lineage_information.txt -e class 

Check out the PDF ``iqtree-clustalo-trimal-70-none-tree.pdf``.

Let's estimate topological conflict for all possible pairs of trees that we've inferred. For a given pair of trees, this is done by drawing quartets of tips from the first tree and check whether this paricular quartet is present in the second tree. Note, that the number of quartets that can be drawn from a given tree increases drastically with the number of tips in the overall tree and in large trees sampling all possibilities may be very time consuming. Therefore we tell our tool to stop the sampling if each tip in the tree was incorporated (on average) in a certain number of radomly drawn unique quartets, say 200. This will not be possible in the current dataset - our tree has just 6 tips - so in this case *estimate-conflict* will actually sample all quartets occuring in our trees and then stop. The proportion of quartets not shared between two trees is taken as an estimate for topological conflict.  We'll do this using 4 computational threads. Output will be written to text files with the prefix ``quartet-stop-200.*``. 

.. note::

	Note that you may specify (optionally) a seed for the random number generator. If you omit this, *estimate-conflict* will pick a random seed and report it to you for future reference and reproducibility. 

.. code-block:: bash

	$ ./phylociraptor util estimate-conflict -o quartet-stop-200 --stopby tipcoverage=200 -t 4 --seed 42


Now, let's plot the results of our conflict estimation to a heatmap.

.. code-block:: bash

	$ ./phylociraptor util plot-heatmap -m quartet-stop-200.similarity_matrix.csv -r quartet-stop-200.treelist.tsv

Let's inspect the result. A value of 1.00 means that, in a given pairwise comparison, 100% of the quartets that were sampled were identical, i.e. present in both trees. 0.96 would indicate that 4% of the sampled quartets were not shared in the pair of trees.

We can also visualise the conflict between two trees (in this case tree T1 vs. T4, see ``quartet-stop-200.treelist.tsv``), given the conflicting quartets that we've inferred.

.. code-block:: bash

	$ ./phylociraptor util plot-conflict -i T1,T4 -l lineage_information.txt -e class -r quartet-stop-200.treelist.tsv -q quartet-stop-200.quartets.csv --seed 42


Another way of summarizing differences between trees is based on a Principal component analysis (PCA) of all tip-2-tip distances of all taxa in each tree. Again we have a utility to calculate and visualise this.

.. code-block:: bash

	$ ./phylociraptor util plot-pca -r quartet-stop-200.treelist.tsv --seed 123 -t 5

The resulting file will be called ``tip2tip-PCA-all.pdf``. Rename it, please, since the next command will overwrite it.

The trees to be included in the comparision are specified in the treelist file. Let's say we'll only want to inspect differences in trees inferred with *ASTRAL*.

.. code-block:: bash

	$ grep "astral" quartet-stop-200.treelist.tsv > quartet-stop-200.treelist.astral.tsv
	$ ./phylociraptor util plot-pca -r quartet-stop-200.treelist.astral.tsv --seed 123 -t 5

The resulting file will be called ``tip2tip-PCA-all.pdf``. Rename it, please, since the next command will overwrite it.

We can also exclude trees inferred with *ASTRAL* and consequently only include maximum-likelihood based trees.

.. code-block:: bash

	$ grep -v "astral" quartet-stop-200.treelist.tsv > quartet-stop-200.treelist.ml.tsv
	$ ./phylociraptor util plot-pca -r quartet-stop-200.treelist.ml.tsv --seed 123 -t 5


The resulting file will be called ``tip2tip-PCA-all.pdf``.

Now we can investigate whether aligner, trimmer or bootstrap cutoffs may have a consistent effect on branch lengths in our trees.


----------------------------------------
A fond farewell
----------------------------------------

We hope you had fun! Give us a shout with feedback! Thanks for joining us today!
