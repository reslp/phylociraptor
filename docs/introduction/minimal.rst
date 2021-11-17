.. role:: bash(code)
    :language: bash

============================
Tutorial - A minimal example
============================

This document guides you through a first example on how to run phylociraptor. This minmal dataset consists of 6 fungal genomes, three belonging to the Ascomycota and three to
Basidiomycota which will be downloaded during the first step of the run. Downstream analyses are carried out on 10 single-copy orthologous genes which will be identified in each genome with a modified BUSCO set.
The so recovered genes will be aligned and trimmed and we will calculate single-gene trees for each gene. After this, we will perform a Maximum-Likelihood analyses on the full dataset and infer a species tree. The whole analysis should take about 1 hour to complete.

------------------------
Dataset and prepatations
------------------------

Input files for this dataset are located in :bash:`data/test_cases/minimal`. To prepare the analysis execute the following command in your phylociraptore base directory:

.. code-block:: bash
        
        $ cp data/test_cases/minimal/config.yaml data/
        $ cp data/test_cases/minimal/minimal.csv data/
        $ cp data/test_cases/minimal/modify_busco.sh .


After that, you can check the :bash:`data/minmal.csv` file to see which species will be included in the analysis or you can view the :bash:`data/config.yaml` file to see
the parameters for each step.

-----------------------
Analysis setup
-----------------------

At first it is necessary to download the genomes which should be analyzed. You can do this with phylociraptor using this command:

.. code-block:: bash

        $ ./phylociraptor setup --verbose -t serial=1

This command runs in serial mode (without job submission to an HPC job scheduling system) and should take about 5 minutes to complete.


.. note::
        
        The flag :bash:`--verbose` in the above command and all the commands below can also be omitted. It is used here to viualize pipeline output on screen.

After this step is finished we will limit the search of orthologous genes to only 10 using a bash script. This is done mainly to reduce computation time. To reduce the search for orthologous genes to 10 run this command:

.. code-block:: bash

        $ ./modify_busco.sh

The downloaded genome assemblies are in :bash:`results/assemblies`. Metadata of downloaded genomes is in :bash:`results/downloaded_genomes`.

--------------------------------------------------
Identify single-copy orthologs
--------------------------------------------------

To infer single-copy orthologs in the six provided genomes run this command:

.. code-block:: bash
        
        $ ./phylociraptor orthology --verbose -t serial=4


.. note::
        
        We will run this command using four threads as indicated with :bash:`serial=4`.

This should take about 5 minutes to complete. Results of this step can be found in :bash:`results/orthology`.


--------------------------------
Filter orthology results
--------------------------------

This step is used to filter orthology results based on how many orthologs where found in each genome. To run execute this command:

.. code-block:: bash
        
        $ ./phylociraptor filter-orthology --verbose -t serial=1


.. note::
       
        In this example nothing will be filtered since the ten selected genes should be present in all genomes.

This step should take less than a minute to complete.

------------------------------
Align single copy orthologs
------------------------------

In this step we will create alignments of the single-copy orthologs recovered in each genome using mafft and clustalo.


.. code-block:: bash

        $ ./phylociraptor align --verbose -t serial=1

This step should take about three minutes to complete. Results will be in :bash:`results/alignments`.


------------------------------
Filter alignments
------------------------------

Now the produced alignments will be trimmed using Aliscore/Alicut annd trimAl to remove poorly aligned regions. This is done for each aligner and trimmer combination, so after this step you will have 4 sets of alignments.

 
.. code-block:: bash

        $ ./phylociraptor filter-align --verbose -t serial=1

Results of this step will be in :bash:`results/alignments`. It will take about four minutes for this step to complete.


----------------------------------
Infer the best substitution models
----------------------------------

Using filtered alignments we will now infer the best substitution model for each alignment and also calculate single gene-trees. This is done for each alignment in every aligner and trimmer combination.


.. code-block:: bash

        $ ./phylociraptor modeltest --verbose -t serial=4

This step should finish in about five minutes and results will be located in :bash:`results/modeltest`

----------------------------------------
Calculate a full Maximum-Likelihood tree
----------------------------------------

Now it is time to calculate full (concatenated) Maximum-Likelihood trees. We will use iqtree and raxml-ng in this step and infer trees for every aligner and trimmer combination. The analysis will be partitioned using the best substitution models inferred with the step above. Addditionally phylociraptor will take only gene-trees above specified bootstrap values. In this example the bootstrap values used are 50, 60 and 70.


.. code-block:: bash

        $ ./phylociraptor mltree --verbose -t serial=4

This step takes about 30 minutes to complete. Results are in :bash:`results/phylogeny-50 results/phylogeny-60 results/phylogeny-70`.


----------------------------------------
Calculate a species tree
----------------------------------------

Finally we will calculate species trees for every aligner and trimmer combination and every bootstrap cutoff value.


.. code-block:: bash

        $ ./phylociraptor speciestree --verbose -t serial=1

This takes about a minute to complete. Results are in :bash:`results/phylogeny-50 results/phylogeny-60 results/phylogeny-70`.


After this is done we will have:

- 6 downloaded fungal genome assemblies
- 10 identified single-copy orthologs in each genome
- 20 multiple sequence alignments using mafft and clustal-o
- 38 trimmed alignments using trimal and aliscore/alicut (2 will be excluded during trimming)
- 38 maximum-likelihood gene trees calculated with iqtree
- 12 concatenated maximum-likelihood phylogenies using 3 different bootstrap-cutoff values
- 12 species trees using 3 different bootstrap-cutoff values



