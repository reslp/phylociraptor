
.. role:: bash(code)
   :language: bash

=============================================
A posteriori analyses with phylociraptor util
=============================================

Phylociraptor can create lots of trees quickly. Which tree is the correct one? There is no good answer to this question, however with functions provided by :bash:`phylociraptor util` it is possible to compare similarity of trees and impose lineage information from NCBI on tree plots.

What is examplained here requires phylociraptor to have successfully produced trees. The provided utilities aim to give an overview about tree similarity and should help to identify problematic or unstable areas in the trees. However they are not thorough comparisons of trees. There are great dedicated software packages for this task. For example: https://ms609.github.io/TreeDist/. 

---------------------------------------------------
Download taxonomic lineage information from NCBI
---------------------------------------------------

It is possible to download lineage information with :bash:`phylociraptor util get-lineage` for all previously downloaded genomes. You can look at the options for this command with :bash:`./phylociraptor util get-lineage -h`.

The output of this command can later be used when plotting trees and to estimate which taxonomic groups have been recovered as monophyletic in different trees.

Here is a full example command:

.. code-block:: bash

    $ ./phylociraptor util get-lineage --genomefile results/statistics/downloaded_genomes_statistics.txt --outfile lineage_information.txt

-----------------------------------------------------------------
Estimate conflicts between trees based on subtrees with four tips
-----------------------------------------------------------------

It can be difficult to compare topological differences between large phylogenetic trees and there are several reasons for that. Trees can differ in the sets of included taxa and many commonly used metrics reduce tree-similarity to a single value. While there is nothing wrong with this approach, it would also be nice to be able to compare where in the tree conflicts are located. Here we adopt a strategy to be able to estimate conflict between trees despite possible differences in the number of tips while also being able to get an idea where the conflicting positions of the tree are. The method is based on comparing trees based on quartets of tips.

Here are a few example commands:

.. code-block:: bash

	$ ./phylociraptor util estimate-conflict --outprefix quartet-5000 -n 5000
	$ ./phylociraptor util estimate-conflict --outprefix quartet-100 --stopby conflicts=100 --seed 120
	$ ./phylociraptor util estimate-conflict --outprefix quartets-conflict --stopby tipcoverage=200 -l lineages.txt --selecttaxa order=Diperta,Hymenoptera -t 20

The first command will estimate conflicts based on 5000 random quartets. Output files will have the prefix quartet-5000.

The second command will estimate conflicts until 100 conflicting quartets have been identified. It also uses a random seed so that the same quartets will be compared during each run.

The third example samples quartets till each tip has been included in quartets at least 200 times. It also uses the lineage information file and only calculates quartets from the orders Diperta and Hymenoptera. It runs this analysis using 20 threads.

-----------------------
Plot phylogenomic trees
-----------------------

With :bash:`phylociraptor util plot-tree` you can visualize your trees. This utility output PDF files with tree images.

Here are some example commands:

.. code-block:: bash

	$ ./phylociraptor util plot-tree
	$ ./phylociraptor util plot-tree -i results/phylogeny-75/iqtree/mafft-aliscore/concat.contree -l lineage_arthropoda.txt -e family

The first command is the simplest plotting command. It create plots for all trees it finds in the :bash:`results/` directory.

The second command will create only a single plot for the iqtree tree with mafft and aliscore alignments. It will highlight different families in color by using the information from the lineage file. When a lineage file and a taxonomic level are specified then a second plot will be produced indicating if the different groups are recovered as monophyletic in the tree.



-----------------------------------------------
Plotting conflicts based on quartet differences
-----------------------------------------------

With this utility it is possible to get a visual idea about conflicts between to trees based on quartets. The plotted tree will look similar to trees plotted with :bash:`phylociraptor util ploz-tree`. The difference is that in this case two trees are plotted instead of one and the conflicts between them are highlighted by thicker branches. The thickness of a branch indicates in how many conflicting quartets this branch is included in. Thicker branches thus mean more conflict.

Again, here are several example commands:

.. code-block:: bash
	
	$ ./phylociraptor util plot-conflict -i T1,T10 --quartetfile quartet-5000.quartets.csv -r quartet-5000.treelist.tsv
	$ ./phylociraptor util plot-conflict -i T2,T23 --quartetfile quartet-5000.quartets.csv -r quartet-5000.treelist.tsv --seed 123 --outgroup Ramazzottius_varieornatus,Hypsibius_dujardini

The first example compares trees 1 and 10. The two parameters --quartetfile and -r are mandatory.

The second command compares trees 2 and 23. It uses a random seed for reproducibility and the trees will be rooted with Ramazzottius_varieornatus and Hypsibius_dujardini before they are plotted.

--------------------------------------------------------------------------------
Plotting similarity of trees based on quartet occurences and tip 2 tip distances
--------------------------------------------------------------------------------

This creates an overview heatmap of quartet conflicts. Additionally this utility can also calculate tip2tip distances for all trees and plot the differences as an annotated PCA for quick comparison of trees.

A few example commands should make it more clear:

.. code-block:: bash
	
	$ ./phylociraptor util plot-similarity -m quartet-5000.similarity_matrix.csv 
	$ ./phylociraptor util plot-similarity -m quartet-5000.similarity_matrix.csv -r quartet-5000.treelist.tsv -s 123 --ndistances 100 -t 8

The first command will only create a heatmap of overall similarity of quartets in percent between trees.

The second command will create the heatmap but will also compute tip2tip distances for 100 pairs of tips using eight threads and a random seed 123. These tip2tip distance matrices will be plotted as a PCA.

-------------------------
Modify the used BUSCO set
-------------------------

It is possible to reduce the downloaded BUSCO set to only a specified number of genes using ``phylociraptor util modify-busco``. This could be helpful if you would like to make quick tests with only a smaller set of genes before running a full analysis.

.. warning::

	For this to work you need to first run ``phylociraptor setup`` so that the specified BUSCO set is available in ``results/orthology/busco/busco_set``.


.. code-block:: bash

	$ ./phylociraptor --debug util modify-busco -g EOG092C3MKB,EOG092C59RA,EOG092C5GXJ -b fungi_odb9
	$ ./phylociraptor --debug util modify-busco -n 20 -b arthropoda_odb9 -s 42

The first command will reduce the the ``fungi_odb9`` set to the three specified genes ``EOG092C3MKB,EOG092C59RA,EOG092C5GXJ``.
The second command reduces the the ``arthropoda_odb9`` set to ``20`` random genes using the random seed ``42`` for reproducibility.

Phylociraptor will keep a backup of the original unmodified BUSCO set in ``results/orthology/busco/busco_set/<SETNAME>_backup``.


