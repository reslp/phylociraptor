.. role:: bash(code)
   :language: bash


=====================================================
Different ways to calculate Maximum-Likelihood trees
=====================================================

Maximum-likelihood trees can be calculated in a variety of ways. The basis typically forms an alignment which can be a supermatrix (concatenated) alignment combining the information of many genes. This alignment is used to infer and optimze parameters of a model which in turn aims to explains the evolutionary history of the included samples to produce a tree. How these models are applied, and which models are used differs widely in the published literature. For example it is possible infer the best fitting evolutionary model (under a specific optimality criterion) for each gene and apply these models in the subsequent maximum-likelihood analysis. It would also be possible to apply a single model to the whole alignment. Or it would be possible to specify a specific model which should be applied to each gene.  


Additionally the included genes could be filtered before combining them to a concatened alignment. Phylociraptor tries to capture the full breadth of approaches typically applied in large-scale phylogenomic analyses. However the number of differnt combinations can be overwhelming. The purpose of this document is therefore to explain how some of the common scenarios to perform maximum-likelihood analyses taking into account different ways to specify models and apply gene filtering can be set up in phylociraptor. 


-------------------------------------------------------
Creating an unpartitioned analysis using a single model
-------------------------------------------------------

The simples way would be to apply a single model to the complete concatenated alignment. This step **does not** require :bash:`phylociraptor modeltest`. Have a look at the relevant section of the config file:

.. code-block:: bash
   :linenos:
   :emphasize-lines: 2,5,17

    genetree_filtering:
    bootstrap_cutoff: [0]

    mltree:
        method: ["iqtree-unpartitioned"]
        threads:
            iqtree: 20
            raxml: 20
            iqtree-unpartitioned: 20
        bootstrap:
            iqtree: 1000
            raxml: 100
            iqtree-unpartitioned: 1000
        options:
            iqtree: ""
            raxml: ""
            iqtree-unpartitioned: "-m C10"

The relevant settings are:   

|Line 2: Make sure to include 0 in the values used for the mean bootstrap cutoff. This is necessary when :bash:`phylociraptor modeltest` has not been run.   
|Line 5: Set the method to :bash:`iqtree-unpartitioned`.   
|Line 17: Specify the desired substitution model in iqtrees syntax using :bash:`-m`. In this case :bash:`-m C10` for a CAT based mixture model.   

.. note::

   You can read more about the available models in IQ-Tree `here <http://iqtree.org/doc/Substitution-Models>`_.


-------------------------------------------------------------------
Creating an unpartitioned analysis while testing for the best model
-------------------------------------------------------------------

This works almost like the above example.


.. code-block:: bash
   :linenos:
   :emphasize-lines: 17

    genetree_filtering:
    bootstrap_cutoff: [0]

    mltree:
        method: ["iqtree-unpartitioned"]
        threads:
            iqtree: 20
            raxml: 20
            iqtree-unpartitioned: 20
        bootstrap:
            iqtree: 1000
            raxml: 100
            iqtree-unpartitioned: 1000
        options:
            iqtree: ""
            raxml: ""
            iqtree-unpartitioned: ""

The only difference to the above example is in line 17. When nothing is specified in this line, phylociraptor will default to :bash:`-m MFP` which includes a full IQ-Tree modeltest.


-----------------------------------------------------------------------------------------------
Creating an unpartitioned analysis while only using genes with a certain mean boostrap support.
-----------------------------------------------------------------------------------------------


.. note::

   This requires that :bash:`phylociraptor modeltest` was run before.


.. code-block:: bash
   :linenos:
   :emphasize-lines: 2,5,17

    genetree_filtering:
    bootstrap_cutoff: [50, 60, 70]

    mltree:
        method: ["iqtree-unpartitioned"]
        threads:
            iqtree: 20
            raxml: 20
            iqtree-unpartitioned: 20
        bootstrap:
            iqtree: 1000
            raxml: 100
            iqtree-unpartitioned: 1000
        options:
            iqtree: ""
            raxml: ""
            iqtree-unpartitioned: "-m C10"


The relevant settings are:   

|Line 2: Specify which mean bootstrap values you would like to use. 
|Line 5: Set the method to :bash:`iqtree-unpartitioned`.   
|Line 17: Specify the desired substitution model in iqtrees syntax using :bash:`-m`. In this case :bash:`-m C10` for a CAT based mixture model.   


.. note::

   Again you could also leave the setting in line 17 blank in which case phylociraptor will use :bash:`-m MFP`(full IQ-Tree modeltest) as default.


-----------------------------------------------------------------------------------------------
Running a partitioned analysis with a fixed model.
-----------------------------------------------------------------------------------------------


.. note::

   This requires that :bash:`phylociraptor modeltest` was run before.


.. code-block:: bash
   :linenos:
   :emphasize-lines: 2,5,15

    genetree_filtering:
    bootstrap_cutoff: [50, 60, 70, 80]

    mltree:
        method: ["iqtree"]
        threads:
            iqtree: 20
            raxml: 20
            iqtree-unpartitioned: 20
        bootstrap:
            iqtree: 1000
            raxml: 100
            iqtree-unpartitioned: 1000
        options:
            iqtree: "-m C10"
            raxml: ""
            iqtree-unpartitioned: ""


The relevant settings are:   

|Line 2: Specify which mean bootstrap values you would like to use. 
|Line 5: Set the method to :bash:`iqtree`.   
|Line 15: Specify the desired substitution model in iqtrees syntax using :bash:`-m`. In this case :bash:`-m C10` for a CAT based mixture model.   


.. note::

   Note that the behavior here is different to an unpartitioned analysis. When you leave the option in line 15 blank, this will default to using the best models from :bash:`phylociraptor modeltest`.


