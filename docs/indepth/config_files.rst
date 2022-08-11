
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
 
================
The config files
================

Phylociraptor uses two main config files for setting analysis options:

* The :bash:`data/config.yaml` file, which controls settings of the running behavior of phylociraptor.
* The :bash:`data/sample.csv` file which contains information about the used samples.

Additionally there is a cluster config file, which control the behavior on HPC clusters with different job submission systems. We provide individual cluster config file for different submission systems.

---------------------
The config.yaml file
---------------------

When you download phylociraptor the file is called :bash:`config.yaml.template` , it is located in the :bash:`data/` directory and you have to first rename it to :bash:`config.yaml`. As you can see from the file extension, the file is in `YAML`_ format. Here is how the complete file looks like (without comments). We also added line numbers for this example, so that is easier to refer to different settings in this documentation.

.. literalinclude:: config.yaml.example
	:language: python
	:linenos:

**Line 1:** Here we specify the samples configuration file path (We will talk about it later). It is important that you give a relative path here. It is best to place the file in your :bash:`data` directory inside your phylociraptor folder.

**Lines 2-6:** Here you can specify the settings for the `BUSCO`_ runs (which will run individually for each genome specified in the file in line 1).  
**Line 3** refers to the BUSCO set path. Your you should give the web address to the BUSCO set. You can look `here <https://busco-archive.ezlab.org/v3/>`_ to see which sets are available. Currently phylociraptor uses BUSCO v3, so make sure you specify the correct download link in the config file.  
**Line 4** sets the number of threads which should be used for each BUSCO run. This will depend on the type and power of your computer environment. Increasing the number of threads can speed up the BUSCO analyses. Usually you can leave this at the standard setting of 8 threads.  

**Line 5** specifies which models `Augustus`_ should uses to predict genes as part of the BUSCO step. Augustus offers many different pretrained species models and it is generally a good idea to use a species which is closely related to the species you want to analyse. It will improve BUSCO performance and results in more (and more accurate) genes available for phylogenomic reconstruction.  

**Line 6** Here you can specify additional parameters which will be passed on to BUSCO. Check the `BUSCO`_ website and manual for possible options. We have not tested all possible parameters and we can not guarantee that all of them will work with phylociraptor.  

**Line 7-12:** This sections controls the filtering behaviour of Orthology inference and alignments. Filtering will take place at different steps in the analysis workflow (-m filter-orthology and -m filter-align). This allows to test the effect of different filtering parameters on the dataset. We will go into additional detail later.  

.. note::
 
   | :bash:`cutoff` should be between 0 (no filtering) and 1 (only keep samples with 100% BUSCO completeness). Both extremes are probably unrealistic.
   | :bash:`minsp` should be between 3 and n (n is the number of samples in your analysis).

**Line 8:** Occasionaly BUSCO will place more than one sequences in files for single-copy BUSCO genes. This flag allows to set the behaviour how phylociraptor handles these cases. Possible options are: "persample" which will remove the duplicated sequences only for the species in which it occurs. This means that for the affected sample, this gene will not be used. "pergene" will discard the whole gene for all samples. This gene will thus not be used for phylogenomic reconstruction at all.  

**Line 9:** The cutoff value sets the minimum proportion of BUSCO completeness (complete single-copy BUSCO genes) a sample needs to have to be included in downstream analysis. This flag allows to exclude samples with poor BUSCO performance eg. due to a suboptimal assembly. The default is 0.5 which will include samples which have at least 50% BUSCO completeness. Your should provide a value between 0 (no-cutoff) to 1 (only genomes with 100% BUSCO completness). Keep in mind that 1 is very unrealistic for most genomes!  

**Line 10:** Not all BUSCO genes will be detected with equal efficiency, thus your alignments for different BUSCO genes will contain a varying number of genes. minsp lets you set the minimum number of sequences need to be in an alignment for it to be included in the phylogenomic reconstruction. The default is 3, which means that at least 3 sequences need to be present for a gene to be included in the analysis. You should use a value between 3 and the total number of samples you have.  

**Line 11:** seq_type lets you specify whether you want to reconstruct trees based on nucleotide "nu" or amino-acid "aa" data. The default is "aa".  

**Line 12:** With min_parsimony_sites you can set how many parsimony informative sites an alignment needs to have to be included in the analysis. The default value is 10.  

**Line 13-16:** In this section you can change settings about the alignment step.

**Line 14:** You can specify the alignment method here. Currently we only support `mafft`_ but we plan support for additional aligners in the future.

**Line 15:** Specify the number of threads for the alignment step. As for BUSCO this depends on your computer environment. The alignment step is not the most time-consuming step, so a smaller number of threads should not impact performance too much. The default value is 8, which is usually fine.

**Line 16:** Here you can pass additional parameters to mafft eg. alignment methods. By default phylociraptor will tell mafft to determine the best alignment strategy automatically (--auto). For additional options see `mafft`_ website.

**Line 17-19:** Here you can set how to control the alignment trimming step, to filter out ambiguous position, misalignments and more. 

**Line 18:** This lets you control the trimming method. Currently `trimal`_ and alicut are the supported options here. To familiarize how these programs work it is a good idea to study the respective websites and manuals.  

**Line 19:** Because trimming can be done in many different ways, you can specify different trimming parameters here which will be passed on to the specified trimming method.  

**Line 20-33** This part of the config file lets you specify the (species)tree reconstruction methods and set different additional parameters if you wish.

**Line 21** Specify the maximum-likelihood tree reconstruction methods here. These should be under quotes and seperated by a space. Currently we support `iqtree`_ and `raxml-ng`_ (specify raxml in the config file). We plan to include additionla methods in the future.

**Line 23** Here you can specify the species tree reconstruction method. Currently we only support `astral`_ but we plan to include additional methods in the future.

**Line 24-29** These settings change the behaviour of iqtree.  

**Line 25:** Lets you change the number of threads for your iqtree run. The way this works is a bit different from other threading settings. For iqtree the number of threads will rever to the maximum number of threads. When phylociraptor is run an a HPC cluster with a job sheduler it will let the sheduler know how many threads should be requested for the iqtree job. Internally iqtree will then determine the optimal number of threads from 1 to number of specified threads (using -nt AUTO -ntmax {threads}).  

**Line 26:** Lets you specify the number of bootstrap replicates. This will use the Ultra-Fast bootstraping (-bb) option of iqtree.  

**Line 27:** Here you can specify the substitution model used by iqtree. The default setting is MFP which will determine the best substitution model for each gene partition as part of the iqtree run. Of you can also specify a specific model here as well, although it is very unrealistic to assume that the same model is optimal for all genes. In case you ran :bash:`phylociraptor -m modeltest` this setting will be ignored, since the best model for each genes was already inferred. To optimize runtimes we recommend to identify the best substitution models with :bash:`phylociraptor -m modeltest` before running :bash:`-m tree`. The modeltesting step is highly parallelized and probably a lot faster than model selection as part of the tree reconstruction.  

**Line 28:** You can also specify additional parameters for iqtree here. Keep in mind however that iqtree offers a myriad of different options and applications, most of which have not been tested with phylociraptor.  

**Line 29:** We provide an extra option to limit the amount of memory iqtree can use. This can be necessary for very large datasets on HPC clusters.  

**Line 31:** Specify the number of threads raxml should use here.  

**Line 32:** Specify the number of bootstrap replicates raxml will calculate here.  

**Line 33:** Here you can pass additional arguments to raxml. Similar to iqtree, raxml offers many different options. We have not tested them all for compatibility with phylociraptor.  

--------------------
The samples CSV file
--------------------

This is the second most important file phylociraptor needs to run and probably the most important file for you. In this file you specify which samples should be analyzed and whether they should be downloaded automatically or if you provide them as local files. Here is an example for a sample file which is part of phylociraptor when you clone it from github.

.. code-block:: bash
   :linenos:

   species,web_local,mode
   Salpingoeca rosetta,web=GCA_000188695.1,
   Coccidioides posadasii,web,
   Sclerophora sanguinea,data/assemblies/Sclerophora_sanguinea_Sclsan1_AssemblyScaffolds_Repeatmasked.fasta,
   Capsaspora owczarzaki,web=GCA_000151315.2,
   Dictyostelium lacteum,web=GCA_001606155.1,
   Paraphelidium tribonemae,data/assemblies/EP00158_Paraphelidium_tribonemae.fasta,protein
   Synchytrium microbalum,web=GCA_006535985.1,
   Nuclearia sp,data/assemblies/Nuclearia_sp_trinity_cdhit-0.95.fasta,transcriptome
   Stereomyxa ramosa,data/assemblies/Stereomyxa_ramosa_trinity_cdhit-0.95.fasta,transcriptome

Let us look at this example in more detail as it contains many different input formats line by line.

**Line 1** is the mandatory header line with the three mandatory columns seperated by comma (``,``).

**Line 2** specifies that a genome for the species ``Salpingoeca rosetta`` should be downloaded (indicated by ``web``). A specific assembly is specified which is ``GCA_000188695.1``. THlast column (``mode``) is empty which means that the default mode (which is ``genome``) should be applied during orthology search.

**Line 3** specifies the taxon ``Coccidioides posadasii`` for which an assembly will be downloaded from NCBI. Phylocirpator will decide which assembly to download in case there are multiple options.

**Line 4**:  for ``Sclerophora sanguinea`` we specify a locally provided genome assembly using a *relative* path.

.. tip::
   
   Only relative paths are supported here. We recommend placing all files fo locally provided assemblies and protein set in ``data/assemblies`` to avoid problems.

**Lines 5,6 and 8** are like **Line 2**. A specific genome assembly will be downloaded automatically.

**Line 7** specifies a locally provided set of predicted proteins for ``Paraphelidium tribonemae``. We need to indicated that this is a protein set in the column ``mode``.

**Line 9 and 10** specify locally provided transcriptome assemblies. We need to indicate that these are transcriptomes in the column ``mode``.

Using files from NCBI Genome Browser
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
 
The samples file is very similar to the TSV file you can download from the `NCBI Genome Browser`_ . This makes it easy to use files downloaded from the Genome Browser with phylociraptor. To make a TSV file you will to make only minimal changes:

1. Make sure the header line is compatible with phylocirpator and that the mandatory columns (``species web_local mode``) are present.

This is how it should look like, if you use the file from NCBI Genome Browser as basis:

.. code-block:: bash

	species,Organism_Groups,Size(Mb),Chromosomes,Organelles,Plasmids,Assemblies,web_local,mode

2. Add information to the last column to let phylocirpator know if a genomes should be downloaded :bash:`web` or is provided locally :bash:`local`.

If you have just downloaded a file from NCBI and you want to reformat the whole file to make it compatible with phylocirpator you can run this command. Make sure to check you header line after this again!

.. code-block:: bash

	perl -p00e 's/\n/,web\n/g' genomes.csv > genomes_added_web.csv

Now you can add your local files to the file at the end if you want. As said above, most of the information in the file is not necessary for phylociraptor, however to maintain the structure of the file you have to make sure to currecntly provide the fild seperators. This is how a line for a locally provided species should look like:

.. code-block:: bash

	Gouania hofrichteri,,,,,,,data/assemblies/Gouania_hofrichteri_4616STDY8352653_CONSENSUS_filtered20.fasta,

As you can see the last column contains the path to the assembly which should be used for this species. phylociraptor will realize that this is a locally provided file and not try to download it. This path has to be relative to the working directory of phylociraptor. To ensure that there are no problems with the paths to locally provided files it is a good idea to place the assemblies inside the :bash:`data/assemblies` directory.


--------------------------
HPC Cluster config files
--------------------------

phylociraptor can be easily scaled to run on large HPC clusters. This has the benefit that many individual tasks can run in parallel which can drastically decrease overall runtime of the analysis.
HPC clusters are usually very different with many individual (cluster specific) settings. Usually, individual jobs on HPC systems are managed by batch job-sheduling software which assigns resources (memory, CPUs and runtime) and starts jobs. Many different job schedulers exist.

We tested phylociraptor on clusters running `SLURM <https://slurm.schedmd.com/>`_ , `SGE <https://www.oracle.com/enterprise-manager/technologies/>`_ and `TORQUE <https://adaptivecomputing.com/cherry-services/torque-resource-manager/>`_.

When used on HPC clusters phylociraptor will automatically create and submit individual jobs for each task (there is no need to use batch submission files and commands such as :bash:`sbatch` or :bash:`qsub`). To make this work, phylociraptor needs to have some specific information on the underlying cluster. This information is provided through cluster configuration files in `YAML <https://en.wikipedia.org/wiki/YAML>`_ format.

Phylociraptor provides cluster configuration files templates for SLURM, SGE and TORQUE clusters. 

.. note::

	Cluster configuration files have to be altered to fit your own computer infrastructure. Usually only minor changes have to be made, most settings should work out of the box.


^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Adjusting the SLURM cluster config template
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Before you make any changes to the file, it is a good idea to make a copy:

.. code-block:: bash

	$ cp data/cluster-config-SLURM.yaml.template data/cluster-config-SLURM.yaml

Here are the first 30 lines of the config file for a SLURM cluster.

.. literalinclude:: cluster-config-SLURM.yaml
        :language: python
        :linenos:
	:emphasize-lines: 2,7,8

Typically it will only be necessary to change the highlighted lines in the example above. QOS (Quality of service) and partition are usually cluster specific and need to be changed. The walltime limit :bash:`time: "72:00:00"` will also depend on your cluster environment.

Make sure to check your cluster's documentation to find the correct values.

Other values (like :bash:`mem`, :bash:`ntasks`, etc. ) don't have to be changed. The provided values have been tested with differently sized datasets on different clusters and should work for the majority of datasets.



