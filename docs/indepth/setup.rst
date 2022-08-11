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
.. _snakemake: https://snakemake.github.io/
.. _singularity: https://sylabs.io/
.. _docker: https://docker.com/

 
=====================
Setup phylocirpator
=====================

----------------
Requirements
----------------

Phylocirpator only has two requirements: 

* `snakemake`_ 6.0.2

* `singularity <https://sylabs.io/>`_ 3.4.1+ or `docker`_ 

All other software used is prepackaged in containers and no installation is needed.

----------------
Installation
----------------

There are three simple options how to install phylocirator:

1. Just download and unpack the latest phylociraptor release.

2. Download the latest development version of phylociraptor using :bash:`git clone`

.. code:: bash

   git clone --recursive git@github.com:reslp/phylociraptor.git

or

.. code:: bash

   git clone --recursive https://github.com/reslp/phylociraptor.git 

3. If you don't have git, you can download the latest version also as ZIP file:

.. code:: bash

   wget https://github.com/reslp/phylociraptor/archive/master.zip
   mv master.zip phylociraptor.zip
   unzip phylociraptor.zip

------------------------------------
How phylociraptor downloads genomes
------------------------------------

After you have downloaded phylociraptor and edited the config files to fit your analysis it is time to run :bash:`phylociraptor setup`. During this step several things happen: First, phylocirpator will read the samples CSV file and start to download genomes in case you have specified them in the config file. 
This is done with a custom python script, to interact with the NCBI servers. Our script includes several fail-safe checks to ensure that the downloaded genome is correct and belongs to the desired species (by determining the NCBI taxon ID of the specified species name in the config file).
In case there are multiple genomes for a species, it will decide which genome to download based on a set of rules implemented in the script: If a genome is labeled "reference" it will download this genomes, if no reference genome is available it will insted look for a "representative genome" (look `here <https://support.nlm.nih.gov/knowledgebase/article/KA-03578/en-us>`_ to learn more about NCBI genome types).

If also no representative genome is available, it will download the latest genome based on publication date of the genome. So that this works properly it is important to **specify the correct name under which a genome is deposited at NCBI**. If you follow the provided instructions above properly, everything should work as expected.  

In some cases phylocirpator may be unable to find a genome. This can be due to several reasons:

1. The species name is mispelled in the config file.
2. The specified species name is ambiguous (eg. when you specify only a genus name and there are genomes of multiple species in this genus available).
3. The deposited genome has a different taxon ID. This can happen if there is no genome of the specified organisms but there is a genome of a symbiont or in another way associated organism which does not have a proper species name.

Downloading many genomes from NCBI will take some time (expect this to take several hours if you have >1000 genomes), however it is not particularly CPU intensive. Make sure you have a stable internet connection though. After genome download is finished phylociraptor will rename all assemblies and place symlinks of each genome it was able to gather into the folder :bash:`results/assemblies`. The files in this folder is what phylocirpator will use for subsequent analysis steps.
In general it is a good idea to run :bash:`./phylociraptor report` after genome download. This will give you an overview of the downloaded genomes and which genomes had a failed download. For more information on the report go here.

-------------------------------
Preparing BUSCO and Augustus
-------------------------------

During setup, phylociraptor will download the specified BUSCO set and place it in the folder :bash:`results/busco_set`. BUSCO relies on `Augustus`_ to predict genes in the genomes. Augustus (and thus BUSCO) uses pre-trained models to improve prediction accuracy. These models are stored in a special directory, to which phylociraptor needs read/write access. We use `software containers <https://www.docker.com/resources/what-container>`_ to ensure maximum portability of phylociraptor. However the used software container system on HPC clusters `Singularity <https://sylabs.io/docs/>`_ does not usually allow write access to directories inside a container. We therefore create a copy of the directory we need to write to inside the results directory. This directory is the Augustus config folder which is located at :bash:`results/augustus_config_path`.

This approach has another benefit for you: If you want you can use **your own pretrained models** which you may have from previous Augustus runs. Just copy your models folder to :bash:`results/augustus_config_path/species` and specify you pretrained species in the config.yaml file.
