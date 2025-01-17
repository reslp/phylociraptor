.. _getting_started-installation:

.. role:: bash(code)
    :language: bash

============
Installation
============

-------------
Requirements
-------------

Phylociraptor extensively uses software containers and thus has only a minimum number of dependencies:

* Linux or MacOS operating system
* globally installed singularity 3.4.1+ 
* installed snakemake 6.0.2+ (eg. in an anaconda environment)

On a HPC cluster (when utilizing parallel job execution):

* SGE, SLURM or TORQUE job scheduling system

.. warning::

	We strongly recommend that you use the snakemake version suggested above (6.0.2). Using other versions could lead to complications.


-------------------------
Obtaining phylocirator
-------------------------

Phylociraptor is available on GitHub. You can download it `here <https://github.com/reslp/phylociraptor>`_ .

The probably best way is to clone the repository directly using git (if available).

.. code-block:: bash

	$ git clone --recursive https://github.com/reslp/phylociraptor.git
	$ cd phylocirpator
	$ ./phylociraptor
	Welcome to phylociraptor, the rapid phylogenomic tree calculator
	
	Usage: phylociraptor <command> <arguments>
	
	Commands:
		setup			Setup pipeline
		orthology		Infer orthologs in a set of genomes
		filter-orthology	Filter orthology results
		align			Create alignments for orthologous genes
		filter-align		Trim and filter alignments
		model			Perform modeltesting
		mltree			Calculate Maximum-Likelihood phylogenomic trees
		speciestree		Calculate gene trees and species tree
		njtree			Calculate Neighbor-Joining tree
                bitree                  Calculate Bayesian-inference phylogenomic trees

		report			Create a HTML report of the run
                check                   Quickly check status of the run
		util			Utilities for a posteriori analyses of trees
	
		-v, --version 		Print version
		-h, --help		Display help
	
	Example:
		To see options for the setup step:
		./phylociraptor setup -h
	
		To run orthology inferrence for a set of genomes on a SLURM cluster:
		./phylociraptor orthology -t slurm -c data/cluster-config-SLURM.yaml
	
	
.. note::

    If you don't have git available, you can also download phylociraptor directly as ZIP file and unpack it to the desired location.

-------------------------------------------
Create a conda environment for snakemake
-------------------------------------------

.. warning::

	It is recommended to install snakemake using conda into it's own environment. Other options, such as loading snakemake as environment module can sometimes cause problems (eg. on HPC systems)

If you don't have conda installed, first look `here <https://docs.conda.io/en/latest/miniconda.html>`_ .

.. code-block:: bash

        $ conda install -n base -c conda-forge mamba
        $ mamba create -c conda-forge -c bioconda -n snakemake snakemake=6.0.2
	$ conda activate snakemake


When you run phylociraptor this environment needs to be activated.

Additional information on how to install snakemake can be found `here <https://snakemake.readthedocs.io/en/stable/getting_started/installation.html>`_ .

-------------------------------------------------------------------------------
Optional: Customize cluster configuration settings to fit your HPC environment
-------------------------------------------------------------------------------

phylociraptor will automatically submit jobs to SLURM, SGE and TORQUE job submission systems using :bash:`sbatch` or :bash:`qsub`. 

For this to work you will probably need to edit the correct cluster configuration file.

The files are :bash:`data/cluster-config-SLURM.yaml.template` and :bash:`data/cluster-config-SGE.yaml.template`. 







