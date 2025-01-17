.. _getting_help-knownproblems::

.. role:: bash(code)
   :language: bash

================================
Known problems
================================

----------------------------------------------------------------
ERROR: singularity image is not in an allowed configured path
----------------------------------------------------------------

This can happen in environments where executing singularity containers is restricted to certain places in the file system.
This is a feature of singularity to make it more secure. It is not a bug of phylociraptor. See `here <https://github.com/hpcng/singularity/issues/458>`_ for details.
To resolve this you need to first figure out where (in the file system) singularity images can be executed. To do this you can run:

.. code:: bash

   $ grep '^limit container paths' /etc/singularity/singularity.conf 
   limit container paths = /scratch, /global, /home


This shows that in this environment singularity images can only be executed when placed in :bash:`/scratch`, :bash:`/global` or :bash:`/home` (including sub directories).
Make sure you place phylociraptor in one of those directories and it should work. 
The file :bash:`/etc/singularity/singularity.conf` may be located in a different place depending on your computational environment.

---------------------------------------------------------------------
Singularity Error: Disk quota exceeded or FATAL: could not open image
---------------------------------------------------------------------

This can happen if there is not enough space in your :bash:`~/.singularity` directory, which is the usual place where singularity puts newly pull container images.
A possible way to solve this is to run this command before running phylociraptor for the first time.

.. code:: bash

	export SINGULARITY_DISABLE_CACHE=true


This command disables the standard storage location of newly pulled singularity containers and uses a system wide tmp directory instead.


----------------------------------------------
Singularity is not installed on my HPC Cluster
----------------------------------------------

Please get in touch with your local HPC support team and ask them to install Singularity for you.

-----------------------------------------------------
ImportError: Unable to import required dependencies:
-----------------------------------------------------

This error indicates that there is something wrong with the python library path. We have seen this happen when snakemake is not installed through conda but loaded as an environment module instead. To fix this problem, please install snakemake into it's own conda environment.

----------------------------------------------------------------------
My jobs fails without error message or log file entry on a HPC cluster
----------------------------------------------------------------------

It can be tricky to diagnose why a job has failed when there is no error message or log file entry. First, it is generally a good idea to doublecheck the log files of the failed job to make sure there is really no Error message in the log file. Please also check the .out files, because some processes will write error messages to stdout and they will end up in the .out rather than the .err files.

We have noticed that the cause of a silently failing job is often lack of memory. One possible solution is therefore to increase the requested amount of memory in the cluster config file. We have seen this work in several occasions e.g. when creating gene-trees for a subsequent species tree.

---------------------------
Error running slurm prolog:
---------------------------

This means that SLURM encountered a problem executing SLURM internal scripts prior to the job execution. There can be several different error codes associated with this problem and it may be necessary to contact your HPC support team to resolve what specifically goes wrong. In our experience it can already be enough to simply resubmit the job for the error to disappear.

-------------------------------------------------------------------
When running phylociraptor I get the error IncompleteFilesException
-------------------------------------------------------------------

An example of how such an error can look like is this:

.. code:: bash
        IncompleteFilesException:
        The files below seem to be incomplete. If you are sure that certain files are not incomplete, mark them as complete with

        snakemake --cleanup-metadata <filenames>

        To re-generate the files rerun your command with the --rerun-incomplete flag.   
        Incomplete files:
        results/checkpoints/gene_trees/mafft-trimal/EOG092C12O0_genetree.done
        results/phylogeny/gene_trees/mafft-trimal/EOG092C12O0/EOG092C12O0_gt.treefile
        results/checkpoints/gene_trees/mafft-trimal/EOG092C0JYF_genetree.done
        results/phylogeny/gene_trees/mafft-trimal/EOG092C0JYF/EOG092C0JYF_gt.treefile
        
Typically this has to do with several individual jobs failing on a HPC cluster without proper error message. We have seen this happening in cases where the job would need more RAM than what was specified in the cluster config file.
You may want to check the log files of the respective jobs to see what the case of the failure was, however we saw several cases where the cluster would not write the reason of the failure to the log file.

There are several possible ways to solve this:

1. Make sure the requested resources match the requirements for your jobs and dataset and change the cluster config file accordingly.

2. Run the step of the pipeline again and add :bash:`--snakemake="--rerun-incomplete"` to your phylociraptor command.

3. As a last resort you can also manually delete all the files listed as incomplete and run the pipeline again.

--------------------------------------------
Segmentation fault (core dumped) in raxml-ng
--------------------------------------------

While this error can have different causes it is known that raxmlng can have problems with alignments that have many more taxa than sites in the alignment. See `here <https://github.com/amkozlov/raxml-ng/issues/122>`_ on the official raxmlng Github Page. In phylociraptor this situation may occur if you have a large dataset while imposing a very high average bootstrap cutoff value. This can lead to only very few alignments for tree calculation while the number of taxa in the alignments is still very high. A solution might be to reduce the mean bootstrap cutoff value, so that more alignments are available for tree calculation. It is a good idea to create a report after the modeltesting step has been run to visualize possible cutoff values.


---------------------------------------
Several of my mafft alignment jobs fail
---------------------------------------

We have encountered this when many long amino-acid sequences should be aligned. In such a case mafft can become quite memory hungry. When phylociraptor is run on a HPC cluster and mafft reaches the memory limit on a node the job will crash. A workaround is to add the flag :bash:`--memsave` to the mafft options.

------------------------------------------
My Maximum Likelihood trees don't complete
------------------------------------------

It is usually a good idea to check the respective log files in the :bash:`log` directory. On HPC systems with batch job submission there is usually a wall-time limit preventing jobs to run only for a maximum amount of time.
Especially with larg phylogenomic datasets this limit can be reached quite easily. The solution to this problem is specific to the cluster environment and you may have to contact your local HPC support team for advice on how to extend the wall time limit.

On SLURM based systems another solution is to add the line :bash:`dependency: "singleton"` to the iqtree or raxmlng (depending on what you use) section in your cluster config file. Next you can run the phylociraptor mltree step several times. The singleton dependency setting will tell SLURM to create a chain of jobs that have the same name. Look `here <https://slurm.schedmd.com/sbatch.html>`_ for additional details.
 

