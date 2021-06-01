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


This command disables the standard stoarge location of newly pulled singulariyt containers and uses a system wide tmp directory instead.


-----------------------------------------------------
ImportError: Unable to import required dependencies:
-----------------------------------------------------

This error indicates that there is something wrong with the python library path. We have seen this happen when snakemake is not installed through conda but loaded as an environment module instead. To fix this problem, please install snakemake into it's own conda environment.

