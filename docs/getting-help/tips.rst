.. _getting_help-faqs::

===============
Tips and Tricks
===============

On this page you can find short answers to frequently asked questions about phylociraptor.

--------------------------
What is phylociraptor?
--------------------------

Phylociraptor is a bioinformatics pipeline to create phylogenomic trees for a specified set of species using different alignment, trimming and tree reconstruction methods. For more info look `here <../introduction/about.html>`_ . 

-----------------------------------
Where can I download phylociraptor?
-----------------------------------

Phylocirpator is freely available on GitHub. You can get it `here <https://github.com/reslp/phylociraptor>`_ .

-------------------------------------------------------------
Why have some of the genomes I specified not been downloaded?
-------------------------------------------------------------

A failed genome download can happen due to serveral reasons. phylociraptor will always check if the taxon id of the specified species and the available genomes match. However there are cases in which they will not:

1. The species name was misspelled
2. The genome is deposited under a different name: If the genome was deposited under an infraspecific name say a strain oder subspecies, you have to provide the full name including infraspecific names.
3. There is a problem with the connection to NCBI.
4. No genome is available for the specified species.

If some of the species you specified are not present, you can always look into the log/genome_download.log file to see what has happened.

----------------------------------------------
Controlling the number of threads in -t serial
----------------------------------------------

By default phylociraptor will use the maximum number of available threads when it is run in :bash:`-t serial` mode. If you would like to change this you can use :bash:`-t serial=4` to restrict it to use only four threads for eaxmple. This only applies to parts which alread use multithreading, such as the tree calculation steps.

