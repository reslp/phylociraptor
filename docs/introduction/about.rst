
============
About
============

Phylociraptor is a bioinformatics pipeline to create **phylogenomic trees** for a specified set of species using different alignment, trimming and tree reconstruction methods. Installation is simple with only two dependencies. 
It is very scalable and runs on Linux/Unix machines, servers as well as HPC clusters. 

----------------------------------------------------------------------
Phylociraptor performs every step of a typical phylogenomic analysis.
----------------------------------------------------------------------

* It automatically downloads genomes available on NCBI and combines them with additional specified genomes provided by the user. 
* It identifies single-copy orthologs for a set of genomes (locally provided or downloaded automatically).
* It creates alignments for single-copy genes.
* It allows to trim and filter alignments in various ways.
* It calculates gene-trees for each alignment.
* It creates a species tree from gene-trees using ASTRAL.
* It creates Neighbor-Joining trees using quicktree.
* It calculates the best substitution model for single-gene alignments.
* It automatically produces all input and runs iqtree and raxml to create concatenated Maximum-Likelihood phylogenies.
* It provides extensive reports for each analysis step combined in a single HTML file which can be viewed in a Web Browser.
* Most steps are highly parallelized and it is possible to get from a list of hundreds of taxa to a phylogenomic tree in a few days. 

---------------------------------------------
Phylociraptor has the following dependencies:
---------------------------------------------

**General:**

* Singularity 3.4.1+ - `https://sylabs.io/ <https://sylabs.io/>`_
* Snakemake 6.0.2 - `https://snakemake.github.io/ <https://snakemake.github.io/>`_

dependencies:
  
	Only Singularity and Snakemake need to be installed. All other software comes in containerized form and no installation is necessary.


**Orthology inference:**

* BUSCO 3.0.2, 5.2.1  - `https://busco.ezlab.org/ <https://busco.ezlab.org/>`_

**Alignment:**

* mafft 7.464 - `https://mafft.cbrc.jp/alignment/software/ <https://mafft.cbrc.jp/alignment/software/>`_
* clustalo 1.2.4 - `http://www.clustal.org/omega/ <http://www.clustal.org/omega/>`_
* muscle 5.1 - `https://drive5.com/muscle5/ <https://drive5.com/muscle5/>`_

**Trimming:**

* trimal 1.4.1 - `http://trimal.cgenomics.org/ <http://trimal.cgenomics.org/>`_
* Aliscore/Alicut 2.31 - `https://www.zfmk.de/en/research/research-centres-and-groups/aliscore <https://www.zfmk.de/en/research/research-centres-and-groups/aliscore>`_ - `https://github.com/PatrickKueck/AliCUT <https://github.com/PatrickKueck/AliCUT>`_
* bmge 1.12 - `https://bioweb.pasteur.fr/packages/pack@BMGE@1.12/ <https://bioweb.pasteur.fr/packages/pack@BMGE@1.12/>`_

**Tree inference:**

* iqtree 2.0.7 - `http://www.iqtree.org/ <http://www.iqtree.org/>`_
* raxml-ng 1.1 - `https://github.com/amkozlov/raxml-ng <https://github.com/amkozlov/raxml-ng>`_
* astral 5.7.1 - `https://github.com/smirarab/ASTRAL <https://github.com/smirarab/ASTRAL>`_
* quicktree 2.5 - `https://github.com/khowe/quicktree <https://github.com/khowe/quicktree>`_


