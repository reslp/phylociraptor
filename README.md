
<p float="left">
  <img src="/docs/images/logo.png" width="200" />
  <img src="/docs/images/phylociraptor.png" width="600" /> 
</p>

![GitHub tag (latest by date)](https://img.shields.io/github/v/tag/reslp/phylociraptor) ![GitHub commit activity](https://img.shields.io/github/commit-activity/m/reslp/phylociraptor) ![GitHub](https://img.shields.io/github/license/reslp/phylociraptor) [![DOI](https://zenodo.org/badge/258200660.svg)](https://zenodo.org/badge/latestdoi/258200660)



# phylociraptor - Rapid phylogenomic tree calculator 

Phylociraptor is a computational framework calculating phylogenomic trees for a specified set of species using different alignment, trimming and tree reconstruction methods. It is very scalable and runs on Linux/Unix machines and servers as well as HPC clusters. Phylociraptor automatically downloads genomes available on NCBI and combines them with additional specified genomes provided by the user. It uses BUSCO or OrthoFinder to identify single-copy orthologs which are filtered, aligned, trimmed and subjected to phylogenetic inference using state of the art software. 

**Want to learn more about**  ...

... [the phylociraptor workflow](#The-phylociraptor-workflow)  
... [prerequisistes and dependencies](#Prerequisites)  
... [installation](#Installing-phylociraptor)  
... [the incorporated software](#Available-tools)  
... [the many customization possibilities](#Customization)  
... [a-posteri analyses of phylogenomic trees](#A-posteriori-analyses)  
... [how to use phylociraptor with docker](#Using-Docker)  
... [how to cite phylociraptor](#Citation)  
... [all the details in the documentation](https://phylociraptor.readthedocs.io/en/latest/)  


## The phylociraptor workflow

**A typical run of phylociraptor looks like this:**

<p float="left">
  <img src="/docs/images/overview.png" width="800" />
</p>

For more information including a short tutorial please refer to the [documentation](https://phylociraptor.readthedocs.io).

### 1. Setup the pipeline:

```
$ ./phylociraptor setup -t local
```

### 2. Identify orthologous genes for all the genomes:

```
$ ./phylociraptor orthology -t local
```

### 3. Filter orthologs using according to settings in the `config.yaml` file:

```
$ ./phylociraptor filter-orthology -t local
```

### 4. Create alignments and trim them:

```
$ ./phylociraptor align -t local
```

### 5. Filter alignments according to settings in the `config.yaml` file:

```
$ ./phylociraptor filter-align -t local
```

### 6. Run modeltesting for individual alignments:

```
$ ./phylociraptor modeltest -t local
```

### 7. Reconstruct phylogenies:

```
$ ./phylociraptor njtree -t local
$ ./phylociraptor mltree -t local
$ ./phylociraptor speciestree -t local
```

### 8. Create a report of the run:

```
$ ./phylociraptor report
```

## Prerequisites
Phylociraptor was designed in such a way that it can run on desktop computers (although this is discouraged), solitary linux servers or large HPC clusters. Depending on the system setup, requirements are different: 

Local computer or solitary server:

- Linux or MacOS operating system
- globally installed singularity 3.4.1+
- installed snakemake 6.0.2+ (best in an anaconda environment)

or 

- Docker (this is still experimental; see information below)

On a HPC cluster:

- installed snakemake 6.0.2+ (best in an anaconda environment)
- globally installed singularity 3.4.1+
- SGE, SLURM or TORQUE job scheduling system

## Installing phylociraptor

1. Create a conda environment with snakemake:
If you don't have conda installed, first look [here](https://docs.conda.io/en/latest/miniconda.html).

```
$ conda install -n base -c conda-forge mamba
$ mamba create -c conda-forge -c bioconda -n snakemake snakemake=6.0.2
$ conda activate snakemake
```

Additional information on how to install snakemake can be found [here](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html).

Now we can proceed with installing phylociraptor.

2. Clone the repository:

```
(snakemake) $ git clone --recursive https://github.com/reslp/phylociraptor.git
(snakemake) $ cd phylociraptor
(snakemake) $ ./phylociraptor

			     Welcome to
           __          __           _                  __            
    ____  / /_  __  __/ /___  _____(_)________ _____  / /_____  _____
   / __ \/ __ \/ / / / / __ \/ ___/ / ___/ __ `/ __ \/ __/ __ \/ ___/
  / /_/ / / / / /_/ / / /_/ / /__/ / /  / /_/ / /_/ / /_/ /_/ / /    
 / .___/_/ /_/\__, /_/\____/\___/_/_/   \__,_/ .___/\__/\____/_/     
/_/          /____/                         /_/                      

	  the rapid phylogenomic tree calculator, ver.1.0.0


Usage: phylociraptor <command> <arguments>

Commands:
	setup			Setup pipeline
	orthology		Infer orthologs in a set of genomes
	filter-orthology	Filter orthology results
	align			Create alignments for orthologous genes
	filter-align		Trim and filter alignments
	modeltest		Calculate gene trees and perform model testing
	mltree			Calculate Maximum-Likelihood phylogenomic trees
	speciestree		Calculate species trees
	njtree			Calculate Neighbor-Joining trees
        bitree                  Calculate Bayesian-inference phylogenomic trees

	report			Create a HTML report of the run
	check			Quickly check status of the run
	util			Utilities for a posteriori analyses of trees

	-v, --version 		Print version
	-h, --help		Display help

Examples:
	To see options for the setup step:
	./phylociraptor setup -h

	To run orthology inferrence for a set of genomes on a SLURM cluster:
	./phylociraptor orthology -t slurm -c data/cluster-config-SLURM.yaml

```

## Available tools

Orthology inference:

- BUSCO 3.0.2, 5.2.1 (https://busco.ezlab.org/)
- OrthoFinder 2.5.5 (https://github.com/davidemms/OrthoFinder)

Alignment:

- Clustal-Omega 1.2.4 (http://www.clustal.org/omega/)
- MAFFT 7.464 (https://mafft.cbrc.jp/alignment/software/)
- MUSCLE 5.1 (https://drive5.com/muscle5/)
- T-Coffee 13.46.0.919e8c6b (https://github.com/cbcrg/tcoffee)
- PRANK v150803 (https://ariloytynoja.github.io/prank-msa/)

Trimming:

- trimAl 1.4.1 (http://trimal.cgenomics.org/)
- Aliscore/Alicut 2.31 (https://www.zfmk.de/en/research/research-centres-and-groups/aliscore; https://github.com/PatrickKueck/AliCUT)
- BMGE 1.12 (https://bioweb.pasteur.fr/packages/pack@BMGE@1.12/)
- ClipKit 2.3.0 (https://github.com/JLSteenwyk/ClipKIT)

Tree inference:

- IQ-Tree 2.0.7 (http://www.iqtree.org/)
- RAxML-NG 1.1 (https://github.com/amkozlov/raxml-ng)
- ASTRAL 5.7.1 (https://github.com/smirarab/ASTRAL)
- Quicktree 2.5 (https://github.com/khowe/quicktree)
- Phylobayes-MPI 1.9 (https://github.com/bayesiancook/pbmpi)


## Customization

For a comprehensive overview of phylociraptors customization option please refer to our [documentation](https://phylociraptor.readthedocs.io).  
**To customize the behavior of the pipeline to fit your needs you can edit the `config.yaml` file in the `data/` folder. Two things are mandatory:**

1. You need to enter the correct name for the data.csv containing the species which should be included in the tree:

```
species: data/data.csv
```

2. Another thing you need is to specify the correct BUSCO set and augustus species:

```
busco:
   set: "fungi_odb9"
   ausgustus_species: anidulans
```

**You will also need to provide a list of genomes which should be used in your analysis. To do this, edit your data.csv file**

The data.csv file should look something like this:

```
$ cat data.csv
species,web_local,mode
Salpingoeca rosetta,web=GCA_000188695.1,
Coccidioides posadasii,web,
Sclerophora sanguinea,data/assemblies/Sclerophora_sanguinea_Sclsan1_AssemblyScaffolds_Repeatmasked.fasta,
Capsaspora owczarzaki,web=GCA_000151315.2,
Dictyostelium lacteum,web=GCA_001606155.1,
Paraphelidium tribonemae,data/assemblies/EP00158_Paraphelidium_tribonemae.fasta,protein
Synchytrium microbalum,web=GCA_006535985.1,
Nuclearia sp,data/assemblies/Nuclearia_sp_trinity_cdhit-0.95.fasta,transcriptome
Aspergillus nidulans,data/assemblies/Asp_nid_full_proteome.fasta,protein
Stereomyxa ramosa,data/assemblies/Stereomyxa_ramosa_trinity_cdhit-0.95.fasta,transcriptome
``` 

The basis of this file can be a CSV file directly downloaded from the [NCBI Genome Browser](https://www.ncbi.nlm.nih.gov/genome/browse#!/overview/). Just mind the changed header and additional column in the example above. The mandatory columns are the `species` and the `web_local` column. The first is the species name and the second specifies whether the genome is provided locally (in which case you should specify the path to the assembly) or not (in which case you should specify web).  
It is important that the species names correspond exactly to the names under which a genome is deposited at NCBI. Therefore it makes sense to use a downloaded file from the NCBI Genome Browser and add local species to them. However, you can also run the pipeline with only your own assemblies without downloading anything. It is also possible to use transcriptome assemblies or sets of proteins. Please refer to the documentation for more information.


## A posteriori analyses

Phylociraptor provides several utilities to investigate tree similarity and plot trees. Please refer to the [documentation](https://phylociraptor.readthedocs.io) for additional details.

## Using Docker

It is also possible to run phylociraptor inside a docker container. This container has singularity and conda installed so the only requirement is that Docker is properly set up. This is still an experimental feature and we have not tested this extensively. We still recommend using phylociraptor as it is described above, especially when you work on a HPC cluster.
In case you would like to try phylociraptor in Docker you just have to subsitute the `phylociraptor` command with `phylociraptor-docker`.

```
$ git clone --recursive https://github.com/reslp/phylociraptor.git
$ cd ./phylociraptor
$ ./phylociraptor-docker
```

## Citation

If you use phylociraptor please cite our preprint, which is available [here](https://www.biorxiv.org/content/10.1101/2023.09.22.558970v1.article-info). 

**Resl Philipp & Hahn Christoph** (2023) *Phylociraptor - A unified computational framework for reproducible phylogenomic inference* (Preprint) bioRxiv doi:10.1101/2023.09.22.558970
