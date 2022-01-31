```         
                      __          __           _                  __           
               ____  / /_  __  __/ /___  _____(_)________ _____  / /_____  _____
              / __ \/ __ \/ / / / / __ \/ ___/ / ___/ __ `/ __ \/ __/ __ \/ ___/
             / /_/ / / / / /_/ / / /_/ / /__/ / /  / /_/ / /_/ / /_/ /_/ / /    
            / .___/_/ /_/\__, /_/\____/\___/_/_/   \__,_/ .___/\__/\____/_/     
           /_/          /____/                         /_/                      
```

![GitHub tag (latest by date)](https://img.shields.io/github/v/tag/reslp/phylociraptor) ![GitHub commit activity](https://img.shields.io/github/commit-activity/m/reslp/phylociraptor) ![GitHub](https://img.shields.io/github/license/reslp/phylociraptor)

# phylociraptor - Rapid phylogenomic tree calculator 

This pipeline creates phylogenomic trees for a specified set of species using different alignment, trimming and tree reconstruction methods. It is very scalable and runs on Linux/Unix machines and servers as well as HPC clusters. Phylociraptor automatically downloads genomes available on NCBI and combines them with additional specified genomes provided by the user. It uses BUSCO to identify single-copy orthologs which are filtered, aligned, trimmed and subjected to phylogenetic inference. 

*One word of caution:* phylociraptor is currently under active development. Different features may change or won't work from time to time before everything is finalized.


## Documentation

This README only provides a quick overview of phylociraptor. An extended documentation can be found here: https://phylociraptor.readthedocs.io


## Prerequisites
Phylociraptor was designed in such a way that it can run on desktop computers (although this is discouraged), solitary linux servers or large HPC clusters. Depending on the system setup, requirements are different: 

Local computer or solitary server:

- Linux or MacOS operating system
- globally installed singularity 3.4.1+
- installed snakemake 6.0.2+ (best in an anaconda environment)

or 

- Docker (this is still experimental; see information belwo)


On a cluster:

- installed snakemake 6.0.2+ (best in an anaconda environment)
- globally installed singularity 3.4.1+
- SGE or SLURM job scheduling system

## Available tools:

Orthology inference:

- BUSCO 3.0.2, 5.2.1 (https://busco.ezlab.org/)

Alignment:

- clustalo 1.2.4 (http://www.clustal.org/omega/)
- mafft 7.464 (https://mafft.cbrc.jp/alignment/software/)

Trimming:

- trimal 1.4.1 (http://trimal.cgenomics.org/)
- Aliscore/Alicut 2.31 (https://www.zfmk.de/en/research/research-centres-and-groups/aliscore; https://github.com/PatrickKueck/AliCUT)

Tree inference:

- iqtree 2.0.7 (http://www.iqtree.org/)
- raxml-ng 1.0 (https://github.com/amkozlov/raxml-ng)
- astral 5.7.1 (https://github.com/smirarab/ASTRAL)

## Getting phylociraptor

1. Clone the repository:

```
$ git clone --recursive https://github.com/reslp/phylociraptor.git
$ cd ./phylociraptor
$ ./phylociraptor

			     Welcome to
           __          __           _                  __            
    ____  / /_  __  __/ /___  _____(_)________ _____  / /_____  _____
   / __ \/ __ \/ / / / / __ \/ ___/ / ___/ __ `/ __ \/ __/ __ \/ ___/
  / /_/ / / / / /_/ / / /_/ / /__/ / /  / /_/ / /_/ / /_/ /_/ / /    
 / .___/_/ /_/\__, /_/\____/\___/_/_/   \__,_/ .___/\__/\____/_/     
/_/          /____/                         /_/                      

	  the rapid phylogenomic tree calculator, ver.0.9.0


Usage: phylociraptor <command> <arguments>

Commands:
	setup			Setup pipeline
	orthology		Infer orthologs in a set of genomes
	filter-orthology	Filter orthology results
	align			Create alignments for orthologous genes
	filter-align		Trim and filter alignments
	modeltest		Calculate gene trees and perform model testing
	mltree			Calculate Maximum-Likelihood phylogenomic trees
	speciestree		Calculate gene trees and species tree
	njtree			Calculate Neighbor-Joining tree

	report			Create a HTML report of the run
	check			Quickly check status of the run

	-v, --version 		Print version
	-h, --help		Display help

Examples:
	To see options for the setup step:
	./phylociraptor setup -h

	To run orthology inferrence for a set of genomes on a SLURM cluster:
	./phylociraptor orthology -t slurm -c data/cluster-config-SLURM.yaml

	To filter alignments overwriting the number of parsimony informative sites set in the config file:
	./phylciraptor filter-align --npars_cutoff 12

```

2. Create a conda environment with snakemake:
If you don't have conda installed, first look [here](https://docs.conda.io/en/latest/miniconda.html).

```
$ conda install -n base -c conda-forge mamba
$ mamba create -c conda-forge -c bioconda -n snakemake snakemake=6.0.2
$ conda activate snakemake
```

Additional information on how to install snakemake can be found [here](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html).

## Running phylociraptor 

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
species,Organism_Groups,Size(Mb),Chromosomes,Organelles,Plasmids,Assemblies,web_local
"Ascobolus_immersus","Eukaryota;Fungi;Ascomycetes",59.5299,0,0,0,1,web
"Ascodesmis_nigricans","Eukaryota;Fungi;Ascomycetes",27.3852,0,0,0,1,web
"Cerataphis_brasiliensis_yeast-like_symbiont","Eukaryota;Fungi;Ascomycetes",25.4655,0,0,0,1,web
"Choiromyces_venosus","Eukaryota;Fungi;Ascomycetes",126.035,0,0,0,1,web
"Endocalyx_cinctus","Eukaryota;Fungi;Ascomycetes",44.7739,0,0,0,1,web
"Margaritispora_aquatica","Eukaryota;Fungi;Ascomycetes",42.5173,0,0,0,1,web
"Morchella_conica","Eukaryota;Fungi;Ascomycetes",52.4255,0,0,0,2,web
"Neolecta_irregularis","Eukaryota;Fungi;Ascomycetes",14.1826,0,0,0,1,web
"Nilaparvata_lugens_yeast-like_symbiont","Eukaryota;Fungi;Ascomycetes",26.8096,0,0,0,1,web
"Pneumocystis_carinii","Eukaryota;Fungi;Ascomycetes",7.66146,0,0,0,1,web
"Pneumocystis_jirovecii","Eukaryota;Fungi;Ascomycetes",8.39624,0,0,0,3,web
"Protomyces_sp._C29","Eukaryota;Fungi;Ascomycetes",11.928,0,0,0,1,web
"Sclerotium_cepivorum","Eukaryota;Fungi;Ascomycetes",56.335,0,0,0,1,web
"Taxomyces_andreanae","Eukaryota;Fungi;Ascomycetes",43.1525,0,0,0,1,web
"Terfezia_boudieri","Eukaryota;Fungi;Ascomycetes",63.2346,0,0,0,1,web
"Amphirosellinia_nigrospora","Eukaryota;Fungi;Ascomycetes",48.1778,0,0,0,1,data/assemblies/assembly.fas
``` 

The basis of this file can be a CSV file directly downloaded from the [NCBI Genome Browser](https://www.ncbi.nlm.nih.gov/genome/browse#!/overview/). Just mind the changed header and additional column in the example above. The mandatory columns are the `species` and the `web_local` column. The first is the species name and the second specifies whether the genome is provided locally (in which case you should specify the path to the assembly) or not (in which case you should specify web). It is important that the species names correspond exactly to the names under which a genome is deposited at NCBI. Therefore it makes sense to use a downloaded file from the NCBI Genome Browser and add local species to them. However, you can also run the pipeline with only your own assemblies without downloading anything.



**A typical run of phylociraptor would look like this:**


1. Setup the pipeline:

```
$ ./phylociraptor setup
```

2. Identify orthologous genes for all the genomes:

```
$ ./phylociraptor orthology -t sge -c data/cluster_config-SGE.yaml
```

3. Filter orthologs using according to settings in the `config.yaml` file:

```
$ ./phylociraptor filter-orthology -t sge -c data/cluster_config-SGE.yaml
```

4. Create alignments and trim them:

```
$ ./phylociraptor align -t sge -c data/cluster_config-SGE.yaml
```

5. Filter alignments according to settings in the `config.yaml`file:

```
$ ./phylociraptor filter-align -t sge -c data/cluster_config-SGE.yaml
```

Optionally you can run extensive model testing for individual alignments. This is done using iqtree. In case you run this step, the next step will use these models. Otherwise phylociraptor will use models specified in the config file.

```
$ ./phylociraptor modeltest -t sge -c data/cluster_config-SGE.yaml
```

6. Reconstruct phylogenies:

```
$ ./phylociraptor njtree -t sge -c data/cluster_config-SGE.yaml
$ ./phylociraptor mltree -t sge -c data/cluster_config-SGE.yaml
$ ./phylociraptor speciestree -t sge -c data/cluster_config-SGE.yaml
```

7. Create a report of the run:

```
$ ./phylociraptor report
```

## Using Docker

It is also possible to run phylociraptor inside a docker container. This container has singularity and conda installed so the only requirement is that Docker is properly set up. This is still an experimental feature and we have not tested this extensively. We still recommend using phylociraptor as it is described above, especially when you work on a HPC cluster.
In case you would like to try phylociraptor in Docker you just have to subsitute the `phylociraptor` command with `phylociraptor-docker`.

```
$ git clone --recursive https://github.com/reslp/phylociraptor.git
$ cd ./phylociraptor
$ ./phylociraptor-docker
```

