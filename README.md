# phylociraptor - Rapid phylogenomic tree calculator 

This pipeline helps to create phylogenomic trees using different alignment, trimming and tree reconstruction methods. It is very scalable and runs single machines and linux clusters. It allows to automatically download user-specified genomes available on NCBI and combines them with newly sequenced genomes. It uses BUSCO to identify single copy orthologs, which are filtered, aligned, trimmed and subjected to phylogenetic inference. 

## Prerequisites
The pipline was designed in such a way that it can run desktop computers (although this is discouraged), solitary linux servers or large HPC clusters. Depending on the system setup, requirements are different: 

Local computer or solitary server:

- Linux or MacOS operating system
- globally installed singularity 3.4.1+ 
- Miniconda or Aanaconda 
- installed snakemake 5.18.3+ (eg. in an anaconda environment)

On a cluster:

- installed snakemake 5.18.3+ (eg. in an anaconda environment)
- globally installed singularity 3.4.1+
- Miniconda or Anaconda
- SGE or SLURM job scheduling system

## Getting phylociraptor

1. Clone the repository:

```
$ git clone https://github.com/reslp/smsi-phylogenomics.git
$ cd smsi-phylogenomics
```

2. Create a conda environment with snakemake:
If you don't have conda installed, first look [here](https://docs.conda.io/en/latest/miniconda.html).

```
$ conda create -c conda-forge -c bioconda -n snakemake snakemake=5.19.3
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
   set: fungi_odb9
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



**A typical run of the pipeline would look like this:**


1. Setup the pipeline:

```
$ ./phylogenomics --setup
```

2. Run BUSCO for all the genomes:

```
$ ./phylogenomics -t sge -c data/cluster_config-sauron.yaml -m busco
```

3. Create alignments and trim them:

```
$ ./phylogenomics -t sge -c data/cluster_config-sauron.yaml -m align
```

4. Reconstruct a phylogeny:

```
$ ./phylogenomics -t sge -c data/cluster_config-sauron.yaml -m tree
```

5. Create a report of the run:

```
$ ./phylogenomics --report
```


After this step, your results directory should look like this:

```
$ ls results/
alignments  assemblies  augustus_config_path  busco  busco_sequences  busco_set  busco_table  checkpoints  downloaded_genomes  filtered_alignments  phylogeny  report  report.txt  trimmed_alignments
```


## Rulegraph

<img src="https://github.com/reslp/smsi-phylogenomics/blob/master/rulegraph.png" height="500">


## Prototyping

Is done in a Docker container which has conda, snakemake and singularity installed:

```
docker run --rm -it --privileged -v $(pwd):/data reslp/smsi_ubuntu:3.4.1
```
