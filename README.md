# smsi-phylogenomics - A containerized pipeline to perform phylogenomics

This snakemake pipeline implements downloading a bunch of genomes for NCBI and performs phylogenomic reconstruction on them. It uses BUSCO to identify single copy orthologs and run iqtree and astral to infer phylogenies.

## Prerequisites
The pipline was designed in such a way that it can run desktop computers (although this is discouraged), solitary linux servers or large HPC clusters. As a result, needed requirements depend.

Local computer or solitary server:

- Linux or MacOS operating system
- Docker (with the possibility to run in privileged mode)
If Docker is not available:
- globally installed singularity 3.4.1+ 
- installed snakemake 5.10.0+ (eg. in an anaconda environment)

On a cluster:

- installed snakemake 5.10.0+ (eg. in an anaconda environment)
- globally installed singularity 3.4.1+
- SGE or SLURM job scheduling system

## Running the pipeline 

**You can customize the behavior of the pipeline by editing the `config.yaml` file in the `data/` folder. Two things are mandatory:**

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

**Edit your data.csv file**

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

After this step, your results directory should look like this:

```
$ ls results/
alignments  assemblies  augustus_config_path  busco  busco_sequences  busco_set  busco_table  checkpoints  downloaded_genomes  filtered_alignments  phylogeny  report.txt  trimmed_alignments
```


## Rulegraph

<img src="https://github.com/reslp/smsi-phylogenomics/blob/master/rulegraph.png" height="500">


## Prototyping

Is done in a Docker container which has conda, snakemake and singularity installed:

```
docker run --rm -it --privileged -v $(pwd):/data reslp/smsi_ubuntu:3.4.1
```
