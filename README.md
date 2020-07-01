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

First specified genomes need to be downloaded:

```
$ ./phylogenomics --setup
```



## Rulegraph

<img src="https://github.com/reslp/smsi-phylogenomics/blob/master/rulegraph.png" height="500">


## Prototyping

Is done in a Docker container which has conda, snakemake and singularity installed:

```
docker run --rm -it --privileged -v $(pwd):/data reslp/smsi_ubuntu:3.4.1
```

## TODO / BUGS:

fixed - modify extract_sequences script so that it ommits empty files (buscos which are missing in all species)

fixed - iqtree rule fails if it tries to include files with less than 3 sequences. Maybe also this should be accounted for in the extract_script somehow. 
