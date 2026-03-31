## This is the small example dataset to test phylociraptor

This small testdata set comprises 20 fungal species. We use this to demonstrate the two principal routes through a phylociraptor workflow with respect to orthology inference:
 A) identification of genes from genomes via BUSCO
 B) identification of suitable genes from protein sets via Orthofinder

Results from this example are also displayed (Figure 2) and discussed in the phylociraptor manuscript.

Run these commands from inside the phylociraptor directory in the given order to reproduce the example.

### Copy the necessary files:

Depending if you want to follow either route A (BUSCO) or B (Orthofinder) you'll need to start from different files.

**Route A** (if you chose this, do not copy relevant files for Route B below)
```
cp data/test_cases/small/config_from_genomes.yaml data/config.yaml
cp data/test_cases/small/small_from_genomes.csv data/
```

**Route B** (if you chose this, do not copy relevant files for Route A above)
```
cp data/test_cases/small/config_from_proteins.yaml data/config.yaml
cp data/test_cases/small/small_from_proteins.csv data/
```

Do inspect the config file (`data/config.yaml`) and sample file (either `data/small_from_genomes.csv` or `data/small_from_proteins.csv`) to make sense of what you are about to get into.

### setup the pipeline:

The setup stage will organise the input data and will perform slightly different tasks depending if you are going for route A or B.

At each stage in the workflow phylociraptor may automatically download required Docker containers with third party software. This usually happens at the beginning and will only happen once.

**Route A**

In this case the setup stage automatically downloads 20 genomes as listed in the sample file `data/small_from_genomes.csv`. Depending on your connection speed this will take about 5 minutes. We allow for the download to be done in three parallel batches. If you want to override this to using only a single CPU, please change to `-t local=1`. 

_IMPORTANT_: This step will download relatively large files from NCBI Genbank. Please make sure you have stable, reasonably fast connection. If this is not the case phylociraptor may abort with an error.

```
./phylociraptor setup -t local=3 # runtime ~ 10 minutes (strongly depending on connection)
```

### Modify the used busco set

We will modify the full BUSCO set that `phylociraptor setup` has downloaded in the previous stage to only 20 genes to speed up computation, using a utility function.

```
./phylociraptor util modify-busco -b fungi_odb10 -n 20 --seed 43
```

This should take only a few seconds. After the script is finished (no errors should show up).

**IMPORTANT**: You have to follow the instructions on screen and modify the config file accordingly to be able to use the modified BUSCO set!

Here's how you make the required change using a `sed` command:
```
sed -i 's/set: "fungi_odb10"/set: "fungi-seed-43-genes-20_odb10"/' data/config.yaml
```

**Route B**

In this case we will start our analyses from protein fasta files as listed in `data/small_from_proteins.csv`. The setup stage will simply organise the input data.

```
./phylociraptor setup -t local=3 # runtime ~ 1 minute
```


### Run orthology to indentify single-copy orthologs:

At this stage we'll identify orthogroups for subsequent phylogenetic inference. The actual phylociraptor call will be the same for route A and B from now on out. Note though, that what phylociraptor will actually do in detail is determined in the config file (per default `data/config.yaml`) and this file differs between route A and B.

**Route A**
For each of the 20 genomes phylociraptor will start a BUSCO run to identify the 20 randomly selected genes from the original BUSCO set in each genome, where possible.

In the config file we have specified to give each BUSCO run 2 CPU cores. Globally, the subsequent steps assume you have at least 4 CPUs on your your local machine at your disposal. This would allow phylociraptor to always run 2 BUSCO jobs in parallel and the whole stage would run for about one hour. If you have fewer resources, tune down the number accordingly, or scale it up if you have more CPUs.. 

Estimated runtime with 4 CPUs is ~ 1 hour.

**Route B**
Here phylociraptor will run Orthofinder (as specified in the config file) and establish orthogroups based on sets of proteins by all-vs-all comparison.

We give phylociraptor access to 4 CPUs for the process, which in this case will be used by Orthofinder at this stage. If you have fewer resources, tune down the number accordingly, or scale it up if you have more CPUs.. If you want to use more CPUs make sure to also change the number of threads for Orthofinder in the config file.

Estimated runtime with 4 CPUs is ~ 10 minutes.

**General comment on HPC clusters**
If you happen to have access to a HPC cluster, you can have phylociraptor automatically distribute jobs to your cluster. You'll need a cluster config file, though. We'll provide one for a SLURM cluster. You may need to adjust it slighly to the particular setup you're working with. Let us know if you need help with that. We provide an example command on a SLURM cluster below.


```
./phylociraptor orthology -t local=4

## SLURM example (after you made sure that your HPC configuration file is fine)
# ./phylociraptor orthology -t slurm -c data/cluster-config-GSC.yaml.template
```

### Filter orthology results:

This stage identifies suitable orthogroups (filter OGs based on missing data, see `data/config.yaml`) and compiles FASTA files for subsequent alignment.

```
./phylociraptor filter-orthology -t local=4 # runtime ~ 30 seconds
```

### Create alignments of single-copy orthologs:

For each orthogroup passing the filtering criteria the next stage will perform multiple sequence alignment (MSA). Phylociraptor currently has 5 aligners implemented. If indicated in the config file (`data/config.yaml`) it will perform MSA using each of the 5 - ["clustalo", "mafft", "muscle", "tcoffee", "prank"]:

```
./phylociraptor align -t local=4 # runtime ~ 1 hour

##SLURM example
#./phylociraptor align -t slurm -c data/cluster-config-GSC.yaml.template
```

### Filter alignments to remove poorly aligned regions:

Each alignment produced in the previous stage will be subjected to alignment trimming using any or all of 5 different strategies, as specified in the config file (`data/config.yaml`) - ["trimal", "aliscore", "bmge", "clipkit", "untrimmed"]:

```
./phylociraptor filter-align -t local=4 # runtime ~ 5 minutes

## SLURM example
#./phylociraptor filter-align -t slurm -c data/cluster-config-GSC.yaml.template
```

### Use filtered alignments to estimate the best substitution model and calculate gene trees

Phylociraptor will identify the best-fitting model of protein evolution and infer a ML tree for each of the trimmed multiple sequence alignments using IQTree:

```
./phylociraptor modeltest -t local=10 # runtime ~ 2 hours

## SLURM example
#./phylociraptor modeltest -t slurm -c data/cluster-config-GSC.yaml.template
```

### Calculate species trees

Phylociraptor will run ASTAL II on the sets of previously inferred ML trees to calculate a species tree. Note that according to the config file (`data/config.yaml`) we'll apply 4 different average bootstrap cutoffs to the initial trees, in this case including all, or only trees with at least 50, 60, or 70 average bootstrap support.

If you stuck to our default settings, the total number of species trees that will be calculated is 100 (5 x aligner * 5 x trimming * 4 x bootstrap cutoff).

```
./phylociraptor speciestree -t local=4 # runtime ~ 1 minute

## SLURM example
#./phylociraptor speciestree -t slurm -c data/cluster-config-GSC.yaml.template
```

### Calculate concatenated Maximum-Likelihood trees

As per the config file (`data/config.yaml`) this stage will infer ML trees based on concatenated supermatrices for each of the aligner/trimmer combination as well as multiple average bootstrap cutoffs, both using IQTree and RaXML-NG. 

Per default, the number of ML trees inferred would be 200 (5 x aligner * 5 x trimming * 4 x bootstrap cutoff * 2 x ML inference).

Note that with 10 CPUs an IQTree inference will take roughly 15 minutes. RaXML-NG takes about 2 hours per dataset (aligner/trimmer/bootstrap cutoff), so be aware that the entire process may easily take multiple days. 

```
./phylociraptor mltree -t local=10 # runtime several days

## SLURM example
#./phylociraptor mltree -t slurm -c data/cluster-config-GSC.yaml.template
```


### Create a HTML report of the run

```
./phylociraptor report
./phylociraptor report --figure
```

### A posteriori analysis of tree

Download NCBI lineage information for tree annotation:

```
./phylociraptor util get-lineage -d data/small.csv -o lineage-info.txt
```

Estimate conflicts between tree:

```
./phylociraptor util estimate-conflict -i all -o tipcov200 -s 43 -l lineage-info.txt -t 8 -b tipcoverage=200
```

Create tree plots as PDFs:

```
./phylociraptor util plot-tree --intrees $(cat tipcov200.treelist.tsv | awk '{print $2}' | tr "\n" "," | sed 's/\(.*\),/\1/') --lineagefile lineage-info.txt --level class --seed 43 --outgroup Mucor_racemosus,Glomus_cerebriforme,Smittium_simulii
```

Plot conflicts between two tree:

```
./phylociraptor util plot-conflict -i T5,T15 -q tipcov200.quartets.csv -r tipcov200.treelist.tsv -s 42 -l lineage-info.txt -e class -g Mucor_racemosus,Glomus_cerebriforme,Smittium_simulii```


