## This is the small example dataset to test phylociraptor mentioned in the phylociraptor manuscript

Run these commands from inside the phylociraptor directory in the given order to reproduce the example.

### Copy the necessary files:

```
cp data/test_cases/small/config.yaml data/
cp data/test_cases/small/small.csv data/
```

### setup the pipeline:

This runs on a single thread (serial=1) and the same machine. It should take about 5 minutes.

```
./phylociraptor setup -t serial=1
```

### Modify the used busco set

We will modify the downloaded BUSCO set to only 20 genes to speed up computation.

```
./phylociraptor util modify-busco -b fungi_odb9 -n 20 --seed 43
```

This should take only a few seconds. After the script is finished (no errors should show up).

**IMPORTANT**: You have to follow the instructions on screen and modify the config file accordingly to be able to use the modified BUSCO set!


### Run orthology to indentify single-copy orthologs:

The subsequent steps use a SLURM cluster. You may have to modify the correpsonding cluster config file accordingly.

```
./phylociraptor orthology -t slurm -c data/cluster-config-GSC.yaml.template
```

### Filter orthology results:

```
./phylociraptor filter-orthology -t local -c data/cluster-config-GSC.yaml.template
```

### Create alignments of single-copy orthologs:

```
./phylociraptor align -t slurm -c data/cluster-config-GSC.yaml.template
```

### Filter alignments to remove poorly aligned regions:

```
./phylociraptor filter-align -t slurm -c data/cluster-config-GSC.yaml.template
```

### Use filtered alignments to estimate the best substitution model and calculate gene trees


```
./phylociraptor modeltest -t slurm -c data/cluster-config-GSC.yaml.template
```

### Calculate concatenated Maximum-Likelihood trees


```
./phylociraptor mltree -t slurm -c data/cluster-config-GSC.yaml.template
```

### Calculate species trees

```
./phylociraptor speciestree -t slurm -c data/cluster-config-GSC.yaml.template
```

### Create a reports of the run

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


