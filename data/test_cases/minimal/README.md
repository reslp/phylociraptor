## This is a minimal example to quickly test if phylociraptor works

Run these commands from inside the phylociraptor directory in the given order to reproduce the minimal example. 

This minimal test dataset consists of six fungal genomes.

Results generated from this run:

- 6 fungal genome assemblies downloaded
- 10 single-copy orthologs identified in each genome (using BUSCO5)
- 18 multiple sequence alignments using mafft and clustal omega. (Could have been 20 but two orthologs are missing in too many genomes)
- 36 trimmed and filtered alignments using trimal and aliscore/alicut (could have been 40, but 4 were eliminated by our filters)
- 36 maximum-likelihood gene trees calculated with iqtree
- 12 concatenated maximum-likelihood phylogenies using 3 different bootstrap-cutoff values
- 12 species trees using 3 different bootstrap-cutoff values
- 1 comprehensive report of our analyses in HTML format

The below commands make use of a single computational thread only. If you have more resources available on your machine you can adjust the part `-t local=1` to however many CPU threads you actually have available, say four threads: `-t local=4`. This will result in parallelization of most of the processes and speed up your analyses.

For some more details about what's happening below, have a look at our tutorial [here](https://phylociraptor.readthedocs.io/en/latest/introduction/minimal.html).

### setup the pipeline:

This runs using only one thread (`-t local=1) on the same machine. Adjust the number of threads to be used if you want (see above). It should take about 5 minutes.

```
./phylociraptor setup --config-file data/test_cases/minimal/config.yaml -t local=1 --verbose
```

### Subsample the BUSCO set

For this minimal example we need to modify the busco set to include only a small number (10) of genes. We do this so the remaining steps complete quickly.
To reduce the busco set to only 10 genes run this command in the phylocirpator directory:

```
./phylociraptor util modify-busco -b fungi_odb10 -n 10 --seed 42
```

This should take only a few seconds. 

To resume with the new, reduced BUSCO set make a copy of the settings file and modify the BUSCO set specified. You can also modify the file manually with a text editor, but here we do it programmatically. Note that this solution works in this case, given the files and BUSCO set as used in the tutorial. Adjust for your own set if needed.

```
cp data/test_cases/minimal/config.yaml data/config.yaml 
sed -i 's/fungi_odb10/fungi-seed-42-genes-10_odb10/' data/config.yaml
```

### Run orthology to indentify single-copy orthologs:

This should take about 15 minutes to complete with a single thread.

```
./phylociraptor orthology --verbose -t local=1
```

### Filter orthology results:

This step should take about one minute to complete

```
./phylociraptor filter-orthology --verbose -t local=1
```

### Create alignments of single-copy orthologs:

With a single thread (`local=1`) this takes about 3 minutes.

```
./phylociraptor align --verbose -t local=1
```

### Filter alignments to remove poorly aligned regions:

This takes about 4 minutes to complete.

```
./phylociraptor filter-align --verbose -t local=1
```

### Use filtered alignments to estimate the best substitution model and calculate gene trees

This takes about 15 minutes to run using 1 thread.

```
./phylociraptor modeltest --verbose -t local=1
```

### Calculate concatenated Maximum-Likelihood trees

This step takes about 15 minutes when using 1 thread.

```
./phylociraptor mltree --verbose -t local=1
```

### Calculate species trees

This step takes 2 minutes to complete using 1 thread:

```
./phylociraptor speciestree --verbose -t local=1
```

### Create a report of the run

This takes about 1 minute:

```
./phylociraptor report
```


