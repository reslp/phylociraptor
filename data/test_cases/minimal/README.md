## This is a minimal example to quickly test if phylociraptor works

Run these commands from inside the phylociraptor directory in the given order to reproduce the minimal example.

This minimal test dataset consists of six fungal genomes.

Results generated from this run:

- 6 downloaded fungal genome assemblies
- 10 identified single-copy orthologs in each genome
- 20 multiple sequence alignments using mafft and clustal-o
- 38 trimmed alignments using trimal and aliscore/alicut (2 will be excluded during trimming)
- 38 maximum-likelihood gene trees calculated with iqtree
- 12 concatenated maximum-likelihood phylogenies using 3 different bootstrap-cutoff values
- 12 species trees using 3 different bootstrap-cutoff values

### Copy the necessary files:

```
cp data/test_cases/minimal/config.yaml data/
cp data/test_cases/minimal/minimal.csv data/
cp data/test_cases/minimal/modify_busco.sh .
```

### setup the pipeline:

This runs on a single thread (serial=1) and the same machine. It should take about 5 minutes.

```
./phylociraptor setup --verbose -t serial=1
```

### Modify the used busco set

For this minimal example we need to modify the busco set to include only a small number (10) of genes. We do this so the remaining steps complete quickly.
To reduce the busco set to only 10 genes run this command in the phylocirpator directory:

```
./modify_busco.sh
```

This should take only a few seconds. After the script is finished (no errors should show up), we can continue with the next step:

### Run orthology to indentify single-copy orthologs:

This step runs using four threads (serial=4). It will run locally and should take about 5 minutes to complete

```
./phylociraptor orthology --verbose -t serial=4
```

### Filter orthology results:

This step should take about one minute to complete

```
./phylociraptor filter-orthology --verbose -t serial=1
```

### Create alignments of single-copy orthologs:

With a single thread (serial=1) this takes about 3 minutes.

```
./phylociraptor align --verbose -t serial=1
```

### Filter alignments to remove poorly aligned regions:

This takes about 4 minutes to complete.

```
./phylociraptor filter-align --verbose -t serial=1
```

### Use filtered alignments to estimate the best substitution model and calculate gene trees

This takes about 5 minutes to run using 4 threads.

```
./phylociraptor modeltest --verbose -t serial=4
```

### Calculate concatenated Maximum-Likelihood trees

This step takes about 30 minutes when using 4 threads.

```
./phylociraptor mltree --verbose -t serial=4
```

### Calculate species trees

This step takes xx minutes to complete using 4 threads:

```
./phylociraptor speciestree --verbose -t serial=4
```

### Create a report of the run

This takes about 1 minute:

```
./phylociraptor report
```


