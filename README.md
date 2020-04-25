# SMSI Phylogenomics

Prototyping is done in:

```
docker run --rm -it --privileged -v $(pwd):/data reslp/smsi_ubuntu:3.4.1
```

TODO / BUGS:

fixed - modify extract_sequences script so that it ommits empty files (buscos which are missing in all species)

fixed - iqtree rule fails if it tries to include files with less than 3 sequences. Maybe also this should be accounted for in the extract_script somehow. 