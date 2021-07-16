.. role:: bash(code)
    :language: bash

==================
Quickstart
==================

-----------------------
Set up your analysis
-----------------------

**To customize the behavior of the pipeline to fit your needs you can edit the `config.yaml` file in the `data/` folder. Two things are mandatory:**

1. You need to enter the correct name for the data.csv containing the species which should be included in the tree:

.. code-block:: bash

	species: data/data.csv


2. Another thing you need is to specify the correct BUSCO set and augustus species:

.. code-block:: bash

	busco:
   		set: fungi_odb9
   		ausgustus_species: anidulans


**You will also need to provide a list of genomes which should be used in your analysis. To do this, edit your data.csv file**

The data.csv file should look something like this:

.. code-block:: bash

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

The basis of this file can be a CSV file directly downloaded from the `NCBI Genome Browser <https://www.ncbi.nlm.nih.gov/genome/browse#!/overview/>`_ . Just mind the changed header and additional column in the example above. The mandatory columns are the :bash:`species` and the :bash:`web_local` column. The first is the species name and the second specifies whether the genome is provided locally (in which case you should specify the path to the assembly) or not (in which case you should specify web). It is important that the species names correspond exactly to the names under which a genome is deposited at NCBI. Therefore it makes sense to use a downloaded file from the NCBI Genome Browser and add local species to them. However, you can also run the pipeline with only your own assemblies without downloading anything.

---------------------
Running phylociraptor
---------------------

**A typical run of phylociraptor would look like this:**

.. note::

	These example commands here assume phylociraptor is run on an SGE cluster.

**1. Setup the pipeline:**

.. code-block:: bash
	
	$ ./phylociraptor setup -t sge


During this step, phylociraptor will download and organize all necessary input data. This includes downloading genome assemblies for NCBI (if there are any), downloading the BUSCO set
specified in the :bash:`config.yaml` file. 

**2. Infer and filter orthologous genes for all the genomes:**

.. code-block:: bash

	$ ./phylociraptor orthology --cluster sge --cluster-config data/cluster-config-SGE.yaml
	$ ./phylociraptor filter-orthology --cluster sge --cluster-config data/cluster-config-SGE.yaml

.. note::

	If you don't specify a cluster configuration file, phylociraptor will try to use the default files, which (most likely) will not work on your cluster.

**3. Create alignments and trim them:**

.. code-block:: bash

	$ ./phylociraptor align -t sge -c data/cluster-config-SGE.yaml

Optionally you can run extensive model testing for individual alignments. This is done using iqtree. In case you run this step, the next step will use these models. Otherwise phylociraptor will use models specified in the config file.

.. code-block:: bash

	$ ./phylociraptor model -t sge -c data/cluster-config-SGE.yaml


**4. Reconstruct a phylogeny:**

.. code-block:: bash

	$ ./phylociraptor tree -t sge -c data/cluster-config-SGE.yaml
	$ ./phylociraptor njtree -t sge -c data/cluster-config-SGE.yaml
	$ ./phylociraptor speciestree -t sge -c data/cluster-config-SGE.yaml


**5. Create a report of the run:**

.. code-block:: bash

	$ ./phylociraptor report


After this step, your results directory should look like this:

.. code-block:: bash

	$ ls results/
	alignments  assemblies  augustus_config_path  busco  busco_sequences  busco_set  busco_table  checkpoints  downloaded_genomes  filtered_alignments  phylogeny  report  report.txt  trimmed_alignments


