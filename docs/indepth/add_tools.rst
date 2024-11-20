
.. role:: bash(code)
   :language: bash

======================================
Expanding phylociraptors functionality
======================================

Phylociraptor comes with a large number of tools already built in. However it is by no way a closed system. It is entirely possible to add additional tools the expand the functionality of the phylociraptor framework. 
As an example, we will walk you through how you can add an additional fictional aligner called MYALIGNER to phylociraptor. 

------------
Prerequisits
------------

Several things are necessary to be able to add your tool of choice:

1. A software container with the executable you would like to use. For example you could search for a container on  `BioContainer <https://biocontainers.pro/>`_ , `Docker Hub <https://hub.docker.com/>`_ or you could create your own container.  
2. Basic knowledge of `Snakemake <https://snakemake.readthedocs.io/en/stable/>`_.

--------------------------
Steps to add a new aligner
--------------------------

The required steps to add a new aligner include:

1. Creating a new rule file or modifying an old one for the aligner and placing it in the ``rules/aligners/`` directory.
2. Creating a new entry in the ``config.yaml.template`` file you would like to use.

Optional steps:

3. Add the container for the software to the ``data/containers.yaml`` file.
4. Expanding the cluster config file to include the new tool.


Now let's look at the steps in detail:

-----------------------------------
1. Creating/Modifying the rule file
-----------------------------------

The easiest way to do this is to look at a rule file of an already present aligner. We can the rule for the `Clustal-O <http://www.clustal.org/omega/>`_ as an example.

First we need to make a copy of the file which we will later modify. 


.. warning::

	In this example we will use the term MYALIGNER as placeholder. If you apply this to your own tool, you will have the substitute MYALIGNER with the actual name of the tool. This must be done for ALL occurences that are mentioned here.



.. code-block:: bash

	$ cp rules/aligners/clustalo.smk rules/aligners/MYALIGNER.smk 
        $ cat rules/aligners/MYALIGNER.smk
        rule clustalo:
		input:
			"results/alignments/full/clustalo.{hash}/parameters.align.clustalo.{hash}.yaml",
			sequence_file = "results/orthology/single-copy-orthologs_deduplicated." + hashes["filter-orthology"]["global"] + "/{busco}_all.fas",
		output:
			alignment = "results/alignments/full/clustalo.{hash}/{busco}_aligned.fas",
		benchmark:
			"results/statistics/benchmarks/align/clustalo_align_{busco}.{hash}.txt"
		log:
			"log/align/clustalo/clustalo_align_{busco}.{hash}.log.txt"
		singularity:
			containers["clustalo"]
		threads:
			int(config["alignment"]["threads"])
		params:
			config["alignment"]["options"]["clustalo"]
		shell:
			"""
			clustalo -i {input.sequence_file} --threads={threads} $(if [[ "{params}" != "None" ]]; then echo {params}; fi) 1> {output.alignment} 2> {log}
			"""

You can see the code which is used to run the Clustal-O aligner. We have to change all occrences of ``clustalo`` in this file to our new aligner ``MYALIGNER``. This is how it should look:


.. literalinclude:: MYALIGNER.smk
        :language: python
        :linenos:
	:emphasize-lines: 12,19


As you can see all occurences of ``clustalo`` have been replaced with MYALIGNER except for two lines which are highlighted. This first line references the file which contains a list of the containers which should be used. Here we will have to add the container for our new aligner. Theoretically it is also possible to use an executable which is locally installed and completely omit the container directive. However we discourage this to maintain reproducibility of analyses. Hence we will change the respective lines to:

.. code-block:: bash
        
        singularity:
                "docker://repository/MYALIGNER:TAG"


Mind you that this is really just a placeholder. This will need to be modified depending on where the image of your container is located.

Next we will change the command line of the program to what MYALIGNER is using. This corresponds to the second highlighted line.


.. code-block:: bash
        
        MYALIGNER -i {input.sequence_file} -t {threads} $(if [[ "{params}" != "None" ]]; then echo {params}; fi) -o {output.alignment}


This is already it for the first part. Since the command above looks a bit different lets quickly explain it:   
The ``-i`` is the hypothetical input flag of MYALIGNER followed by ``{input.sequence_file}``. This is by convention a placeholder for the unaligned input file. When phylociraptor is run, this will be replaced with anactual file path. Make sure that your write this exactly as it is giving (incl. curly brackets) otherwise phylociraptor will report an error message. After that, we specify how many threads should be used with the ``-t`` flag followed by another placeholder ``{threads}``. This information is taken fom the ``config.yaml`` file.  

The next part is this: ``$(if [[ "{params}" != "None" ]]; then echo {params}; fi)``. This section adds additional parameters to the command in case something has been specified in the ``config.yaml`` file. 

Finally with specify where the alignment should be written with the ``-o`` flag of MYALIGNER and another placeholder ``{output.alignment}``.

This is all there is to it. Just make sure to spell all the placeholder correctly and don't forget the curly brackets and everything should work.

Let us have another look at the final file:


.. literalinclude:: MYALIGNER_final.smk
        :language: python
        :linenos:


------------------------------
Modifying the config yaml file
------------------------------

The last small step is to modify the ``config.yaml`` file to include a section for additional parameters for the new aligner. This is found in the section ``align``. Let us have a look:

.. code-block:: bash

        alignment:
            method: ["clustalo","mafft","muscle"] #, "clustalo", "muscle"]
            threads: 2
            options:
                mafft: "--quiet --auto"
                clustalo: ""
                muscle: ""


Under ``options`` we have to add a new entry for our new aligner:


.. code-block:: bash

        alignment:
            method: ["clustalo","mafft","muscle", "MYALIGNER"] #, "clustalo", "muscle"]
            threads: 2
            options:
                mafft: "--quiet --auto"
                clustalo: ""
                muscle: ""
                MYALIGNER: ""


Note how we have added MYALIGNER to the list of aligners in the ``method`` section.

Between the quotes in the new option section we could add additional parameters which should be used by MYALIGNER. This information is transferred automatically to all alignment jobs for MYALIGNER when phylociraptor is run.

This is already it. You have now successfully added an additional aligner to phylociraptor. 

