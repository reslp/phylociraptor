rule report:
	input:
		busco = rules.busco.output,
		busco_table = rules.extract_busco_table.output.busco_table
	output:
		busco_table = "results/report/busco_table.pdf",
		busco_summary = "results/report/busco_summary.pdf"
	conda:
		"../envs/report.yml"
	singularity:
		"docker://continuumio/miniconda3:4.7.10"
	params:
		width = config["report"]["width"],
		height = config["report"]["height"],
		busco_set = config["busco"]["set"]
	shell:
		"""
		# First extract information on buscos:
		echo "species\tcomplete\tsingle_copy\tduplicated\tfragmented\tmissing\ttotal" > results/report/busco_summary.txt
		for file in $(ls results/busco/*/run_busco/short_summary_busco.txt);do  name=$(echo $file | sed 's#results/busco/##' | sed 's#/run_busco/short_summary_busco.txt##'); printf $name; cat $file | grep -P '\t\d' | awk -F "\t" '{{printf "\t"$2}}' | awk '{{print}}'; done >> results/report/busco_summary.txt
		Rscript bin/report.R {params.busco_set} {params.width} {params.height} {input.busco_table} results/report/busco_summary.txt {output.busco_table} {output.busco_summary}
		"""
