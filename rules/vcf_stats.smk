#generate stats for vcf files after phasing
rule vcf_stats_initial:
	output:
		os.path.join(config.get("paths").get("base_out"),config.get("rules").get("stats").get("out_dir"),"{vcf_name}_initial.stats")
	input:
		os.path.join(BASE_OUT,config.get("rules").get("mergeReapplyVQSR").get("out_dir"),"{vcf_name}.vcf.gz"),
		os.path.join(BASE_OUT,config.get("rules").get("mergeReapplyVQSR").get("out_dir"),"{vcf_name}.vcf.gz.tbi")
		# rules.mergeReapplyVQSR.output[0],
		# rules.mergeReapplyVQSR.output[1]
	params:
		bcftools=config['BCFTOOLS']
	log:
		config["paths"]["log_dir"] + "/{vcf_name}_stats_initial.log",
		config["paths"]["log_dir"] + "/{vcf_name}_stats_initial.e"
	benchmark:
		config["paths"]["benchmark"] + "/{vcf_name}_stats_initial.tsv"
	resources:
		mem_mb=5000
	envmodules:
		"bcftools/1.14"
	message: """ VCF stats after phasing """
	shell:
		"""
		{params.bcftools} stats -v {input[0]} > {output} 2> {log[1]}
		"""