#generate stats for vcf files after phasing
rule vcf_stats_initial:
	wildcard_constraints:
		vcf_name='\s+_MERGED'
	output:
		os.path.join(config.get("paths").get("base_out"),config.get("rules").get("stats").get("out_dir"),"{vcf_name}_initial.stats")
	input:
		os.path.join(BASE_OUT,config.get("rules").get("mergeReapplyVQSR").get("out_dir"),"{vcf_name}_VQSLODrefilter.vcf.gz"),
		os.path.join(BASE_OUT,config.get("rules").get("mergeReapplyVQSR").get("out_dir"),"{vcf_name}_VQSLODrefilter.vcf.gz.tbi")
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
	message: """ VCF stats after VQSR reapply """
	shell:
		"""
		{params.bcftools} stats -v {input[0]} > {output} 2> {log[1]}
		"""

rule vcf_stats_afterHWE95Clean:
	output:
		os.path.join(config.get("paths").get("base_out"),config.get("rules").get("stats").get("out_dir"),"{vcf_name}_HWE95call.stats")
	input:
		# os.path.join(BASE_OUT,config.get("rules").get("mergeReapplyVQSR").get("out_dir"),"{vcf_name}.vcf.gz"),
		# os.path.join(BASE_OUT,config.get("rules").get("mergeReapplyVQSR").get("out_dir"),"{vcf_name}.vcf.gz.tbi")
		vcf=rules.cleanMissingHwe.output[0],
		vcf_index=rules.cleanMissingHwe.output[0]		

		# rules.mergeReapplyVQSR.output[0],
		# rules.mergeReapplyVQSR.output[1]
	params:
		bcftools=config['BCFTOOLS']
	log:
		config["paths"]["log_dir"] + "/{vcf_name}_stats_HWE95call.log",
		config["paths"]["log_dir"] + "/{vcf_name}_stats_HWE95call.e"
	benchmark:
		config["paths"]["benchmark"] + "/{vcf_name}_stats_HWE95call.tsv"
	resources:
		mem_mb=5000
	envmodules:
		"bcftools/1.14"
	message: """ VCF stats after first rough clean """
	shell:
		"""
		{params.bcftools} stats -v {input[0]} > {output} 2> {log[1]}
		"""

