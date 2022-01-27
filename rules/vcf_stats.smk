#generate stats for vcf files after phasing
rule vcf_stats_initial:
	output:
		os.path.join(config.get("paths").get("base_out"),config.get("rules").get("stats").get("out_dir"),"{vcf_name}_initial.stats")
	input:
		os.path.join(BASE_OUT,config.get("rules").get("mergeReapplyVQSR").get("out_dir"),"{vcf_name}_VQSLODrefilter.vcf.gz"),
		os.path.join(BASE_OUT,config.get("rules").get("mergeReapplyVQSR").get("out_dir"),"{vcf_name}_VQSLODrefilter.vcf.gz.tbi"),
		rules.getSamples.output[0]
	params:
		bcftools=config['BCFTOOLS']
	log:
		config["paths"]["log_dir"] + "/{vcf_name}_stats_initial.log",
		config["paths"]["log_dir"] + "/{vcf_name}_stats_initial.e"
	benchmark:
		config["paths"]["benchmark"] + "/{vcf_name}_stats_initial.tsv"
	resources:
		mem_mb=10000
	envmodules:
		"bcftools/1.14"
	message: """ VCF stats after VQSR reapply """
	shell:
		"""
		{params.bcftools} stats -S {input[2]} -v {input[0]} > {output} 2> {log[1]}
		"""

rule vcf_stats_afterHWE95Clean:
	output:
		# os.path.join(config.get("paths").get("base_out"),config.get("rules").get("stats").get("out_dir"),"{vcf_name}_{chr}_HWE95call.stats")
		os.path.join(config.get("paths").get("base_out"),config.get("rules").get("stats").get("out_dir"),"{vcf_name}_HWE95call.stats")
	input:
		# os.path.join(BASE_OUT,config.get("rules").get("mergeReapplyVQSR").get("out_dir"),"{vcf_name}.vcf.gz"),
		# os.path.join(BASE_OUT,config.get("rules").get("mergeReapplyVQSR").get("out_dir"),"{vcf_name}.vcf.gz.tbi")
		vcf=rules.cleanMissingHwe.output[0],
		vcf_index=rules.cleanMissingHwe.output[0],		
		samples=rules.getSamples.output[0]
	params:
		bcftools=config['BCFTOOLS']
	log:
		config["paths"]["log_dir"] + "/{vcf_name}-stats_HWE95call.log",
		config["paths"]["log_dir"] + "/{vcf_name}-stats_HWE95call.e"
	benchmark:
		config["paths"]["benchmark"] + "/{vcf_name}_stats_HWE95call.tsv"
	resources:
		mem_mb=10000
	envmodules:
		"bcftools/1.14"
	message: """ VCF stats after first rough clean """
	shell:
		"""
		{params.bcftools} stats -S {input.samples} -v {input.vcf} > {output[0]} 2> {log[1]}
		"""

