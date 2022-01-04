#generate stats for vcf files after phasing
rule vcf_stats_initial:
	output:
		os.path.join(config.get("patha").get("base_out"),config.get("rules").get("stats").get("out_dir"),"all.{vcf_name}.stats")
	input:
		rules.mergeReapplyVQSR.output[0],
		rules.mergeReapplyVQSR.output[1]
	params:
		bcftools=config['BCFTOOLS']
	log:
		config["files_path"]["log_dir"] + "/vcf_stats_initial.log",
		config["files_path"]["log_dir"] + "/vcf_stats_initial.e"
	benchmark:
		config["files_path"]["benchmark"] + "/vcf_stats_initial.tsv"
	resources:
		mem_mb=5000
	envmodules:
		"bcftools/1.14"
	message: """ VCF stats after phasing """
	shell:
		"""
		{params.bcftools} stats -v {input[0]} > {output} 2> {log[1]}
		"""