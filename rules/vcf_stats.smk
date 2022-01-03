#generate stats for vcf files after phasing
rule vcf_stats_final:
	output:
		os.path.join(config.get("files_path").get("base_joint_call_path"),config.get("rules").get("phase").get("out_dir"),"all.{interval_name}.stats")
	input:
		rules.phase.output[0],rules.phase.output[1]
	params:
		bcftools=config['BCFTOOLS']
	log:
		config["files_path"]["log_dir"] + "/{interval_name}-vcf_stats_phase.log",
		config["files_path"]["log_dir"] + "/{interval_name}-vcf_stats_phase.e"
	benchmark:
		config["files_path"]["benchmark"] + "/{interval_name}_vcf_stats_phase.tsv"
	resources:
		mem_mb=5000
	envmodules:
    	"bcftools/1.14"
	message: """ VCF stats after phasing """
	shell:
		"""
		{params.bcftools} stats -v {input[0]} > {output} 2> {log[1]}
		"""