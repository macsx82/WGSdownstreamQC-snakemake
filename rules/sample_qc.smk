#module containing sample based QC for WGS data

#extract singleton count for each sample
rule singletons:
	# wildcard_constraints:
	# 	vcf_name=''
	output:
		expand(os.path.join(BASE_OUT,config.get("rules").get("singletons").get("out_dir"), "{vcf_name}_singletons.{ext}"), ext=["singletons", "log"])
	input:
		vcf=os.path.join(BASE_OUT,config.get("rules").get("mergeReapplyVQSR").get("out_dir"),"{vcf_name}.vcf.gz"),
		vcf_index=os.path.join(BASE_OUT,config.get("rules").get("mergeReapplyVQSR").get("out_dir"),"{vcf_name}.vcf.gz.tbi")
	params:
		bcftools=config['BCFTOOLS'],
		vcftools=config['VCFTOOLS'],
		tmp=config.get("paths").get("tmp"),
		out_prefix=os.path.join(BASE_OUT,config.get("rules").get("singletons").get("out_dir"), "{vcf_name}_singletons")
	log:
		config["paths"]["log_dir"] + "/{vcf_name}-singletons.log",
		config["paths"]["log_dir"] + "/{vcf_name}-singletons.e"
	threads: 1
	resources:
		mem_mb=5000
	benchmark:
		config["paths"]["benchmark"] + "/{vcf_name}_singletons.tsv"
	envmodules:
		"vcftools/0.1.16"
	shell:
		"""
		{params.vcftools} --gzvcf {input} --singletons --out {params.out_prefix}
		"""

#we may need data from varcall pipeline
#or we can use Read Depth by Sample
# rule sampleDP:
# 	output:
# 	input:
# 	params:
# 		bcftools=config['BCFTOOLS'],
# 		vcftools=config['VCFTOOLS'],
# 		tmp=os.path.join(BASE_OUT,config.get("paths").get("tmp")),
# 		out_prefix=
# 	log:
# 		config["paths"]["log_dir"] + "/{vcf_name}-sampleDP.log",
# 		config["paths"]["log_dir"] + "/{vcf_name}-sampleDP.e"
# 	threads: 1
# 	resources:
# 		mem_mb=5000
# 	benchmark:
# 		config["paths"]["benchmark"] + "/{vcf_name}_sampleDP.tsv"
# 	envmodules:
# 		"vcftools/0.1.16"
# 	shell:
# 		"""
# 		{params.vcftools} --gzvcf {input} --depth --out {params.out_prefix}
# 		"""

# rule coverage:
# 	output:
# 	input:
# 	params:
# 	log:
# 	threads:
# 	resources:
# 	benchmark:
# 	shell:

# #het rate ruel: vcftool or custom script
# rule SampleHetRate:
# 	output:
# 	input:
# 	params:
# 		bcftools=config['BCFTOOLS'],
# 		vcftools=config['VCFTOOLS'],
# 		tmp=os.path.join(BASE_OUT,config.get("paths").get("tmp")),
# 		out_prefix=
# 	log:
# 		config["paths"]["log_dir"] + "/{vcf_name}-hetRate.log",
# 		config["paths"]["log_dir"] + "/{vcf_name}-hetRate.e"
# 	threads: 1
# 	resources:
# 		mem_mb=5000
# 	benchmark:
# 		config["paths"]["benchmark"] + "/{vcf_name}_hetRate.tsv"
# 	envmodules:
# 		"vcftools/0.1.16"
# 	shell:
# 		"""
# 		{params.vcftools} --gzvcf {input} --het --out {params.out_prefix}
# 		"""

# rule SampleGetHetRateOut:
# 	output:
# 	input:
# 		rules.SampleHetRate.output[0]
# 	params:
# 	log:
# 	threads:
# 	resources:
# 	benchmark:
# 	run:
# 		get_het_sample_outliers(input[0], output[0])


# #Missing rate rule: plink or vcftools or bcftools command
# rule SampleMissingRate:
# 	output:
# 	input:
# 	params:
# 	log:
# 	threads:
# 	resources:
# 	benchmark:
# 	shell:


# rule PCA:
# 	output:
# 	input:
# 	params:
# 	log:
# 	threads:
# 	resources:
# 	benchmark:
# 	shell:



