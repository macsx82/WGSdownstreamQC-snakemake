#module containing sample based QC for WGS data

#extract singleton count for each sample
rule singletons:
	output:
		# expand(os.path.join(BASE_OUT,config.get("rules").get("singletons").get("out_dir"), "{{vcf_name}}_singletons.{ext}"), ext=["singletons", "log"])
		os.path.join(BASE_OUT,config.get("rules").get("singletons").get("out_dir"), "{vcf_name}_singletons.singletons")
	input:
		# vcf=os.path.join(BASE_OUT,config.get("rules").get("mergeReapplyVQSR").get("out_dir"),"{vcf_name}.vcf.gz"),
		# vcf_index=os.path.join(BASE_OUT,config.get("rules").get("mergeReapplyVQSR").get("out_dir"),"{vcf_name}.vcf.gz.tbi")
		vcf=rules.cleanMissingHwe.output[0],
		vcf_index=rules.cleanMissingHwe.output[0]
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
		{params.vcftools} --gzvcf {input.vcf} --singletons --out {params.out_prefix} 1> {log[0]} 2> {log[1]}
		"""

#het rate rule: first get the data with vcftools
rule SampleHetRate:
	output:
		expand(os.path.join(BASE_OUT,config.get("rules").get("SampleHetRate").get("out_dir"), "{{vcf_name}}_het.{ext}"), ext=["het"])
	input:
		# vcf=os.path.join(BASE_OUT,config.get("rules").get("mergeReapplyVQSR").get("out_dir"),"{vcf_name}.vcf.gz"),
		# vcf_index=os.path.join(BASE_OUT,config.get("rules").get("mergeReapplyVQSR").get("out_dir"),"{vcf_name}.vcf.gz.tbi")
		vcf=rules.cleanMissingHwe.output[0],
		vcf_index=rules.cleanMissingHwe.output[0]
	params:
		bcftools=config['BCFTOOLS'],
		vcftools=config['VCFTOOLS'],
		tmp=os.path.join(BASE_OUT,config.get("paths").get("tmp")),
		out_prefix=os.path.join(BASE_OUT,config.get("rules").get("SampleHetRate").get("out_dir"), "{vcf_name}_het")
	log:
		config["paths"]["log_dir"] + "/{vcf_name}-het.log",
		config["paths"]["log_dir"] + "/{vcf_name}-het.e"
	threads: 1
	resources:
		mem_mb=5000
	benchmark:
		config["paths"]["benchmark"] + "/{vcf_name}_het.tsv"
	envmodules:
		"vcftools/0.1.16"
	shell:
		"""
		{params.vcftools} --gzvcf {input.vcf} --het --out {params.out_prefix} 1> {log[0]} 2> {log[1]}
		"""

#het rate rule: get vcftools result and extract the het rate for plotting
rule SampleGetHetRateOut:
	output:
		os.path.join(BASE_OUT,config.get("rules").get("SampleHetRate").get("out_dir"), "{vcf_name}_hetRate.txt")
	input:
		rules.SampleHetRate.output[0]
	params:
		vcftools=config['VCFTOOLS']
	log:
		config["paths"]["log_dir"] + "/{vcf_name}-hetRate.log",
		config["paths"]["log_dir"] + "/{vcf_name}-hetRate.e"
	threads: 1
	resources:
		mem_mb=5000
	benchmark:
		config["paths"]["benchmark"] + "/{vcf_name}_hetRate.tsv"
	run:
		get_het_sample_outliers(input[0], output[0])

#we may need data from varcall pipeline
#or we can use Read Depth by Sample
rule sampleDP:
	output:
		expand(os.path.join(BASE_OUT,config.get("rules").get("coverage").get("out_dir"), "{{vcf_name}}_dp.{ext}"), ext=["idepth"])
	input:
		# vcf=os.path.join(BASE_OUT,config.get("rules").get("mergeReapplyVQSR").get("out_dir"),"{vcf_name}.vcf.gz"),
		# vcf_index=os.path.join(BASE_OUT,config.get("rules").get("mergeReapplyVQSR").get("out_dir"),"{vcf_name}.vcf.gz.tbi")
		vcf=rules.cleanMissingHwe.output[0],
		vcf_index=rules.cleanMissingHwe.output[0]
	params:
		bcftools=config['BCFTOOLS'],
		vcftools=config['VCFTOOLS'],
		tmp=os.path.join(BASE_OUT,config.get("paths").get("tmp")),
		out_prefix=os.path.join(BASE_OUT,config.get("rules").get("coverage").get("out_dir"), "{vcf_name}_dp")
	log:
		config["paths"]["log_dir"] + "/{vcf_name}-sampleDP.log",
		config["paths"]["log_dir"] + "/{vcf_name}-sampleDP.e"
	threads: 1
	resources:
		mem_mb=5000
	benchmark:
		config["paths"]["benchmark"] + "/{vcf_name}_sampleDP.tsv"
	envmodules:
		"vcftools/0.1.16"
	shell:
		"""
		{params.vcftools} --gzvcf {input.vcf} --depth --out {params.out_prefix} 1> {log[0]} 2> {log[1]}
		"""

# rule coverage:
# 	output:
# 	input:
# 	params:
# 	log:
# 	threads:
# 	resources:
# 	benchmark:
# 	shell:


#Missing rate rule
rule SampleMissingRate:
	output:
		expand(os.path.join(BASE_OUT,config.get("rules").get("SampleMissingRate").get("out_dir"), "{{vcf_name}}_missing.{ext}"), ext=["imiss"])
	input:
		# vcf=os.path.join(BASE_OUT,config.get("rules").get("mergeReapplyVQSR").get("out_dir"),"{vcf_name}.vcf.gz"),
		# vcf_index=os.path.join(BASE_OUT,config.get("rules").get("mergeReapplyVQSR").get("out_dir"),"{vcf_name}.vcf.gz.tbi")
		vcf=rules.cleanMissingHwe.output[0],
		vcf_index=rules.cleanMissingHwe.output[0]		
	params:
		bcftools=config['BCFTOOLS'],
		vcftools=config['VCFTOOLS'],
		tmp=os.path.join(BASE_OUT,config.get("paths").get("tmp")),
		out_prefix=os.path.join(BASE_OUT,config.get("rules").get("SampleMissingRate").get("out_dir"), "{vcf_name}_missing")
	log:
		config["paths"]["log_dir"] + "/{vcf_name}-sampleMissing.log",
		config["paths"]["log_dir"] + "/{vcf_name}-sampleMissing.e"
	threads: 1
	resources:
		mem_mb=5000
	benchmark:
		config["paths"]["benchmark"] + "/{vcf_name}_sampleMissing.tsv"
	envmodules:
		"vcftools/0.1.16"
	shell:
		"""
		{params.vcftools} --gzvcf {input.vcf} --missing-indv --out {params.out_prefix} 1> {log[0]} 2> {log[1]}
		"""
