#this set of rules is devised in two sections:
# 1) the first rule will perform the automatic removal of all sites with missing rate and hwe values exceeding some thresholds
rule cleanMissingHwe:
	output:
		os.path.join(BASE_OUT,config.get("rules").get("cleanMissingHwe").get("out_dir"), "{vcf_name}_HWE95call.vcf.gz"),
		os.path.join(BASE_OUT,config.get("rules").get("cleanMissingHwe").get("out_dir"), "{vcf_name}_HWE95call.vcf.gz.tbi"),
		os.path.join(BASE_OUT,config.get("rules").get("cleanMissingHwe").get("out_dir"), "{vcf_name}_HWE95call.removed.sites"),
		os.path.join(BASE_OUT,config.get("rules").get("cleanMissingHwe").get("out_dir"), "{vcf_name}_HWE95call.log")

	input:
		vcf=os.path.join(BASE_OUT,config.get("rules").get("mergeReapplyVQSR").get("out_dir"),"{vcf_name}.vcf.gz"),
		vcf_index=os.path.join(BASE_OUT,config.get("rules").get("mergeReapplyVQSR").get("out_dir"),"{vcf_name}.vcf.gz.tbi")
	params:
		vcftools=config['VCFTOOLS'],
		bcftools=config['BCFTOOLS'],
		hwe_thr=gonfig.get(rules).get("cleanMissingHwe").get("hwe_thr"),
		missing_thr=config.get(rules).get("cleanMissingHwe").get("missing_thr"),
		out_prefix=os.path.join(BASE_OUT,config.get("rules").get("cleanMissingHwe").get("out_dir"), "{vcf_name}_HWE95call")
	log:
		config["paths"]["log_dir"] + "/{vcf_name}-cleanMissingHwe.log",
		config["paths"]["log_dir"] + "/{vcf_name}-cleanMissingHwe.e"
	threads: 1
	resources:
		mem_mb=5000
	benchmark:
		config["paths"]["benchmark"] + "/{vcf_name}_cleanMissingHwe.tsv"
	envmodules:
		"vcftools/0.1.16",
		"bcftools/1.14"
	shell:
		"""
		{params.vcftools} --gzvcf {input.vcf} --hwe {params.hwe_thr} --max-missing {params.missing_thr} --recode --recode-INFO-all -c | {params.bcftools} view -O z -o {output[0]} 1> {log[0]} 2> {log[1]}
		{params.vcftools} --gzvcf {input.vcf} --hwe {params.hwe_thr} --max-missing {params.missing_thr} --removed-sites --out {params.out_prefix} 1>> {log[0]} 2>> {log[1]}
		{params.bcftools} index -t {output[0]}
		"""

# 2) set of rules to calculate het rate and missing rate lists to be used as summary data and removal lists
#het rate rule: first get the data with vcftools
rule VariantsHetRate:
	output:
		expand(os.path.join(BASE_OUT,config.get("rules").get("VariantsHetRate").get("out_dir"), "{{vcf_name}}_hwe.{ext}"), ext=["hwe", "log"])
	input:
		# vcf=os.path.join(BASE_OUT,config.get("rules").get("mergeReapplyVQSR").get("out_dir"),"{vcf_name}.vcf.gz"),
		# vcf_index=os.path.join(BASE_OUT,config.get("rules").get("mergeReapplyVQSR").get("out_dir"),"{vcf_name}.vcf.gz.tbi")
		vcf=rules.cleanMissingHwe.output[0],
		vcf_index=rules.cleanMissingHwe.output[0]
	params:
		vcftools=config['VCFTOOLS'],
		out_prefix=os.path.join(BASE_OUT,config.get("rules").get("VariantsHetRate").get("out_dir"), "{vcf_name}_hwe")
	log:
		config["paths"]["log_dir"] + "/{vcf_name}-hwe.log",
		config["paths"]["log_dir"] + "/{vcf_name}-hwe.e"
	threads: 1
	resources:
		mem_mb=5000
	benchmark:
		config["paths"]["benchmark"] + "/{vcf_name}_hwe.tsv"
	envmodules:
		"vcftools/0.1.16"
	shell:
		"""
		{params.vcftools} --gzvcf {input.vcf} --hardy --out {params.out_prefix} 1> {log[0]} 2> {log[1]}
		"""

#het rate rule: get vcftools result and extract the het rate for plotting
# rule VariantsGetHetRateOut:
# 	output:
# 		os.path.join(BASE_OUT,config.get("rules").get("VariantsHetRate").get("out_dir"), "{vcf_name}_hetRate.txt")
# 	input:
# 		rules.VariantsHetRate.output[0]
# 	params:
# 		vcftools=config['VCFTOOLS']
# 	log:
# 		config["paths"]["log_dir"] + "/{vcf_name}-hetHwe.log",
# 		config["paths"]["log_dir"] + "/{vcf_name}-hetHwe.e"
# 	threads: 1
# 	resources:
# 		mem_mb=5000
# 	benchmark:
# 		config["paths"]["benchmark"] + "/{vcf_name}_hetHwe.tsv"
# 	run:
# 		get_het_sample_outliers(input[0], output[0])

#Missing rate rule
# rule VariantsMissingRate:
# 	output:
# 		expand(os.path.join(BASE_OUT,config.get("rules").get("SampleMissingRate").get("out_dir"), "{{vcf_name}}_missing.{ext}"), ext=["imiss", "log"])
# 	input:
# 		vcf=os.path.join(BASE_OUT,config.get("rules").get("mergeReapplyVQSR").get("out_dir"),"{vcf_name}.vcf.gz"),
# 		vcf_index=os.path.join(BASE_OUT,config.get("rules").get("mergeReapplyVQSR").get("out_dir"),"{vcf_name}.vcf.gz.tbi")
# 	params:
# 		bcftools=config['BCFTOOLS'],
# 		vcftools=config['VCFTOOLS'],
# 		tmp=os.path.join(BASE_OUT,config.get("paths").get("tmp")),
# 		out_prefix=os.path.join(BASE_OUT,config.get("rules").get("SampleMissingRate").get("out_dir"), "{vcf_name}_missing")
# 	log:
# 		config["paths"]["log_dir"] + "/{vcf_name}-sampleMissing.log",
# 		config["paths"]["log_dir"] + "/{vcf_name}-sampleMissing.e"
# 	threads: 1
# 	resources:
# 		mem_mb=5000
# 	benchmark:
# 		config["paths"]["benchmark"] + "/{vcf_name}_sampleMissing.tsv"
# 	envmodules:
# 		"vcftools/0.1.16"
# 	shell:
# 		"""
# 		{params.vcftools} --gzvcf {input.vcf} --missing-indv --out {params.out_prefix} 1> {log[0]} 2> {log[1]}
# 		"""