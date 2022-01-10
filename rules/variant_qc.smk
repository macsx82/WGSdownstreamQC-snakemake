#this set of rules is devised in two sections:
# 1) the first rule will perform the automatic removal of all sites with missing rate and hwe values exceeding some thresholds
rule cleanMissingHwe:
	output:
		os.path.join(BASE_OUT,config.get("rules").get("cleanMissingHwe").get("out_dir"), "{vcf_name}_HWE95call.vcf.gz"),
		os.path.join(BASE_OUT,config.get("rules").get("cleanMissingHwe").get("out_dir"), "{vcf_name}_HWE95call.vcf.gz.tbi"),
		os.path.join(BASE_OUT,config.get("rules").get("cleanMissingHwe").get("out_dir"), "{vcf_name}_HWE95call.removed.sites"),
		os.path.join(BASE_OUT,config.get("rules").get("cleanMissingHwe").get("out_dir"), "{vcf_name}_HWE95call.log")
	input:
		vcf=os.path.join(BASE_OUT,config.get("rules").get("mergeReapplyVQSR").get("out_dir"),"{vcf_name}_VQSLODrefilter.vcf.gz"),
		vcf_index=os.path.join(BASE_OUT,config.get("rules").get("mergeReapplyVQSR").get("out_dir"),"{vcf_name}_VQSLODrefilter.vcf.gz.tbi")
	params:
		vcftools=config['VCFTOOLS'],
		bcftools=config['BCFTOOLS'],
		hwe_thr=config.get("rules").get("cleanMissingHwe").get("hwe_thr"),
		missing_thr=config.get("rules").get("cleanMissingHwe").get("missing_thr"),
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
		({params.vcftools} --gzvcf {input.vcf} --hwe {params.hwe_thr} --max-missing {params.missing_thr} --recode --recode-INFO-all -c | {params.bcftools} view -O z -o {output[0]}) 1> {log[0]} 2> {log[1]}
		({params.vcftools} --gzvcf {input.vcf} --hwe {params.hwe_thr} --max-missing {params.missing_thr} --removed-sites --out {params.out_prefix}) 1>> {log[0]} 2>> {log[1]}
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

# #het rate rule: get vcftools result and extract the het rate for plotting
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
rule VariantsMissingRate:
	output:
		expand(os.path.join(BASE_OUT,config.get("rules").get("VariantsMissingRate").get("out_dir"), "{{vcf_name}}_missing.{ext}"), ext=["lmiss", "log"])
	input:
		# vcf=os.path.join(BASE_OUT,config.get("rules").get("mergeReapplyVQSR").get("out_dir"),"{vcf_name}.vcf.gz"),
		# vcf_index=os.path.join(BASE_OUT,config.get("rules").get("mergeReapplyVQSR").get("out_dir"),"{vcf_name}.vcf.gz.tbi")
		vcf=rules.cleanMissingHwe.output[0],
		vcf_index=rules.cleanMissingHwe.output[0]
	params:
		bcftools=config['BCFTOOLS'],
		vcftools=config['VCFTOOLS'],
		tmp=os.path.join(BASE_OUT,config.get("paths").get("tmp")),
		out_prefix=os.path.join(BASE_OUT,config.get("rules").get("VariantsMissingRate").get("out_dir"), "{vcf_name}_missing")
	log:
		config["paths"]["log_dir"] + "/{vcf_name}-variantsMissing.log",
		config["paths"]["log_dir"] + "/{vcf_name}-variantsMissing.e"
	threads: 1
	resources:
		mem_mb=5000
	benchmark:
		config["paths"]["benchmark"] + "/{vcf_name}_variantsMissing.tsv"
	envmodules:
		"vcftools/0.1.16"
	shell:
		"""
		{params.vcftools} --gzvcf {input.vcf} --missing-site --out {params.out_prefix} 1> {log[0]} 2> {log[1]}
		"""

#extract AF from the  study population data
rule getPopAF:
	output:
		os.path.join(BASE_OUT,config.get("rules").get("getPopAF").get("out_dir"), "{vcf_name}_af.txt")
	input:
		vcf=rules.cleanMissingHwe.output[0],
		vcf_index=rules.cleanMissingHwe.output[0]
	params:
		bcftools=config['BCFTOOLS']
	log:
		config["paths"]["log_dir"] + "/{vcf_name}-getPopAF.log",
		config["paths"]["log_dir"] + "/{vcf_name}-getPopAF.e"
	threads: 2
	resources:
		mem_mb=5000
	benchmark:
		config["paths"]["benchmark"] + "/{vcfname}_getPopAF.tsv"
	envmodules:
		"bcftools/1.14"
	shell:
		"""
		({params.bcftools} +fill-tags {input.vcf} -- -t all | {params.bcftools} query -f "%CHROM\_%POS\_%REF\_%ALT\t%CHROM\t%POS\t%ID\t%REF\t%ALT\t%AC\t%AN\t%AF\t%MAF\n" -o {output}) 1> {log[0]} 2> {log[1]}
		"""

#merge pop data with data from TGP and EUR only data
rule comparePopAF:
	output:
		expand(os.path.join(BASE_OUT,config.get("rules").get("comparePopAF").get("out_dir"), "{{vcf_name}}_{ext_ref}_af_extrDiff.txt"), ext_ref=list(config.get("rules").get("comparePopAF").get("ref_pops").keys())),
		expand(os.path.join(BASE_OUT,config.get("rules").get("comparePopAF").get("out_dir"), "{{vcf_name}}_{ext_ref}_af.pdf"), ext_ref=list(config.get("rules").get("comparePopAF").get("ref_pops").keys()))
	input:
		wgs_table=rules.getPopAF.output[0]
	params:
		ext_tables=config.get("rules").get("comparePopAF").get("ref_pops"),
		out_prefix=os.path.join(BASE_OUT,config.get("rules").get("comparePopAF").get("out_dir")),
		# out_tab=os.path.join(BASE_OUT,config.get("rules").get("comparePopAF").get("out_dir"), "{{vcf_name}}_{ext_ref}_af_extrDiff.txt"),
		# out_plot=os.path.join(BASE_OUT,config.get("rules").get("comparePopAF").get("out_dir"), "{{vcf_name}}_{ext_ref}_af.pdf")
	log:
		config["paths"]["log_dir"] + "/{vcf_name}_{ext_pop}-comparePopAF.log",
		config["paths"]["log_dir"] + "/{vcf_name}_{ext_pop}-comparePopAF.e"
	threads: 1
	resources:
		mem_mb=10000
	benchmark:
		config["paths"]["benchmark"] + "/{vcfname}_{ext_pop}_comparePopAF.tsv"
	envmodules:
		"bcftools/1.14"
	run:
		for ext_table in ext_tables.keys():
			outname_tab=params.out_prefix + "/{vcf_name}_" + ext_table + "_af_extrDiff.txt"
			outname_plot=params.out_prefix + "/{vcf_name}_" + ext_table + "_af.pdf"
			af_diff(input.wgs_table, ext_table, outname_tab, outname_plot)