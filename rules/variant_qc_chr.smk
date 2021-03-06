#this set of rules is devised in two sections:
# 1) the first rule will perform the automatic removal of all sites with missing rate and hwe values exceeding some thresholds
rule cleanMissingHwe:
	output:
		os.path.join(BASE_OUT,config.get("rules").get("cleanMissingHwe").get("out_dir"), "{vcf_name}_{chr}_HWE95call.vcf.gz"),
		os.path.join(BASE_OUT,config.get("rules").get("cleanMissingHwe").get("out_dir"), "{vcf_name}_{chr}_HWE95call.vcf.gz.tbi"),
	input:
		vcf=os.path.join(BASE_OUT,config.get("rules").get("mergeReapplyVQSR").get("out_dir"),"{vcf_name}_VQSLODrefilter.vcf.gz"),
		vcf_index=os.path.join(BASE_OUT,config.get("rules").get("mergeReapplyVQSR").get("out_dir"),"{vcf_name}_VQSLODrefilter.vcf.gz.tbi")
	params:
		bcftools=config['BCFTOOLS'],
		hwe_thr=config.get("rules").get("cleanMissingHwe").get("hwe_thr"),
		missing_thr=config.get("rules").get("cleanMissingHwe").get("missing_thr"),
		# out_prefix=os.path.join(BASE_OUT,config.get("rules").get("cleanMissingHwe").get("out_dir"), "{vcf_name}_HWE95call")
	log:
		config["paths"]["log_dir"] + "/{vcf_name}-{chr}-cleanMissingHwe.log",
		config["paths"]["log_dir"] + "/{vcf_name}-{chr}-cleanMissingHwe.e"
	threads: 2
	resources:
		mem_mb=5000
	benchmark:
		config["paths"]["benchmark"] + "/{vcf_name}_{chr}_cleanMissingHwe.tsv"
	envmodules:
		"vcftools/0.1.16",
		"bcftools/1.14"
	shell:
		"""
		({params.bcftools} view -r {wildcards.chr} {input.vcf}|{params.bcftools} +fill-tags -- -t all,F_MISSING,HWE | bcftools view -e "HWE < {params.hwe_thr} | F_MISSING > {params.missing_thr}" -O z -o {output[0]}) 1> {log[0]} 2> {log[1]}
		{params.bcftools} index -t {output[0]} 1>> {log[0]} 2>> {log[1]}
		"""

rule cleanMissingHweList:
	output:
		os.path.join(BASE_OUT,config.get("rules").get("cleanMissingHwe").get("out_dir"), "{vcf_name}_{chr}_HWE95call.removed.sites")
	input:
		vcf=os.path.join(BASE_OUT,config.get("rules").get("mergeReapplyVQSR").get("out_dir"),"{vcf_name}_VQSLODrefilter.vcf.gz"),
		vcf_index=os.path.join(BASE_OUT,config.get("rules").get("mergeReapplyVQSR").get("out_dir"),"{vcf_name}_VQSLODrefilter.vcf.gz.tbi")
	params:
		vcftools=config['VCFTOOLS'],
		hwe_thr=config.get("rules").get("cleanMissingHwe").get("hwe_thr"),
		missing_thr=config.get("rules").get("cleanMissingHwe").get("missing_thr"),
		# out_prefix=os.path.join(BASE_OUT,config.get("rules").get("cleanMissingHwe").get("out_dir"), "{vcf_name}_HWE95call")
	log:
		config["paths"]["log_dir"] + "/{vcf_name}-{chr}-cleanMissingHweList.log",
		config["paths"]["log_dir"] + "/{vcf_name}-{chr}-cleanMissingHweList.e"
	threads: 1
	resources:
		mem_mb=5000
	benchmark:
		config["paths"]["benchmark"] + "/{vcf_name}_{chr}_cleanMissingHweList.tsv"
	envmodules:
		"bcftools/1.14"
	shell:
		"""
		({params.bcftools} view -r {wildcards.chr} {input.vcf} |{params.bcftools} +fill-tags -- -t all,F_MISSING,HWE | {params.bcftools} view -i "HWE < {params.hwe_thr} | F_MISSING > {params.missing_thr}" | {params.bcftools} view -G -O z -o {output[0]}) 1> {log[0]} 2> {log[1]}
		"""


# 2) set of rules to calculate het rate and missing rate lists to be used as summary data and removal lists
#het rate rule: first get the data with vcftools
rule VariantsHetRate:
	output:
		os.path.join(BASE_OUT,config.get("rules").get("VariantsHetRate").get("out_dir"), "{vcf_name}_{chr}_hwe.hwe")
	input:
		vcf=rules.cleanMissingHwe.output[0],
		vcf_index=rules.cleanMissingHwe.output[1]
	params:
		vcftools=config['VCFTOOLS'],
		out_prefix=os.path.join(BASE_OUT,config.get("rules").get("VariantsHetRate").get("out_dir"), "{vcf_name}_{chr}_hwe")
	log:
		config["paths"]["log_dir"] + "/{vcf_name}-{chr}-hwe.log",
		config["paths"]["log_dir"] + "/{vcf_name}-{chr}-hwe.e"
	threads: 1
	resources:
		mem_mb=5000
	benchmark:
		config["paths"]["benchmark"] + "/{vcf_name}_{chr}_hwe.tsv"
	envmodules:
		"vcftools/0.1.16"
	shell:
		"""
		{params.vcftools} --gzvcf {input.vcf} --chr {wildcards.chr} --hardy --out {params.out_prefix} 1> {log[0]} 2> {log[1]}
		"""

#het rate rule: get vcftools result and extract the het rate for plotting
rule VariantsGetHetRateOut:
	output:
		os.path.join(BASE_OUT,config.get("rules").get("VariantsHetRate").get("out_dir"), "{vcf_name}_{chr}_ToRemHetRate.txt")
	input:
		rules.VariantsHetRate.output[0]
	params:
		vcftools=config['VCFTOOLS'],
		exc_het_thr=float(config.get('rules').get('VariantsGetHetRateOut').get('exc_het_pval_thr'))
	log:
		config["paths"]["log_dir"] + "/{vcf_name}-{chr}-toRemHet.log",
		config["paths"]["log_dir"] + "/{vcf_name}-{chr}-toRemHet.e"
	threads: 1
	resources:
		mem_mb=5000
	benchmark:
		config["paths"]["benchmark"] + "/{vcf_name}_{chr}_toRemHet.tsv"
	run:
		logger = logging.getLogger('logging_test')
		fh = logging.FileHandler(str(log[1]))
		fh.setLevel(logging.INFO)
		formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
		fh.setFormatter(formatter)
		logger.addHandler(fh)
		try: 
			logger.info('Starting operation!')
			# do something
			get_het_hwe_variants_outliers(input[0],params.exc_het_thr,output[0])
			logger.info('Ended!')
		except Exception as e: 
			logger.error(e, exc_info=True)

#Missing rate rule
rule VariantsMissingRate:
	output:
		os.path.join(BASE_OUT,config.get("rules").get("VariantsMissingRate").get("out_dir"), "{vcf_name}_{chr}_missing.lmiss")
	input:
		vcf=rules.cleanMissingHwe.output[0],
		vcf_index=rules.cleanMissingHwe.output[1]
	params:
		bcftools=config['BCFTOOLS'],
		vcftools=config['VCFTOOLS'],
		tmp=os.path.join(BASE_OUT,config.get("paths").get("tmp")),
		out_prefix=os.path.join(BASE_OUT,config.get("rules").get("VariantsMissingRate").get("out_dir"), "{vcf_name}_{chr}_missing")
	log:
		config["paths"]["log_dir"] + "/{vcf_name}-{chr}-variantsMissing.log",
		config["paths"]["log_dir"] + "/{vcf_name}-{chr}-variantsMissing.e"
	threads: 1
	resources:
		mem_mb=5000
	benchmark:
		config["paths"]["benchmark"] + "/{vcf_name}_{chr}_variantsMissing.tsv"
	envmodules:
		"vcftools/0.1.16"
	shell:
		"""
		{params.vcftools} --gzvcf {input.vcf} --chr {wildcards.chr} --missing-site --out {params.out_prefix} 1> {log[0]} 2> {log[1]}
		"""

#extract AF from the  study population data
rule getPopAF:
	output:
		os.path.join(BASE_OUT,config.get("rules").get("getPopAF").get("out_dir"), "{vcf_name}_{chr}_af.txt")
	input:
		vcf=rules.cleanMissingHwe.output[0],
		vcf_index=rules.cleanMissingHwe.output[1]
	params:
		bcftools=config['BCFTOOLS']
	log:
		config["paths"]["log_dir"] + "/{vcf_name}-{chr}-getPopAF.log",
		config["paths"]["log_dir"] + "/{vcf_name}-{chr}-getPopAF.e"
	threads: 2
	resources:
		mem_mb=5000
	benchmark:
		config["paths"]["benchmark"] + "/{vcf_name}_{chr}_getPopAF.tsv"
	envmodules:
		"bcftools/1.14"
	shell:
		"""
		({params.bcftools} +fill-tags {input.vcf} -- -t all | {params.bcftools} query -f "%CHROM\_%POS\_%REF\_%ALT\t%CHROM\t%POS\t%ID\t%REF\t%ALT\t%AC\t%AN\t%AF\t%MAF\n" -o {output}) 1> {log[0]} 2> {log[1]}
		"""

#merge pop data with data from TGP and EUR only data
rule comparePopAF:
	output:
		os.path.join(BASE_OUT,config.get("rules").get("comparePopAF").get("out_dir"), "{vcf_name}_{ref_pop}_{chr}_af_extrDiff.txt"),
		os.path.join(BASE_OUT,config.get("rules").get("comparePopAF").get("out_dir"), "{vcf_name}_{ref_pop}_{chr}_af.png"),
		os.path.join(BASE_OUT,config.get("rules").get("comparePopAF").get("out_dir"), "{vcf_name}_{ref_pop}_{chr}_af_extrDiff.pdf")
	input:
		wgs_table=rules.getPopAF.output[0]
	params:
		ext_table=lambda wildcards: config.get("rules").get("comparePopAF").get("ref_pops").get(wildcards.ref_pop),
		out_prefix=os.path.join(BASE_OUT,config.get("rules").get("comparePopAF").get("out_dir"))
	log:
		config["paths"]["log_dir"] + "/{vcf_name}-{ref_pop}-{chr}-comparePopAF.log",
		config["paths"]["log_dir"] + "/{vcf_name}-{ref_pop}-{chr}-comparePopAF.e"
	threads: 1
	resources:
		mem_mb=10000
	benchmark:
		config["paths"]["benchmark"] + "/{vcf_name}_{ref_pop}_{chr}_comparePopAF.tsv"
	run:
		outname_tab=params.out_prefix + "/"+ wildcards.vcf_name + "_" + wildcards.ref_pop + "_" + wildcards.chr +"_af_extrDiff.txt"
		outname_plot=params.out_prefix + "/"+ wildcards.vcf_name + "_" + wildcards.ref_pop + "_" + wildcards.chr +"_af.png"
		logger = logging.getLogger('logging_test')
		fh = logging.FileHandler(str(log[1]))
		fh.setLevel(logging.INFO)
		formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
		fh.setFormatter(formatter)
		logger.addHandler(fh)
		try: 
			logger.info('Starting operation!')
			# do something
			# af_diff(input.wgs_table, params.ext_table, outname_tab, outname_plot,output[2])
			af_diff(input.wgs_table, params.ext_table, output[0], output[1],output[2])
			logger.info('Ended!')
		except Exception as e: 
			logger.error(e, exc_info=True)
