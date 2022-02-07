#this set of rules is devised in two sections:
# 1) the first rule will perform the automatic removal of all sites with missing rate and hwe values exceeding some thresholds
rule cleanMissingHwe:
	output:
		os.path.join(BASE_OUT,config.get("rules").get("cleanMissingHwe").get("out_dir"), "{vcf_name}_HWE95call.vcf.gz"),
		os.path.join(BASE_OUT,config.get("rules").get("cleanMissingHwe").get("out_dir"), "{vcf_name}_HWE95call.vcf.gz.tbi"),
	input:
		vcf=os.path.join(BASE_OUT,config.get("rules").get("mergeReapplyVQSR").get("out_dir"),"{vcf_name}_VQSLODrefilter.vcf.gz"),
		vcf_index=os.path.join(BASE_OUT,config.get("rules").get("mergeReapplyVQSR").get("out_dir"),"{vcf_name}_VQSLODrefilter.vcf.gz.tbi")
	params:
		bcftools=config['BCFTOOLS'],
		hwe_thr=config.get("rules").get("cleanMissingHwe").get("hwe_thr"),
		missing_thr=config.get("rules").get("cleanMissingHwe").get("missing_thr"),
		# out_prefix=os.path.join(BASE_OUT,config.get("rules").get("cleanMissingHwe").get("out_dir"), "{vcf_name}_HWE95call")
	log:
		config["paths"]["log_dir"] + "/{vcf_name}-cleanMissingHwe.log",
		config["paths"]["log_dir"] + "/{vcf_name}-cleanMissingHwe.e"
	threads: 3
	resources:
		mem_mb=5000
	benchmark:
		config["paths"]["benchmark"] + "/{vcf_name}_cleanMissingHwe.tsv"
	envmodules:
		"bcftools/1.14"
	shell:
		"""
		({params.bcftools} +fill-tags {input.vcf} -- -t all,F_MISSING,HWE | {params.bcftools} view -e "HWE < {params.hwe_thr}"| {params.bcftools} view -e "F_MISSING > {params.missing_thr}" -O z -o {output[0]}) 1> {log[0]} 2> {log[1]}
		{params.bcftools} index -t {output[0]} 1>> {log[0]} 2>> {log[1]}
		"""


rule cleanMissingHweList:
	output:
		os.path.join(BASE_OUT,config.get("rules").get("cleanMissingHwe").get("out_dir"), "{vcf_name}_HWE95call.removed.sites.vcf.gz")
	input:
		vcf=os.path.join(BASE_OUT,config.get("rules").get("mergeReapplyVQSR").get("out_dir"),"{vcf_name}_VQSLODrefilter.vcf.gz"),
		vcf_index=os.path.join(BASE_OUT,config.get("rules").get("mergeReapplyVQSR").get("out_dir"),"{vcf_name}_VQSLODrefilter.vcf.gz.tbi")
	params:
		bcftools=config['BCFTOOLS'],
		hwe_thr=config.get("rules").get("cleanMissingHwe").get("hwe_thr"),
		missing_thr=config.get("rules").get("cleanMissingHwe").get("missing_thr"),
		# out_prefix=os.path.join(BASE_OUT,config.get("rules").get("cleanMissingHwe").get("out_dir"), "{vcf_name}_HWE95call")
	log:
		config["paths"]["log_dir"] + "/{vcf_name}-cleanMissingHweList.log",
		config["paths"]["log_dir"] + "/{vcf_name}-cleanMissingHweList.e"
	threads: 1
	resources:
		mem_mb=5000
	benchmark:
		config["paths"]["benchmark"] + "/{vcf_name}_cleanMissingHweList.tsv"
	envmodules:
		"bcftools/1.14"
	shell:
		"""
		({params.bcftools} +fill-tags {input.vcf} -- -t all,F_MISSING,HWE | {params.bcftools} view -e "HWE < {params.hwe_thr}"| {params.bcftools} view -G -e "F_MISSING > {params.missing_thr}" -O z -o {output[0]}) 1> {log[0]} 2> {log[1]}
		"""


# 2) set of rules to calculate het rate and missing rate lists to be used as summary data and removal lists
#het rate rule: first get the data with vcftools
rule VariantsHetRate:
	output:
		os.path.join(BASE_OUT,config.get("rules").get("VariantsHetRate").get("out_dir"), "{vcf_name}_hwe.hwe")
	input:
		vcf=rules.cleanMissingHwe.output[0],
		vcf_index=rules.cleanMissingHwe.output[1]
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
		"bcftools/1.14"
	shell:
		"""
		(echo -e "CHROM\tPOS\tID\tREF\tALT\tP_HET_EXCESS";{params.bcftools} +fill-tags {input.vcf} -- -t all,ExcHet | {params.bcftools} query -f "%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO/ExcHet\n") > {output[0]} 1> {log[0]} 2> {log[1]}
		"""

#rule to get the outliers variants to remove, based on the excess of heterozygosity pvalue calculated with vcftools
# we set the exclusion threshold in the config file
rule VariantsGetHetRateOut:
	output:
		os.path.join(BASE_OUT,config.get("rules").get("VariantsHetRate").get("out_dir"), "{vcf_name}_ToRemHetRate.txt")
	input:
		rules.VariantsHetRate.output[0]
	params:
		vcftools=config['VCFTOOLS'],
		exc_het_thr=float(config.get('rules').get('VariantsGetHetRateOut').get('exc_het_pval_thr'))
	log:
		config["paths"]["log_dir"] + "/{vcf_name}-toRemHet.log",
		config["paths"]["log_dir"] + "/{vcf_name}-toRemHet.e"
	threads: 1
	resources:
		mem_mb=5000
	benchmark:
		config["paths"]["benchmark"] + "/{vcf_name}toRemHet.tsv"
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
		os.path.join(BASE_OUT,config.get("rules").get("VariantsMissingRate").get("out_dir"), "{vcf_name}_missing.lmiss")
	input:
		vcf=rules.cleanMissingHwe.output[0],
		vcf_index=rules.cleanMissingHwe.output[1]
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
		"bcftools/1.14"
	shell:
		"""
		samples_n=$({params.bcftools} query -l {input.vcf})
		(echo -e "CHR\tPOS\tID\tREF\tALT\tF_MISS\tN_MISS\tN_DATA";{params.bcftools} +fill-tags {input.vcf} -- -t all,NS,F_MISSING | {params.bcftools} query -f "%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO/F_MISSING\t%INFO/AN\n"| awk -v snum=${{samples_n}} '{FS="\t";OFS="\t"}{print $1,$2,$3,$4,$5,$6,snum*2-$7,snum*2}') > {output[0]} 1> {log[0]} 2> {log[1]}
		"""
		# {params.vcftools} --gzvcf {input.vcf} --missing-site --out {params.out_prefix} 1> {log[0]} 2> {log[1]}

#extract AF from the  study population data
rule getPopAF:
	output:
		os.path.join(BASE_OUT,config.get("rules").get("getPopAF").get("out_dir"), "{vcf_name}_af.txt")
	input:
		vcf=rules.cleanMissingHwe.output[0],
		vcf_index=rules.cleanMissingHwe.output[1]
	params:
		bcftools=config['BCFTOOLS']
	log:
		config["paths"]["log_dir"] + "/{vcf_name}-getPopAF.log",
		config["paths"]["log_dir"] + "/{vcf_name}-getPopAF.e"
	threads: 2
	resources:
		mem_mb=5000
	benchmark:
		config["paths"]["benchmark"] + "/{vcf_name}_getPopAF.tsv"
	envmodules:
		"bcftools/1.14"
	shell:
		"""
		({params.bcftools} +fill-tags {input.vcf} -- -t all | {params.bcftools} query -f "%CHROM\_%POS\_%REF\_%ALT\t%CHROM\t%POS\t%ID\t%REF\t%ALT\t%AC\t%AN\t%AF\t%MAF\n" -o {output}) 1> {log[0]} 2> {log[1]}
		"""

#merge pop data with data from TGP and EUR only data
rule comparePopAF:
	output:
		os.path.join(BASE_OUT,config.get("rules").get("comparePopAF").get("out_dir"), "{vcf_name}_{ref_pop}_af_extrDiff.txt"),
		os.path.join(BASE_OUT,config.get("rules").get("comparePopAF").get("out_dir"), "{vcf_name}_{ref_pop}_af.png"),
		os.path.join(BASE_OUT,config.get("rules").get("comparePopAF").get("out_dir"), "{vcf_name}_{ref_pop}_af_extrDiff.pdf")
	input:
		wgs_table=rules.getPopAF.output[0]
	params:
		ext_table=lambda wildcards: config.get("rules").get("comparePopAF").get("ref_pops").get(wildcards.ref_pop),
		out_prefix=os.path.join(BASE_OUT,config.get("rules").get("comparePopAF").get("out_dir"))
	log:
		config["paths"]["log_dir"] + "/{vcf_name}_{ref_pop}-comparePopAF.log",
		config["paths"]["log_dir"] + "/{vcf_name}_{ref_pop}-comparePopAF.e"
	threads: 1
	resources:
		mem_mb=10000
	benchmark:
		config["paths"]["benchmark"] + "/{vcf_name}_{ref_pop}_comparePopAF.tsv"
	run:
		outname_tab=params.out_prefix + "/"+ wildcards.vcf_name + "_" + wildcards.ref_pop + "_af_extrDiff.txt"
		outname_plot=params.out_prefix + "/"+ wildcards.vcf_name + "_" + wildcards.ref_pop + "_af.png"
		logger = logging.getLogger('logging_test')
		fh = logging.FileHandler(str(log[1]))
		fh.setLevel(logging.INFO)
		formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
		fh.setFormatter(formatter)
		logger.addHandler(fh)
		try: 
			logger.info('Starting operation!')
			# do something
			af_diff(input.wgs_table, params.ext_table, outname_tab, outname_plot,output[2])
			logger.info('Ended!')
		except Exception as e: 
			logger.error(e, exc_info=True)
