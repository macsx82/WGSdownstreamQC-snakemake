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
		logger = logging.getLogger('logging_test')
		fh = logging.FileHandler(str(log[1]))
		fh.setLevel(logging.DEBUG)
		formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
		fh.setFormatter(formatter)
		logger.addHandler(fh)
		try: 
			logger.info('Starting operation!')
			# do something
			get_het_sample_outliers(input[0], output[0])
			logger.info('Ended!')
		except Exception as e: 
			logger.error(e, exc_info=True)

#plot het rate per sample
rule PlotHetRateSample:
	output:
		os.path.join(BASE_OUT,config.get("rules").get("SampleHetRate").get("out_dir"), "{vcf_name}_hetRate.pdf")
	input:
		rules.SampleGetHetRateOut.output[0]
	params:
		vcftools=config['VCFTOOLS']
	log:
		config["paths"]["log_dir"] + "/{vcf_name}-PlotHetRateSample.log",
		config["paths"]["log_dir"] + "/{vcf_name}-PlotHetRateSample.e"
	threads: 1
	resources:
		mem_mb=5000
	benchmark:
		config["paths"]["benchmark"] + "/{vcf_name}_PlotHetRateSample.tsv"
	run:
		logger = logging.getLogger('logging_test')
		fh = logging.FileHandler(str(log[1]))
		fh.setLevel(logging.DEBUG)
		formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
		fh.setFormatter(formatter)
		logger.addHandler(fh)
		try: 
			logger.info('Starting operation!')
			# do something
			plot_het_rate_sample(input[0], output[0])
			logger.info('Ended!')
		except Exception as e: 
			logger.error(e, exc_info=True)

#we may need data from varcall pipeline
#or we can use Read Depth by Sample
rule sampleDP:
	output:
		os.path.join(BASE_OUT,config.get("rules").get("coverage").get("out_dir"), "{vcf_name}_dp.idepth")
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

#plot singletons vs coverage
rule SingCovPlot:
	output:
		os.path.join(BASE_OUT,config.get("rules").get("SingCovPlot").get("out_dir"), "{vcf_name}_SingCov.pdf"),
		os.path.join(BASE_OUT,config.get("rules").get("SingCovPlot").get("out_dir"), "{vcf_name}_SingCov.txt")
	input:
		sample_coverage=rules.sampleDP.output[0],
		sample_singletons=rules.singletons.output[0]
	params:
	log:
		config["paths"]["log_dir"] + "/{vcf_name}-SingCovPlot.log",
		config["paths"]["log_dir"] + "/{vcf_name}-SingCovPlot.e"
	threads: 1
	resources:
		mem_mb=5000
	benchmark:
		config["paths"]["benchmark"] + "/{vcf_name}_SingCovPlot.tsv"
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
			plot_sing_vs_cov(input.sample_singletons, input.sample_coverage, output[0],output[1])
			logger.info('Ended!')
		except Exception as e: 
			logger.error(e, exc_info=True)

#plot het rate per sample by average coverage
rule PlotHetRateSampleCov:
	output:
		os.path.join(BASE_OUT,config.get("rules").get("PlotHetRateSampleCov").get("out_dir"), "{vcf_name}_hetRateByCov.pdf")
	input:
		rules.SampleGetHetRateOut.output[0],
		rules.sampleDP.output[0]
	params:
		sex_table=config.get('paths').get('sex_table')
	log:
		config["paths"]["log_dir"] + "/{vcf_name}-PlotHetRateSampleCov.log",
		config["paths"]["log_dir"] + "/{vcf_name}-PlotHetRateSampleCov.e"
	threads: 1
	resources:
		mem_mb=5000
	benchmark:
		config["paths"]["benchmark"] + "/{vcf_name}_PlotHetRateSampleCov.tsv"
	run:
		logger = logging.getLogger('logging_test')
		fh = logging.FileHandler(str(log[1]))
		fh.setLevel(logging.DEBUG)
		formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
		fh.setFormatter(formatter)
		logger.addHandler(fh)
		try: 
			logger.info('Starting operation!')
			# do something
			plot_het_rate_vs_coverage(input[0], input[1],params.sex_table, output[0])
			logger.info('Ended!')
		except Exception as e: 
			logger.error(e, exc_info=True)
