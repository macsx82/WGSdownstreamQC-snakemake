#Keep all the explicit plotting rules here, to keep things organized
#Plot het rate per sample
rule PlotHetRateSample:
	output:
		os.path.join(BASE_OUT,config.get("rules").get("SampleHetRate").get("out_dir"), "{vcf_name}_hetRate.pdf"),
		os.path.join(BASE_OUT,config.get("rules").get("SampleHetRate").get("out_dir"), "{vcf_name}_hetRate_SEQ.pdf"),
	input:
		rules.SampleGetHetRateOut.output[0]
	params:
		vcftools=config['VCFTOOLS'],
		outplot_prefix=os.path.join(BASE_OUT,config.get("rules").get("SampleHetRate").get("out_dir"), "{vcf_name}_hetRate"),
		manifest_table=config.get('paths').get('manifest_table')
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
			plot_het_rate_sample(input[0], params.manifest_table,params.outplot_prefix)
			logger.info('Ended!')
		except Exception as e: 
			logger.error(e, exc_info=True)

# #Plot het rate per sample
# rule PlotHetRateSampleChr:
# 	output:
# 		os.path.join(BASE_OUT,config.get("rules").get("SampleHetRate").get("out_dir"), "{vcf_name}_{chr}_hetRate.pdf"),
# 		os.path.join(BASE_OUT,config.get("rules").get("SampleHetRate").get("out_dir"), "{vcf_name}_{chr}_hetRate_SEQ.pdf"),
# 	input:
# 		# os.path.join(BASE_OUT,config.get("rules").get("SampleHetRate").get("out_dir"), "{vcf_name}_{chr}_hetRate.txt")
# 		rules.SampleGetHetRateOut.output[0]
# 	params:
# 		vcftools=config['VCFTOOLS'],
# 		outplot_prefix=os.path.join(BASE_OUT,config.get("rules").get("SampleHetRate").get("out_dir"), "{vcf_name}_{chr}_hetRate")
# 	log:
# 		config["paths"]["log_dir"] + "/{vcf_name}-{chr}-PlotHetRateSample.log",
# 		config["paths"]["log_dir"] + "/{vcf_name}-{chr}-PlotHetRateSample.e"
# 	threads: 1
# 	resources:
# 		mem_mb=5000
# 	benchmark:
# 		config["paths"]["benchmark"] + "/{vcf_name}_{chr}_PlotHetRateSample.tsv"
# 	run:
# 		logger = logging.getLogger('logging_test')
# 		fh = logging.FileHandler(str(log[1]))
# 		fh.setLevel(logging.DEBUG)
# 		formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
# 		fh.setFormatter(formatter)
# 		logger.addHandler(fh)
# 		try: 
# 			logger.info('Starting operation!')
# 			# do something
# 			plot_het_rate_sample(input[0], params.outplot_prefix)
# 			logger.info('Ended!')
# 		except Exception as e: 
# 			logger.error(e, exc_info=True)


#plot singletons vs coverage
rule SingCovPlot:
	output:
		os.path.join(BASE_OUT,config.get("rules").get("SamplePlots").get("out_dir"), "{vcf_name}_SingCov.pdf"),
		os.path.join(BASE_OUT,config.get("rules").get("SamplePlots").get("out_dir"), "{vcf_name}_SingCov.txt")
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
		os.path.join(BASE_OUT,config.get("rules").get("SamplePlots").get("out_dir"), "{vcf_name}_hetRateByCov_sex.pdf"),
		os.path.join(BASE_OUT,config.get("rules").get("SamplePlots").get("out_dir"), "{vcf_name}_hetRateByCov_cohort.pdf")
	input:
		rules.SampleGetHetRateOut.output[0],
		rules.sampleDP.output[0]
	params:
		manifest_table=config.get('paths').get('manifest_table'),
		plot_prefix=os.path.join(BASE_OUT,config.get("rules").get("SamplePlots").get("out_dir"), "{vcf_name}_hetRateByCov")
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
			plot_het_rate_vs_coverage(input[0], input[1],params.manifest_table, params.plot_prefix)
			logger.info('Ended!')
		except Exception as e: 
			logger.error(e, exc_info=True)

#plot het rate per sample by Missing sites fraction
rule PlotHetRateSampleMissing:
	output:
		os.path.join(BASE_OUT,config.get("rules").get("SamplePlots").get("out_dir"), "{vcf_name}_hetRateByMiss_sex.pdf"),
		os.path.join(BASE_OUT,config.get("rules").get("SamplePlots").get("out_dir"), "{vcf_name}_hetRateByMiss_cohort.pdf")
	input:
		rules.SampleGetHetRateOut.output[0],
		rules.SampleMissingRate.output[0]
	params:
		manifest_table=config.get('paths').get('manifest_table'),
		plot_prefix=os.path.join(BASE_OUT,config.get("rules").get("SamplePlots").get("out_dir"), "{vcf_name}_hetRateByMiss")
	log:
		config["paths"]["log_dir"] + "/{vcf_name}-PlotHetRateSampleMiss.log",
		config["paths"]["log_dir"] + "/{vcf_name}-PlotHetRateSampleMiss.e"
	threads: 1
	resources:
		mem_mb=5000
	benchmark:
		config["paths"]["benchmark"] + "/{vcf_name}_PlotHetRateSampleMiss.tsv"
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
			plot_het_rate_vs_missing(input[0], input[1],params.manifest_table, params.plot_prefix)
			logger.info('Ended!')
		except Exception as e: 
			logger.error(e, exc_info=True)

#plot het rate per sample by singleton count
rule PlotHetRateSampleSing:
	output:
		os.path.join(BASE_OUT,config.get("rules").get("SamplePlots").get("out_dir"), "{vcf_name}_hetRateBySing_sex.pdf"),
		os.path.join(BASE_OUT,config.get("rules").get("SamplePlots").get("out_dir"), "{vcf_name}_hetRateBySing_cohort.pdf")
	input:
		rules.SampleGetHetRateOut.output[0],
		rules.singletons.output[0]
	params:
		manifest_table=config.get('paths').get('manifest_table'),
		plot_prefix=os.path.join(BASE_OUT,config.get("rules").get("SamplePlots").get("out_dir"), "{vcf_name}_hetRateBySing")
	log:
		config["paths"]["log_dir"] + "/{vcf_name}-PlotHetRateSampleSing.log",
		config["paths"]["log_dir"] + "/{vcf_name}-PlotHetRateSampleSing.e"
	threads: 1
	resources:
		mem_mb=5000
	benchmark:
		config["paths"]["benchmark"] + "/{vcf_name}_PlotHetRateSampleSing.tsv"
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
			plot_het_rate_vs_singletons(input[0], input[1],params.manifest_table, params.plot_prefix)
			logger.info('Ended!')
		except Exception as e: 
			logger.error(e, exc_info=True)



#plot het rate per sample by singleton count
rule PlotHetRateSampleNRD:
	output:
		os.path.join(BASE_OUT,config.get("rules").get("SamplePlots").get("out_dir"), "{vcf_name}_hetRateByNRD_sex.pdf"),
		os.path.join(BASE_OUT,config.get("rules").get("SamplePlots").get("out_dir"), "{vcf_name}_hetRateByNRD_cohort.pdf"),
		os.path.join(BASE_OUT,config.get("rules").get("SamplePlots").get("out_dir"), "{vcf_name}_hetRateByNRD_seq.pdf")
	input:
		rules.SampleGetHetRateOut.output[0],
		rules.NRDbySample.output[0]
	params:
		manifest_table=config.get('paths').get('manifest_table'),
		plot_prefix=os.path.join(BASE_OUT,config.get("rules").get("SamplePlots").get("out_dir"), "{vcf_name}_hetRateByNRD")
	log:
		config["paths"]["log_dir"] + "/{vcf_name}-PlotHetRateSampleNRD.log",
		config["paths"]["log_dir"] + "/{vcf_name}-PlotHetRateSampleNRD.e"
	threads: 1
	resources:
		mem_mb=5000
	benchmark:
		config["paths"]["benchmark"] + "/{vcf_name}_PlotHetRateSampleNRD.tsv"
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
			plot_het_rate_vs_nrd(input[0], input[1],params.manifest_table, params.plot_prefix)
			logger.info('Ended!')
		except Exception as e: 
			logger.error(e, exc_info=True)


#plot PCA
rule kingPCAplot:
	output:
		plot_pca=os.path.join(BASE_OUT,config.get("rules").get("kingPCA").get("out_dir"), "{vcf_name}_pca.pdf"),
		plot_pcaproj=os.path.join(BASE_OUT,config.get("rules").get("kingPCA").get("out_dir"), "{vcf_name}_pca_projection_on_1000GP.pdf")
	input:
		pc=os.path.join(BASE_OUT,config.get("rules").get("kingPCA").get("out_dir"), "{vcf_name}_cleaned.LD0.3_kingpcapc.txt"),
		proj_pc=os.path.join(BASE_OUT,config.get("rules").get("kingPCA").get("out_dir"), "{vcf_name}_cleaned.LD0.3_kingpcaprojpc.txt"),
		proj_popref=os.path.join(BASE_OUT,config.get("rules").get("kingPCA").get("out_dir"), "{vcf_name}_cleaned.LD0.3_kingpcaproj_popref.txt")
	params:
		plot_dir=os.path.join(BASE_OUT,config.get("rules").get("kingPCA").get("out_dir")),
		tgRefBed=config.get("paths").get("1000G_ref_for_king"),
		scripts=config.get('paths').get('scripts')
	log:
		config["paths"]["log_dir"] + "/{vcf_name}_kingPCAplot.log",
		config["paths"]["log_dir"] + "/{vcf_name}_kingPCAplot.e"
	threads: 1
	resources:
		mem_mb=5000
	benchmark:
		config["paths"]["benchmark"] + "/{vcf_name}_kingPCAplot.tsv"
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
			PCAplots(params.plot_dir,input.pc,input.proj_pc,input.proj_popref,wildcards.vcf_name)
			logger.info('Ended!')
		except Exception as e: 
			logger.error(e, exc_info=True)
