#this set of rules is separated, since we will need also to get some additional data, like genotypes for all our existing samples

#1) Remove multiallelic sites from the VCF
rule VcfMultiClean:
	output:
		os.path.join(BASE_OUT, config.get('rules').get('VcfMultiClean').get('out_dir'), "{vcf_name}_SNPSMultiClean.vcf.gz"),
		os.path.join(BASE_OUT, config.get('rules').get('VcfMultiClean').get('out_dir'), "{vcf_name}_SNPSMultiClean.vcf.gz.tbi"),
	input:
		rules.cleanMissingHwe.output[0],
		rules.cleanMissingHwe.output[1],
	params:
		bcftools=config['BCFTOOLS']
	log:
		config["paths"]["log_dir"] + "/{vcf_name}-VcfMultiClean.log",
		config["paths"]["log_dir"] + "/{vcf_name}-VcfMultiClean.e"
	threads: 2
	resources:
		mem_mb=5000
	benchmark:
		config["paths"]["benchmark"] + "/{vcf_name}_VcfMultiClean.tsv"
	envmodules:
		"bcftools/1.14"
	shell:
		"""
		({params.bcftools} norm -m+both {input[0]} | bcftools view -m2 -M2 -v snps -O z -o {output[0]}) 1> {log[0]} 2> {log[1]}
		{params.bcftools} index -t {output[0]}
		"""

#get the common samples from the array data
rule VcfArrayCommonSamples:
	output:
		os.path.join(BASE_OUT, config.get('rules').get('VcfMultiClean').get('out_dir'), "{vcf_name}_VcfArrayCommonSamples.vcf.gz"),
		os.path.join(BASE_OUT, config.get('rules').get('VcfMultiClean').get('out_dir'), "{vcf_name}_VcfArrayCommonSamples.vcf.gz.tbi"),
		temp(os.path.join(BASE_OUT, config.get('rules').get('VcfMultiClean').get('out_dir'), "{vcf_name}_VcfArraySamples.txt"))
	input:
		snp_array=config.get('paths').get('snp_array_data'),
		wgs_samples=rules.getSamples.output[0]
		# rules.VcfMultiClean.output[0],
	params:
		bcftools=config['BCFTOOLS']
	log:
		config["paths"]["log_dir"] + "/{vcf_name}-VcfArrayCommonSamples.log",
		config["paths"]["log_dir"] + "/{vcf_name}-VcfArrayCommonSamples.e"
	threads: 1
	resources:
		mem_mb=5000
	benchmark:
		config["paths"]["benchmark"] + "/{vcf_name}_VcfArrayCommonSamples.tsv"
	envmodules:
		"bcftools/1.14"
	shell:
		"""
		({params.bcftools} view -S {input.wgs_samples} --force-samples -O z -o {output[0]} {input.snp_array}) 1> {log[0]} 2> {log[1]}
		{params.bcftools} index -t {output[0]}
		{params.bcftools} query -l {output[0]} > {output[2]}
		"""

#get the common samples from the wgs data
rule VcfWgsArrayCommon:
	output:
		os.path.join(BASE_OUT, config.get('rules').get('VcfMultiClean').get('out_dir'), "{vcf_name}_VcfWgsArrayCommon.vcf.gz"),
		os.path.join(BASE_OUT, config.get('rules').get('VcfMultiClean').get('out_dir'), "{vcf_name}_VcfWgsArrayCommon.vcf.gz.tbi"),
	input:
		vcf_file=rules.VcfMultiClean.output[0],
		snp_array=rules.VcfArrayCommonSamples.output[0],
		samples_list=rules.VcfArrayCommonSamples.output[2]
	params:
		bcftools=config['BCFTOOLS']
	log:
		config["paths"]["log_dir"] + "/{vcf_name}-VcfWgsArrayCommon.log",
		config["paths"]["log_dir"] + "/{vcf_name}-VcfWgsArrayCommon.e"
	threads: 1
	resources:
		mem_mb=5000
	benchmark:
		config["paths"]["benchmark"] + "/{vcf_name}_VcfWgsArrayCommon.tsv"
	envmodules:
		"bcftools/1.14"
	shell:
		"""
		({params.bcftools} view -R {input.snp_array} -S {input.samples_list} --force-samples -O z -o {output[0]} {input.vcf_file}) 1> {log[0]} 2> {log[1]}
		{params.bcftools} index -t {output[0]}
		"""

#2) generate concordance stats (among others) using bcftools by chromosome
rule NRDstats:
	output:
		# os.path.join(BASE_OUT, config.get('rules').get('NRD').get('out_dir'), "{vcf_name}_{chr}_NRDR.txt"),
		os.path.join(BASE_OUT, config.get('rules').get('NRD').get('out_dir'), "{vcf_name}_NRDR.txt"),
	input:
		snp_array=rules.VcfArrayCommonSamples.output[0],
		vcf=rules.VcfWgsArrayCommon.output[0],
		vcf_index=rules.VcfWgsArrayCommon.output[1],
		samples=rules.VcfArrayCommonSamples.output[2]
	params:
		bcftools=config['BCFTOOLS']
	log:
		config["paths"]["log_dir"] + "/{vcf_name}-NRD.log",
		config["paths"]["log_dir"] + "/{vcf_name}-NRD.e"
	threads: 1
	resources:
		mem_mb=10000
	benchmark:
		config["paths"]["benchmark"] + "/{vcf_name}_NRD.tsv"
	envmodules:
		"bcftools/1.14"
	shell:
		"""
		{params.bcftools} stats {input.snp_array} {input.vcf} -S {input.samples} --verbose > {output[0]} 2> {log[1]}
		"""

#3) get only info for NRD by sites and samples, for each chromosome
rule getNRDbySiteAndSamples:
	output:
		os.path.join(BASE_OUT, config.get('rules').get('NRD').get('out_dir'), "{vcf_name}_NRDRsites.txt"),
		os.path.join(BASE_OUT, config.get('rules').get('NRD').get('out_dir'), "{vcf_name}_NRDRsamples.txt")
	input:
		rules.NRDstats.output[0]
	params:
		bcftools=config['BCFTOOLS']
	log:
		config["paths"]["log_dir"] + "/{vcf_name}-getNRDbySiteAnSamples.log",
		config["paths"]["log_dir"] + "/{vcf_name}-getNRDbySiteAnSamples.e"
	threads: 1
	resources:
		mem_mb=10000
	benchmark:
		config["paths"]["benchmark"] + "/{vcf_name}_getNRDbySiteAnSamples.tsv"
	envmodules:
		"bcftools/1.14"
	shell:
		"""
		egrep "^PSD" {input[0]} > {output[0]}
		egrep "^GCsS" {input[0]} > {output[1]}
		"""

#3) Process the splitted information to calculate the complete tables per sample and per site
# this is an AGGREGATION rule, since we need a genome wide value for concordance and non reference discordance
rule NRDbySample:
	output:
		os.path.join(BASE_OUT, config.get('rules').get('NRD').get('out_dir'), "{out_name}_NRDRsamples.txt"),
	input:
		# all_samples_NRD=expand(os.path.join(BASE_OUT, config.get('rules').get('NRD').get('out_dir'), "{{vcf_name}}_{chr}_NRDRsamples.txt"), chr=autosomal)
		all_samples_NRD=expand(os.path.join(BASE_OUT, config.get('rules').get('NRD').get('out_dir'), "{vcf_name}_NRDRsamples.txt"), vcf_name=out_prefix_autosomal)
	params:
		bcftools=config['BCFTOOLS']
	log:
		config["paths"]["log_dir"] + "/{out_name}-NRDbySample.log",
		config["paths"]["log_dir"] + "/{out_name}-NRDbySample.e"
	threads: 1
	resources:
		mem_mb=5000
	benchmark:
		config["paths"]["benchmark"] + "/{out_name}_NRDbySample.tsv"
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
			get_NRD_by_sample(input.all_samples_NRD, output[0])
			logger.info('Ended!')
		except Exception as e: 
			logger.error(e, exc_info=True)

#simple AGGREGATION rule to concat all chr NRD sites values and add a fla
rule NRDbySite:
	output:
		os.path.join(BASE_OUT, config.get('rules').get('NRD').get('out_dir'), "{out_name}_NRDRsites.txt"),
	input:
		all_sites_NRD=expand(os.path.join(BASE_OUT, config.get('rules').get('NRD').get('out_dir'), "{vcf_name}_NRDRsites.txt"), vcf_name=out_prefix_autosomal)
	log:
		config["paths"]["log_dir"] + "/{out_name}-NRDbySite.log",
		config["paths"]["log_dir"] + "/{out_name}-NRDbySite.e"
	threads: 1
	resources:
		mem_mb=2000
	benchmark:
		config["paths"]["benchmark"] + "/{out_name}_NRDbySite.tsv"
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
			get_NRD_by_site(input.all_sites_NRD, output[0])
			logger.info('Ended!')
		except Exception as e: 
			logger.error(e, exc_info=True)
