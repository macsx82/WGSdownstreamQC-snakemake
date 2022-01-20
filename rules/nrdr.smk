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

#2) generate concordance stats (among others) using bcftools
rule NRDstats:
	output:
		os.path.join(BASE_OUT, config.get('rules').get('NRDR').get('out_dir'), "{vcf_name}_{chr}_NRDR.txt"),
	input:
		snp_array=config.get('paths').get('snp_array_data'),
		vcf=rules.VcfMultiClean.output[0],
		vcf_index=rules.VcfMultiClean.output[1],
		samples=rules.getSamples.output[0]
	params:
		bcftools=config['BCFTOOLS']
	log:
		config["paths"]["log_dir"] + "/{vcf_name}-{chr}-NRD.log",
		config["paths"]["log_dir"] + "/{vcf_name}-{chr}-NRD.e"
	threads: 1
	resources:
		mem_mb=10000
	benchmark:
		config["paths"]["benchmark"] + "/{vcf_name}_{chr}_NRD.tsv"
	envmodules:
		"bcftools/1.14"
	shell:
		"""
		{params.bcftools} stats -r {wildcards.chr} {input.snp_array} {input.vcf} -S {input.samples} --verbose > {output[0]} 2> {log[1]}
		"""

#3) get only info for NRD by sites
rule getNRDbySiteAndSamples:
	output:
		os.path.join(BASE_OUT, config.get('rules').get('NRDR').get('out_dir'), "{vcf_name}_{chr}_NRDRsites.txt"),
		os.path.join(BASE_OUT, config.get('rules').get('NRDR').get('out_dir'), "{vcf_name}_{chr}_NRDRsamples.txt")
	input:
		rules.NRDstats.output[0]
	params:
		bcftools=config['BCFTOOLS']
	log:
		config["paths"]["log_dir"] + "/{vcf_name}-{chr}-getNRDbySiteAnSamples.log",
		config["paths"]["log_dir"] + "/{vcf_name}-{chr}-getNRDbySiteAnSamples.e"
	threads: 1
	resources:
		mem_mb=10000
	benchmark:
		config["paths"]["benchmark"] + "/{vcf_name}_{chr}_getNRDbySiteAnSamples.tsv"
	envmodules:
		"bcftools/1.14"
	shell:
		"""
		egrep "^GCsS" {input[0]} > {output[0]}
		egrep "^PSD" {input[0]} > {output[1]}
	
		"""

#3) Process the splitted information to calculate the complete tables per sample and per site
# rule getNRDbySample:
# 	output:
# 		os.path.join(BASE_OUT, config.get('rules').get('NRD').get('out_dir'), "{vcf_name}_{chr}_NRDRsamples.txt"),
# 	input:
# 		rules.NRDstats.output[0]
# 	params:
# 		bcftools=config['BCFTOOLS']
# 	log:
# 		config["paths"]["log_dir"] + "/{vcf_name}-{chr}-getNRDbySample.log",
# 		config["paths"]["log_dir"] + "/{vcf_name}-{chr}-getNRDbySample.e"
# 	threads: 1
# 	resources:
# 		mem_mb=10000
# 	benchmark:
# 		config["paths"]["benchmark"] + "/{vcf_name}_{chr}_getNRDbySample.tsv"
# 	envmodules:
# 		"bcftools/1.14"
# 	shell:
# 		"""
# 		{params.bcftools} stats {input.snp_array} {input.vcf} -S {input.samples} --verbose > {output[0]} 2> {log[1]}
# 		"""
