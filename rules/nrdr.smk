#this set of rules is separated, since we will need also to get some additional data, like genotypes for all our existing samples

#1) convert to vcf the SNP array data
rule Plink2Vcf:
	output:
		os.path.join(BASE_OUT, config.get('rules').get('Plink2Vcf').get('out_dir'), "{vcf_name}_array_data.vcf.gz"),
		os.path.join(BASE_OUT, config.get('rules').get('Plink2Vcf').get('out_dir'), "{vcf_name}_array_data.vcf.gz.tbi"),
	input:
		config.get('paths').get('snp_array_data')
	params:
		plink=config['PLINK']
	log:
		config["paths"]["log_dir"] + "/{vcf_name}-Plink2Vcf.log",
		config["paths"]["log_dir"] + "/{vcf_name}-Plink2Vcf.e"
	threads: 1
	resources:
		mem_mb=5000
	benchmark:
		config["paths"]["benchmark"] + "/{vcf_name}_Plink2Vcf.tsv"
	envmodules:
		"plink/1.90"
	shell:
		"""
		{params.plink} --file {input[0]} --recode 
		"""

#2) generate concordance stats using bcftools

rule NRD:
	output:
		os.path.join(BASE_OUT, config.get('rules').get('Plink2Vcf').get('out_dir'), "{vcf_name}_NRDR.txt"),
	input:
		snp_array=config.get('paths').get('snp_array_data'),
		vcf=rules.cleanMissingHwe.output[0],
		vcf_index=rules.cleanMissingHwe.output[1]
