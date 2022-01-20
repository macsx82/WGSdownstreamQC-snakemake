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
		{params.bcftools} norm -m+both {input[0]} | bcftools view -m2 -M2 -v snps -O z -o {output[0]}
		{params.bcftools} index -t {output[0]}
		"""

#2) generate concordance stats using bcftools
# rule NRD:
# 	output:
# 		os.path.join(BASE_OUT, config.get('rules').get('Plink2Vcf').get('out_dir'), "{vcf_name}_NRDR.txt"),
# 	input:
# 		snp_array=config.get('paths').get('snp_array_data'),
# 		vcf=rules.cleanMissingHwe.output[0],
# 		vcf_index=rules.cleanMissingHwe.output[1]
