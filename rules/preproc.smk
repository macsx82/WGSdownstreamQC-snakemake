#module containing preprocessing rules for the last downstream QC on WGS data

#Rule to generate a sample list to be used in other rules
rule getSamples:
	output:
		os.path.join(BASE_OUT,config.get("rules").get("getSamples").get("out_dir"), PROJECT_NAME + "_SamplesList.txt")
	input:
		MAIN_VCF_INPUT
	params:
		bcftools=config.get("BCFTOOLS")
	log:
		config["paths"]["log_dir"] + "/getSamples.log",
		config["paths"]["log_dir"] + "/getSamples.e"
	benchmark:
		config["paths"]["benchmark"] + "/getSamples.tsv"
	envmodules:
		"bcftools/1.14"
	resources:
		mem_mb=5000
	message: """ Get samples list """
	shell:
		"""
		{params.bcftools} query -l {input} > {output[0]}
		"""

    
#Rule to apply a more stringen VQSLOD filter, if needed, to snps
rule reapplyVQSRsnps:
	output:
		os.path.join(BASE_OUT,config.get("rules").get("reapplyVQSRsnps").get("out_dir"), PROJECT_NAME + "_snps_VQSLODrefilter.vcf.gz"),
		os.path.join(BASE_OUT,config.get("rules").get("reapplyVQSRsnps").get("out_dir"), PROJECT_NAME + "_snps_VQSLODrefilter.vcf.gz.tbi")
	input:
		MAIN_VCF_INPUT
	params:
		vqslod_thr=config.get("rules").get("reapplyVQSRsnps").get("VQSLOD_thr"),
		bcftools_bin=config.get("BCFTOOLS")
	log:
		config["paths"]["log_dir"] + "/reapplyVQSRsnps.log",
		config["paths"]["log_dir"] + "/reapplyVQSRsnps.e"
	benchmark:
		config["paths"]["benchmark"] + "/reapplyVQSRsnps.tsv"
	envmodules:
		"bcftools/1.14"
	resources:
		mem_mb=5000
	message: """ Filter SNPs """
	shell:
		"""
		{params.bcftools_bin} view -i 'VQSLOD >= {params.vqslod_thr}' -v snps {input} -O z -o {output[0]}
		{params.bcftools_bin} index -t {output[0]}
		"""

#Rule to apply a more stringen VQSLOD filter, if needed, to indels
rule reapplyVQSRindels:
	output:
		os.path.join(BASE_OUT,config.get("rules").get("reapplyVQSRindels").get("out_dir"), PROJECT_NAME + "_indels_VQSLODrefilter.vcf.gz"),
		os.path.join(BASE_OUT,config.get("rules").get("reapplyVQSRindels").get("out_dir"), PROJECT_NAME + "_indels_VQSLODrefilter.vcf.gz.tbi")
	input:
		MAIN_VCF_INPUT
	params:
		vqslod_thr=config.get("rules").get("reapplyVQSRindels").get("VQSLOD_thr"),
		bcftools_bin=config.get("BCFTOOLS")
	log:
		config["paths"]["log_dir"] + "/reapplyVQSRindels.log",
		config["paths"]["log_dir"] + "/reapplyVQSRindels.e"
	benchmark:
		config["paths"]["benchmark"] + "/reapplyVQSRindels.tsv"
	envmodules:
		"bcftools/1.14"
	resources:
		mem_mb=5000
	message: """ Filter INDELs """
	shell:
		"""
		{params.bcftools_bin} view -i 'VQSLOD >= {params.vqslod_thr}' -v indels {input} -O z -o {output[0]}
		{params.bcftools_bin} index -t {output[0]}
		"""

#Rule to concat again the refiltered data
rule mergeReapplyVQSR:
	output:
		# temp(os.path.join(BASE_OUT,config.get("rules").get("mergeReapplyVQSR").get("out_dir"), "{vcf_name}_concatVQSLODrefilter.vcf.gz")),
		# os.path.join(BASE_OUT,config.get("rules").get("mergeReapplyVQSR").get("out_dir"), "{vcf_name}_concatVQSLODrefilter.vcf.gz"),
		os.path.join(BASE_OUT,config.get("rules").get("mergeReapplyVQSR").get("out_dir"), "{vcf_name}_VQSLODrefilter.vcf.gz"),
		os.path.join(BASE_OUT,config.get("rules").get("mergeReapplyVQSR").get("out_dir"), "{vcf_name}_VQSLODrefilter.vcf.gz.tbi")
	input:
		vcf_snps=rules.reapplyVQSRsnps.output[0],
		vcf_indels=rules.reapplyVQSRindels.output[0]
	params:
		bcftools=config.get("BCFTOOLS"),
		tmp=config.get("paths").get("tmp"),
		ref_genome=resolve_single_filepath(*references_abs_path(), config.get("genome_fasta"))
	log:
		config["paths"]["log_dir"] + "/{vcf_name}-mergeReapplyVQSR.log",
		config["paths"]["log_dir"] + "/{vcf_name}-mergeReapplyVQSR.e"
	benchmark:
		config["paths"]["benchmark"] + "/{vcf_name}_mergeReapplyVQSR.tsv"
	envmodules:
		"bcftools/1.14"
	resources:
		mem_mb=10000
	threads: 3
	message: """ Merge back refiltered data """
	shell:
		"""
		temp=$(mktemp -u -d -p {params.tmp})
		({params.bcftools} concat -a {input.vcf_snps} {input.vcf_indels} | {params.bcftools} sort -T ${{temp}} | {params.bcftools} norm -f {params.ref_genome} -O z -o {output[0]}) > {log[0]} 2> {log[1]}
        {params.bcftools} index -t {output[0]}
		"""
