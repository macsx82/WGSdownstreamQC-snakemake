#module containing preprocessing rules for the last downstream QC on WGS data

#Rule to apply a more stringen VQSLOD filter, if needed, to snps
rule reapplyVQSRsnps:
	output:
		os.path.join(BASE_OUT,config.get("rules").get("reapplyVQSRsnps").get("out_dir"), PROJECT_NAME + "_snps_VQSLODrefilter.vcf.gz")
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
		{params.bcftools_bin} view -i 'VQSLOD >= {params.vqslod_thr}' -v snps {input} -O z -o {output}

		"""

#Rule to apply a more stringen VQSLOD filter, if needed, to indels
rule reapplyVQSRindels:
	output:
		os.path.join(BASE_OUT,config.get("rules").get("reapplyVQSRindels").get("out_dir"), PROJECT_NAME + "_indels_VQSLODrefilter.vcf.gz")
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
		{params.bcftools_bin} view -i 'VQSLOD >= {params.vqslod_thr}' -v indels {input} -O z -o {output}
		"""

#Rule to concat again the refiltered data
rule mergeReapplyVQSR:
	output:
		os.path.join(BASE_OUT,config.get("rules").get("mergeReapplyVQSR").get("out_dir"), PROJECT_NAME + "_MERGED_VQSLODrefilter.vcf.gz"),
		os.path.join(BASE_OUT,config.get("rules").get("mergeReapplyVQSR").get("out_dir"), PROJECT_NAME + "_MERGED_VQSLODrefilter.vcf.gz.tbi")
	input:
		vcf_snps=rules.reapplyVQSRsnps.output,
		vcf_indels=rules.reapplyVQSRindels.output
	params:
		bcftools=config.get("BCFTOOLS"),
		tmp=config.get("paths").get("tmp"),
		ref_genome=resolve_single_filepath(*references_abs_path(), config.get("genome_fasta"))
	log:
		config["paths"]["log_dir"] + "/mergeReapplyVQSR.log",
		config["paths"]["log_dir"] + "/mergeReapplyVQSR.e"
	benchmark:
		config["paths"]["benchmark"] + "/mergeReapplyVQSR.tsv"
	envmodules:
		"bcftools/1.14"
	resources:
		mem_mb=10000
	message: """ Merge back refiltered data """
	shell:
		"""
		temp=$(mktemp -u -d -p {params.tmp})
		{params.bcftools} concat {input.vcf_snps} {input.vcf_indels}| {params.bcftools} sort -T ${{temp}} | {params.bcftools} norm -f {params.ref_genome} -O z -o {output[0]} > {log[0]} 2> {log[1]}
        {params.bcftools} index -t {output[0]}
		"""



