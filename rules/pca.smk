#set of luse to generate pca plots to identify outliers
#We can reuse some of the rules implementerd by margherita.francescatto in the SNP array qc pipeline.

rule cleanPlinkFilesPCA:
	output:
		vcf=os.path.join(BASE_OUT,config.get("rules").get("kingPCA").get("out_dir"), "{vcf_name}_cleaned.LD0.3.vcf.gz"),
		bed=os.path.join(BASE_OUT,config.get("rules").get("kingPCA").get("out_dir"), "{vcf_name}_cleaned.LD0.3.bed"),
		bim=os.path.join(BASE_OUT,config.get("rules").get("kingPCA").get("out_dir"), "{vcf_name}_cleaned.LD0.3.bim"),
		fam=os.path.join(BASE_OUT,config.get("rules").get("kingPCA").get("out_dir"), "{vcf_name}_cleaned.LD0.3.fam")
	input:
		vcf=rules.cleanMissingHwe.output[0]
	params:
		plink=config.get("PLINK"),
		bcftools=config.get("BCFTOOLS"),
		prefix=os.path.join(BASE_OUT,config.get("rules").get("kingPCA").get("out_dir"), "{vcf_name}_cleaned.LD0.3")
	log:
		config["paths"]["log_dir"] + "/{vcf_name}-cleanPlinkFilesPCA.log",
		config["paths"]["log_dir"] + "/{vcf_name}-cleanPlinkFilesPCA.e"
	threads: 1
	resources:
		mem_mb=5000
	benchmark:
		config["paths"]["benchmark"] + "/{vcf_name}_cleanPlinkFilesPCA.tsv"
	envmodules:
		"plink/1.90",
		"bcftools/1.14"
	shell:
		"""
		{params.bcftools} +prune -m 0.3 -N "rand" -O z -o {output.vcf} {input.vcf} 1> {log[0]} 2> {log[1]}
		{params.plink} --vcf {output.vcf} --keep-allele-order --double-id --biallelic-only 'strict' --snps-only --autosome --make-bed --out {params.prefix} 1>> {log[0]} 2>> {log[1]}
		"""

### Rule from https://gitlab.burlo.trieste.it/marghi/snparrayqc_snakemake/blob/master/Snakefile
### Calculate and plot a PCA projection of the study cohort onto a reference population, 1000Genomes Project phase 3 data, for visual identification of outliers in the study cohort.
rule kingPCA:
	output:
		pc=os.path.join(BASE_OUT,config.get("rules").get("kingPCA").get("out_dir"), "{vcf_name}_cleaned.LD0.3_kingpcapc.txt"),
		proj_pc=os.path.join(BASE_OUT,config.get("rules").get("kingPCA").get("out_dir"), "{vcf_name}_cleaned.LD0.3_kingpcaprojpc.txt"),
		proj_dist=os.path.join(BASE_OUT,config.get("rules").get("kingPCA").get("out_dir"), "{vcf_name}_cleaned.LD0.3_kingpcaproj_Dist.txt"),
		proj_popref=os.path.join(BASE_OUT,config.get("rules").get("kingPCA").get("out_dir"), "{vcf_name}_cleaned.LD0.3_kingpcaproj_popref.txt"),
		plot_pca=os.path.join(BASE_OUT,config.get("rules").get("kingPCA").get("out_dir"), "{vcf_name}_pca.png"),
		plot_pcaproj=os.path.join(BASE_OUT,config.get("rules").get("kingPCA").get("out_dir"), "{vcf_name}_pca_projection_on_1000GP.png")
		# pc="{w_dir}/01_intermediate_files_to_delete/11_pca/{data_name}_cleaned.LD0.3_kingpcapc.txt",
		# proj_pc="{w_dir}/01_intermediate_files_to_delete/11_pca/{data_name}_cleaned.LD0.3_kingpcaprojpc.txt",
		# proj_dist="{w_dir}/01_intermediate_files_to_delete/11_pca/{data_name}_cleaned.LD0.3_kingpcaproj_Dist.txt",
		# proj_popref="{w_dir}/01_intermediate_files_to_delete/11_pca/{data_name}_cleaned.LD0.3_kingpcaproj_popref.txt",
		# plot_pca="{w_dir}/00_release/plots/{data_name}_pca.png",
		# plot_pcaproj="{w_dir}/00_release/plots/{data_name}_pca_projection_on_1000GP.png"
	input:
		ibed=os.path.join(BASE_OUT,config.get("rules").get("kingPCA").get("out_dir"), "{vcf_name}_cleaned.LD0.3.bed")
	params:
		pca_pref=os.path.join(BASE_OUT,config.get("rules").get("kingPCA").get("out_dir"), "{vcf_name}_cleaned.LD0.3_kingpca"),
		proj_pref=os.path.join(BASE_OUT,config.get("rules").get("kingPCA").get("out_dir"), "{vcf_name}_cleaned.LD0.3_kingpcaproj"),
		plot_dir=os.path.join(BASE_OUT,config.get("rules").get("kingPCA").get("out_dir")),
		tgRefBed=config.get("paths").get("1000G_ref_for_king"),
		king=config['KING'],
		scripts=config.get('paths').get('scripts')
	log:
		config["paths"]["log_dir"] + "/{vcf_name}_kingPCA.log",
		config["paths"]["log_dir"] + "/{vcf_name}_kingPCA.e"
	threads: 1
	resources:
		mem_mb=5000
	envmodules:
		"R/4.0.5"
	benchmark:
		config["paths"]["benchmark"] + "/{vcf_name}_kingPCA.tsv"
	shell:
		"""
		{params.king} -b {input.ibed} --mds --prefix {params.pca_pref} 1> {log[0]} 2> {log[1]}
		{params.king} -b {params.tgRefBed},{input.ibed} --projection --mds --prefix {params.proj_pref} 1>> {log[0]} 2>> {log[1]}
		"""
		# Rscript --no-save {params.scripts}/PCA.R {wildcards.vcf_name} {params.plot_dir} {params.pca_pref} {params.proj_pref} 1>> {log[0]} 2>> {log[1]}
# #plot PCA
# rule kingPCAplot:
# 	output:
# 		# pc=os.path.join(BASE_OUT,config.get("rules").get("kingPCA").get("out_dir"), "{vcf_name}_cleaned.LD0.3_kingpcapc.txt"),
# 		# proj_pc=os.path.join(BASE_OUT,config.get("rules").get("kingPCA").get("out_dir"), "{vcf_name}_cleaned.LD0.3_kingpcaprojpc.txt"),
# 		# proj_dist=os.path.join(BASE_OUT,config.get("rules").get("kingPCA").get("out_dir"), "{vcf_name}_cleaned.LD0.3_kingpcaproj_Dist.txt"),
# 		# proj_popref=os.path.join(BASE_OUT,config.get("rules").get("kingPCA").get("out_dir"), "{vcf_name}_cleaned.LD0.3_kingpcaproj_popref.txt"),
# 		plot_pca=os.path.join(BASE_OUT,config.get("rules").get("kingPCA").get("out_dir"), "{vcf_name}_pca.png"),
# 		plot_pcaproj=os.path.join(BASE_OUT,config.get("rules").get("kingPCA").get("out_dir"), "{vcf_name}_pca_projection_on_1000GP.png")
# 	input:
# 		ibed=os.path.join(BASE_OUT,config.get("rules").get("kingPCA").get("out_dir"), "{vcf_name}_cleaned.LD0.3.bed")
# 	params:
# 		pca_pref=os.path.join(BASE_OUT,config.get("rules").get("kingPCA").get("out_dir"), "{vcf_name}_cleaned.LD0.3_kingpca"),
# 		proj_pref=os.path.join(BASE_OUT,config.get("rules").get("kingPCA").get("out_dir"), "{vcf_name}_cleaned.LD0.3_kingpcaproj"),
# 		plot_dir=os.path.join(BASE_OUT,config.get("rules").get("kingPCA").get("out_dir")),
# 		tgRefBed=config.get("paths").get("1000G_ref_for_king"),
# 		scripts=config.get('paths').get('scripts')
# 	log:
# 		config["paths"]["log_dir"] + "/{vcf_name}_kingPCA.log",
# 		config["paths"]["log_dir"] + "/{vcf_name}_kingPCA.e"
# 	threads: 1
# 	resources:
# 		mem_mb=5000
# 	envmodules:
# 		"R/4.0.5"
# 	benchmark:
# 		config["paths"]["benchmark"] + "/{vcf_name}_kingPCA.tsv"
# 	run:
# 		"""
# 		{params.king} -b {input.ibed} --mds --prefix {params.pca_pref} 1> {log[0]} 2> {log[1]}
# 		{params.king} -b {params.tgRefBed},{input.ibed} --projection --mds --prefix {params.proj_pref} 1>> {log[0]} 2>> {log[1]}
# 		Rscript --no-save {params.scripts}/PCA.R {wildcards.vcf_name} {params.plot_dir} {params.pca_pref} {params.proj_pref} 1>> {log[0]} 2>> {log[1]}
# 		"""
