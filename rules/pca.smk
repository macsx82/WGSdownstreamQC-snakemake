#set of luse to generate pca plots to identify outliers
#We can reuse some of the rules implementerd by margherita.francescatto in the SNP array qc pipeline.

rule cleanPlinkFilesPCA:
	output:
		pc="{w_dir}/01_intermediate_files_to_delete/11_pca/{data_name}_cleaned.LD0.3_kingpcapc.txt",
		proj_pc="{w_dir}/01_intermediate_files_to_delete/11_pca/{data_name}_cleaned.LD0.3_kingpcaprojpc.txt",
		proj_dist="{w_dir}/01_intermediate_files_to_delete/11_pca/{data_name}_cleaned.LD0.3_kingpcaproj_Dist.txt",
		proj_popref="{w_dir}/01_intermediate_files_to_delete/11_pca/{data_name}_cleaned.LD0.3_kingpcaproj_popref.txt",
		plot_pca="{w_dir}/00_release/plots/{data_name}_pca.png",
		plot_pcaproj="{w_dir}/00_release/plots/{data_name}_pca_projection_on_1000GP.png"
	input:
		vcf=rules.cleanMissingHwe.output[0]
	params:
		plink=config.get("PLINK"),
	log:
		stde="{w_dir}/00_release/logs/20_{data_name}_kingPCA.stderr",
		stdo="{w_dir}/00_release/logs/20_{data_name}_kingPCA.stdout",
		plot_stde="{w_dir}/00_release/logs/rplot_{data_name}_kingPCA.stderr",
		plot_stdo="{w_dir}/00_release/logs/rplot_{data_name}_kingPCA.stdout"
	threads:
	resources:
	benchmark:
	shell:
		"""
		{params.plink}
		"""

### Rule from https://gitlab.burlo.trieste.it/marghi/snparrayqc_snakemake/blob/master/Snakefile
### Calculate and plot a PCA projection of the study cohort onto a reference population, 1000Genomes Project phase 3 data, for visual identification of outliers in the study cohort.
rule kingPCA:
	input:
		ibed="{w_dir}/00_release/output/{data_name}_cleaned.LD0.3.bed"
	output:
		pc="{w_dir}/01_intermediate_files_to_delete/11_pca/{data_name}_cleaned.LD0.3_kingpcapc.txt",
		proj_pc="{w_dir}/01_intermediate_files_to_delete/11_pca/{data_name}_cleaned.LD0.3_kingpcaprojpc.txt",
		proj_dist="{w_dir}/01_intermediate_files_to_delete/11_pca/{data_name}_cleaned.LD0.3_kingpcaproj_Dist.txt",
		proj_popref="{w_dir}/01_intermediate_files_to_delete/11_pca/{data_name}_cleaned.LD0.3_kingpcaproj_popref.txt",
		plot_pca="{w_dir}/00_release/plots/{data_name}_pca.png",
		plot_pcaproj="{w_dir}/00_release/plots/{data_name}_pca_projection_on_1000GP.png"
	params:
		pca_pref="{w_dir}/01_intermediate_files_to_delete/11_pca/{data_name}_cleaned.LD0.3_kingpca",
		proj_pref="{w_dir}/01_intermediate_files_to_delete/11_pca/{data_name}_cleaned.LD0.3_kingpcaproj",
		plot_dir="{w_dir}/00_release/plots/",
		tgRefBed=config.get("paths").get("1000G_ref_for_king")
	log:
		stde="{w_dir}/00_release/logs/20_{data_name}_kingPCA.stderr",
		stdo="{w_dir}/00_release/logs/20_{data_name}_kingPCA.stdout",
		plot_stde="{w_dir}/00_release/logs/rplot_{data_name}_kingPCA.stderr",
		plot_stdo="{w_dir}/00_release/logs/rplot_{data_name}_kingPCA.stdout"
	threads:
	resources:
	benchmark:
	shell:
		"""
		{king_path}/king -b {input.ibed} --mds --prefix {params.pca_pref} > {log.stdo} 2> {log.stde}
		{king_path}/king -b {params.tgRefBed},{input.ibed} --projection --mds --prefix {params.proj_pref} >> {log.stdo} 2>> {log.stde}
		Rscript --no-save scripts/PCA.R {data_name} {params.plot_dir} {params.pca_pref} {params.proj_pref} > {log.plot_stdo} 2> {log.plot_stde}
		"""
