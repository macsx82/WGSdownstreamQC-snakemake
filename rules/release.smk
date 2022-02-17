#rule to generate a single folder with all data needed for the filtering step
rule release:
	output:
		os.path.join(BASE_OUT,config.get("rules").get("release").get("out_dir"), "{out_name}_cleaned.LD0.3_kingpcaprojpc.txt"),
		os.path.join(BASE_OUT,config.get("rules").get("release").get("out_dir"), "{out_name}_hetRate.txt"),
		os.path.join(BASE_OUT,config.get("rules").get("release").get("out_dir"), "{out_name}_missing_ALL.imiss"),
		os.path.join(BASE_OUT,config.get("rules").get("release").get("out_dir"), "{out_name}_SingCov.txt"),
		expand(os.path.join(BASE_OUT,config.get("rules").get("release").get("out_dir"), "{vcf_name}_{ref_pop}_af_extrDiff.txt"), vcf_name=out_prefix,ref_pop=ref_pop)
	input:
		rules.kingPCA.output.proj_pc,
		rules.SampleGetHetRateOut.output[0],
		rules.collectSampleMissingRate.output[0],
		rules.SingCovPlot.output[1],
		#variants
		expand(os.path.join(BASE_OUT,config.get("rules").get("comparePopAF").get("out_dir"), "{vcf_name}_{ref_pop}_af_extrDiff.txt"),vcf_name=out_prefix,ref_pop=ref_pop)
	params:
		base_out=os.path.join(BASE_OUT,config.get("rules").get("release").get("out_dir"))
	log:
		config["paths"]["log_dir"] + "/Release_filter_files.log",
		config["paths"]["log_dir"] + "/Release_filter_files.e"
	threads: 1
	resources:
		mem_mb=5000
	benchmark:
		config["paths"]["benchmark"] + "/Release_filter_files.tsv"
	run:
		for c_out in input:
			#get file basename
			outname=os.path.basename(c_out)
			outpath=os.path.join(params.base_out, outname)
			cp_cmd="cp %s %s" %(c_out, outpath)


if SNP_DATA != "NONE" :
	#rule to collect the NRD output files to move in the release folder. This rule will be executed only if SNP array data is provided
	rule release_nrd:
		output:
			os.path.join(BASE_OUT, config.get('rules').get('release').get('out_dir'), "{out_name}_NRDRsamples.txt"),
			os.path.join(BASE_OUT, config.get('rules').get('release').get('out_dir'), "{out_name}_NRDRsites.txt"),
			expand(os.path.join(BASE_OUT,config.get("rules").get("release").get("out_dir"), "{vcf_name}_ARRAY_af_extrDiff.txt"), vcf_name=out_prefix_autosomal)
		input:
			rules.NRDbySample.output[0],
			#variants
			rules.NRDbySite.output[0],
			expand(os.path.join(BASE_OUT,config.get("rules").get("getArrayPopAF").get("out_dir"), "{vcf_name}_ARRAY_af_extrDiff.txt"), vcf_name=out_prefix_autosomal)
		params:
			base_out=os.path.join(BASE_OUT,config.get("rules").get("release").get("out_dir"))
		log:
			config["paths"]["log_dir"] + "/release_nrd_files.log",
			config["paths"]["log_dir"] + "/release_nrd_files.e"
		threads: 1
		resources:
			mem_mb=5000
		benchmark:
			config["paths"]["benchmark"] + "/release_nrd_files.tsv"
		run:
			for c_out in input:
				#get file basename
				outname=os.path.basename(c_out)
				outpath=os.path.join(params.base_out, outname)
				cp_cmd="cp %s %s" %(c_out, outpath)

