#Snakefile for the pipeline to perform downstream QC on WGS data
#
#
# 16/12/2021
#
# Author: massimiliano [dot] Cocca [at] burlo [dot] trieste [dot] it
#import libraries
import pandas as pd
import pathlib
import io
import os
import os.path
import re
import logging
from snakemake.exceptions import print_exception, WorkflowError
from snakemake.utils import validate, min_version

##### set minimum snakemake version #####
min_version("5.32.0")

include_prefix="rules"

#define some global variables from the config file
BASE_OUT=config["paths"]["base_out"]
MAIN_VCF_INPUT=config["paths"]["input_vcf"]
chroms=config["chrs"]
autosomal=[ x for x in config["chrs"] if x != "chrX" ]
PROJECT_NAME=config["project_name"]
SNP_DATA=config.get("paths").get("snp_array_data")
# out_prefix=PROJECT_NAME + "_MERGED"
#here we are defining outprefix as a list, since we want to parallelize everything and work by cromosome
#We should be able to switch flawlessly from a chr based to a merged based pipeline,just modifying this parameter
out_prefix=[PROJECT_NAME + "_" + chrom + "_MERGED" for chrom in chroms]
out_prefix_autosomal=[PROJECT_NAME + "_" + chrom + "_MERGED" for chrom in autosomal]

print(SNP_DATA)
### path to resources needed for plots
tgRefBed = config['paths']['1000G_ref_for_king']
ref_pop=list(config.get("rules").get("comparePopAF").get("ref_pops").keys())

##### global wildcard costraints #####
wildcard_constraints:
    vcf_name = "\w+_MERGED"
##### functions #####
include:
    "scripts/functions.py"
include:
    "scripts/PCA.py"
include:
    "scripts/plots.py"
    
##### local rules #####
# localrules: all

##### Target rules #####
if SNP_DATA != "NONE" :
    rule all:
        input:
            #define target input for all first steps needed (basically, initial data filtering to perform other comparisons)
            expand(os.path.join(BASE_OUT,config.get("rules").get("cleanMissingHwe").get("out_dir"), "{vcf_name}_HWE95call.{ext}"),vcf_name=out_prefix, ext=["vcf.gz", "vcf.gz.tbi", "removed.sites.vcf.gz"]),
            expand(os.path.join(BASE_OUT, config.get('rules').get('VcfMultiClean').get('out_dir'), "{vcf_name}_VcfWgsArrayCommon.{ext}"),vcf_name=out_prefix, ext=["vcf.gz", "vcf.gz.tbi"]),
            #variants qc rules
            expand(os.path.join(BASE_OUT,config.get("rules").get("VariantsMissingRate").get("out_dir"), "{vcf_name}_missing.{ext}"), ext=["lmiss"],vcf_name=out_prefix),
            expand(os.path.join(BASE_OUT,config.get("rules").get("comparePopAF").get("out_dir"), "{vcf_name}_{ref_pop}_af_extrDiff.txt"),vcf_name=out_prefix,ref_pop=ref_pop),
            expand(os.path.join(BASE_OUT,config.get("rules").get("VariantsHetRate").get("out_dir"), "{vcf_name}_ToRemHetRate.txt"),vcf_name=out_prefix),
            #samples qc rules
            expand(os.path.join(BASE_OUT,config.get("rules").get("coverage").get("out_dir"), "{out_name}_dp.{ext}"), ext=["idepth"],out_name=PROJECT_NAME),
            expand(os.path.join(BASE_OUT,config.get("rules").get("singletons").get("out_dir"), "{out_name}_singletons_ALL.singletons"),out_name=PROJECT_NAME),
            expand(os.path.join(BASE_OUT,config.get("rules").get("SampleMissingRate").get("out_dir"), "{out_name}_missing_ALL.imiss"),out_name=PROJECT_NAME),
            expand(os.path.join(BASE_OUT,config.get("rules").get("SampleHetRate").get("out_dir"), "{vcf_name}_hetRate.txt"),vcf_name=PROJECT_NAME),
            #PCA rules
            expand(os.path.join(BASE_OUT,config.get("rules").get("kingPCA").get("out_dir"), "{vcf_name}_cleaned.LD0.3.{ext}"), ext=["vcf.gz","vcf.gz.tbi"],vcf_name=out_prefix),
            expand(os.path.join(BASE_OUT,config.get("rules").get("kingPCA").get("out_dir"), "{out_name}_cleaned.LD0.3.{ext}"), ext=["vcf.gz","bed","bim","fam"],out_name=PROJECT_NAME),
            expand(os.path.join(BASE_OUT,config.get("rules").get("kingPCA").get("out_dir"), "{out_name}_cleaned.LD0.3_kingpca{ext}"),ext=["pc.txt","projpc.txt","proj_Dist.txt","proj_popref.txt"],out_name=PROJECT_NAME),

            # #plots
            expand(os.path.join(BASE_OUT,config.get("rules").get("kingPCA").get("out_dir"), "{out_name}_{ext}"),ext=["pca.pdf","pca_projection_on_1000GP.pdf"],out_name=PROJECT_NAME),
            expand(os.path.join(BASE_OUT,config.get("rules").get("comparePopAF").get("out_dir"), "{vcf_name}_{ref_pop}_af.png"),vcf_name=out_prefix, ref_pop=ref_pop),
            expand(os.path.join(BASE_OUT,config.get("rules").get("comparePopAF").get("out_dir"), "{vcf_name}_{ref_pop}_af_extrDiff.pdf"),vcf_name=out_prefix, ref_pop=ref_pop),
            expand(os.path.join(BASE_OUT,config.get("rules").get("SampleHetRate").get("out_dir"), "{out_name}_hetRate.pdf"),out_name=PROJECT_NAME),
            expand(os.path.join(BASE_OUT,config.get("rules").get("SamplePlots").get("out_dir"), "{out_name}_SingCov.{ext}"),out_name=PROJECT_NAME,ext=["pdf","txt"]),
            expand(os.path.join(BASE_OUT,config.get("rules").get("SamplePlots").get("out_dir"), "{out_name}_hetRateByCov_{group}.pdf"),out_name=PROJECT_NAME,group=['sex','cohort']),
            expand(os.path.join(BASE_OUT,config.get("rules").get("SamplePlots").get("out_dir"), "{out_name}_hetRateByMiss_{group}.pdf"),out_name=PROJECT_NAME,group=['sex','cohort']),
            expand(os.path.join(BASE_OUT,config.get("rules").get("SamplePlots").get("out_dir"), "{out_name}_hetRateBySing_{group}.pdf"),out_name=PROJECT_NAME,group=['sex','cohort']),

            #nrdr rules and plots, to work only on AUTOSOMAL chromosomes
            expand(os.path.join(BASE_OUT, config.get('rules').get('NRD').get('out_dir'), "{vcf_name}_NRDR.txt"), vcf_name=out_prefix_autosomal),
            expand(os.path.join(BASE_OUT, config.get('rules').get('NRD').get('out_dir'), "{vcf_name}_NRDRsites.txt"), vcf_name=out_prefix_autosomal),
            expand(os.path.join(BASE_OUT, config.get('rules').get('NRD').get('out_dir'), "{vcf_name}_NRDRsamples.txt"), vcf_name=out_prefix_autosomal),
            expand(os.path.join(BASE_OUT, config.get('rules').get('NRD').get('out_dir'), "{out_name}_NRDR{subset}.txt"), out_name=PROJECT_NAME, subset=['samples','sites']),
            expand(os.path.join(BASE_OUT,config.get("rules").get("SamplePlots").get("out_dir"), "{out_name}_hetRateByNRD_{group}.pdf"),out_name=PROJECT_NAME,group=['sex','cohort','seq']),

            expand(os.path.join(BASE_OUT,config.get("rules").get("getArrayPopAF").get("out_dir"), "{vcf_name}_ARRAY_af_extrDiff.txt"), vcf_name=out_prefix_autosomal),
            expand(os.path.join(BASE_OUT,config.get("rules").get("getArrayPopAF").get("out_dir"), "{vcf_name}_ARRAY_af.png"), vcf_name=out_prefix_autosomal),
            expand(os.path.join(BASE_OUT,config.get("rules").get("getArrayPopAF").get("out_dir"), "{vcf_name}_ARRAY_af_extrDiff.pdf"), vcf_name=out_prefix_autosomal),
            #stats
            expand(os.path.join(config.get("paths").get("base_out"),config.get("rules").get("stats").get("out_dir"),"{vcf_name}_initial.stats"), vcf_name=out_prefix),
            expand(os.path.join(config.get("paths").get("base_out"),config.get("rules").get("stats").get("out_dir"),"{vcf_name}_HWE95call.stats"), vcf_name=out_prefix),
            #release files to use for filtering
            expand(os.path.join(BASE_OUT,config.get("rules").get("release").get("out_dir"), "{out_name}_cleaned.LD0.3_kingpcaprojpc.txt"),out_name=PROJECT_NAME),
            expand(os.path.join(BASE_OUT,config.get("rules").get("release").get("out_dir"), "{out_name}_hetRate.txt"),out_name=PROJECT_NAME),
            expand(os.path.join(BASE_OUT,config.get("rules").get("release").get("out_dir"), "{out_name}_missing_ALL.imiss"),out_name=PROJECT_NAME),
            expand(os.path.join(BASE_OUT,config.get("rules").get("release").get("out_dir"), "{out_name}_SingCov.txt"),out_name=PROJECT_NAME),
            expand(os.path.join(BASE_OUT,config.get("rules").get("release").get("out_dir"), "{vcf_name}_{ref_pop}_af_extrDiff.txt"), vcf_name=out_prefix,ref_pop=ref_pop),
            expand(os.path.join(BASE_OUT, config.get('rules').get('release').get('out_dir'), "{out_name}_NRDRsamples.txt"),out_name=PROJECT_NAME),
            expand(os.path.join(BASE_OUT, config.get('rules').get('release').get('out_dir'), "{out_name}_NRDRsites.txt"),out_name=PROJECT_NAME),
            expand(os.path.join(BASE_OUT,config.get("rules").get("release").get("out_dir"), "{vcf_name}_ARRAY_af_extrDiff.txt"), vcf_name=out_prefix_autosomal)

else :
    rule all:
        input:
            #define target input for all first steps needed (basically, initial data filtering to perform other comparisons)
            expand(os.path.join(BASE_OUT,config.get("rules").get("cleanMissingHwe").get("out_dir"), "{vcf_name}_HWE95call.{ext}"),vcf_name=out_prefix, ext=["vcf.gz", "vcf.gz.tbi", "removed.sites.vcf.gz"]),
            #variants qc rules
            expand(os.path.join(BASE_OUT,config.get("rules").get("VariantsMissingRate").get("out_dir"), "{vcf_name}_missing.{ext}"), ext=["lmiss"],vcf_name=out_prefix),
            expand(os.path.join(BASE_OUT,config.get("rules").get("comparePopAF").get("out_dir"), "{vcf_name}_{ref_pop}_af_extrDiff.txt"),vcf_name=out_prefix,ref_pop=ref_pop),
            expand(os.path.join(BASE_OUT,config.get("rules").get("VariantsHetRate").get("out_dir"), "{vcf_name}_ToRemHetRate.txt"),vcf_name=out_prefix),
            #samples qc rules
            expand(os.path.join(BASE_OUT,config.get("rules").get("coverage").get("out_dir"), "{out_name}_dp.{ext}"), ext=["idepth"],out_name=PROJECT_NAME),
            expand(os.path.join(BASE_OUT,config.get("rules").get("singletons").get("out_dir"), "{out_name}_singletons_ALL.singletons"),out_name=PROJECT_NAME),
            expand(os.path.join(BASE_OUT,config.get("rules").get("SampleMissingRate").get("out_dir"), "{out_name}_missing_ALL.imiss"),out_name=PROJECT_NAME),
            expand(os.path.join(BASE_OUT,config.get("rules").get("SampleHetRate").get("out_dir"), "{vcf_name}_hetRate.txt"),vcf_name=PROJECT_NAME),
            #PCA rules
            expand(os.path.join(BASE_OUT,config.get("rules").get("kingPCA").get("out_dir"), "{vcf_name}_cleaned.LD0.3.{ext}"), ext=["vcf.gz","vcf.gz.tbi"],vcf_name=out_prefix),
            expand(os.path.join(BASE_OUT,config.get("rules").get("kingPCA").get("out_dir"), "{out_name}_cleaned.LD0.3.{ext}"), ext=["vcf.gz","bed","bim","fam"],out_name=PROJECT_NAME),
            expand(os.path.join(BASE_OUT,config.get("rules").get("kingPCA").get("out_dir"), "{out_name}_cleaned.LD0.3_kingpca{ext}"),ext=["pc.txt","projpc.txt","proj_Dist.txt","proj_popref.txt"],out_name=PROJECT_NAME),
            # #plots
            expand(os.path.join(BASE_OUT,config.get("rules").get("kingPCA").get("out_dir"), "{out_name}_{ext}"),ext=["pca.pdf","pca_projection_on_1000GP.pdf"],out_name=PROJECT_NAME),
            expand(os.path.join(BASE_OUT,config.get("rules").get("comparePopAF").get("out_dir"), "{vcf_name}_{ref_pop}_af.png"),vcf_name=out_prefix, ref_pop=ref_pop),
            expand(os.path.join(BASE_OUT,config.get("rules").get("comparePopAF").get("out_dir"), "{vcf_name}_{ref_pop}_af_extrDiff.pdf"),vcf_name=out_prefix, ref_pop=ref_pop),
            expand(os.path.join(BASE_OUT,config.get("rules").get("SampleHetRate").get("out_dir"), "{out_name}_hetRate.pdf"),out_name=PROJECT_NAME),
            expand(os.path.join(BASE_OUT,config.get("rules").get("SamplePlots").get("out_dir"), "{out_name}_SingCov.{ext}"),out_name=PROJECT_NAME,ext=["pdf","txt"]),
            expand(os.path.join(BASE_OUT,config.get("rules").get("SamplePlots").get("out_dir"), "{out_name}_hetRateByCov_{group}.pdf"),out_name=PROJECT_NAME,group=['sex','cohort']),
            expand(os.path.join(BASE_OUT,config.get("rules").get("SamplePlots").get("out_dir"), "{out_name}_hetRateByMiss_{group}.pdf"),out_name=PROJECT_NAME,group=['sex','cohort']),
            expand(os.path.join(BASE_OUT,config.get("rules").get("SamplePlots").get("out_dir"), "{out_name}_hetRateBySing_{group}.pdf"),out_name=PROJECT_NAME,group=['sex','cohort']),

            #stats
            expand(os.path.join(config.get("paths").get("base_out"),config.get("rules").get("stats").get("out_dir"),"{vcf_name}_initial.stats"), vcf_name=out_prefix),
            expand(os.path.join(config.get("paths").get("base_out"),config.get("rules").get("stats").get("out_dir"),"{vcf_name}_HWE95call.stats"), vcf_name=out_prefix),
            #release files to use for filtering
            expand(os.path.join(BASE_OUT,config.get("rules").get("release").get("out_dir"), "{out_name}_cleaned.LD0.3_kingpcaprojpc.txt"),out_name=PROJECT_NAME),
            expand(os.path.join(BASE_OUT,config.get("rules").get("release").get("out_dir"), "{out_name}_hetRate.txt"),out_name=PROJECT_NAME),
            expand(os.path.join(BASE_OUT,config.get("rules").get("release").get("out_dir"), "{out_name}_missing_ALL.imiss"),out_name=PROJECT_NAME),
            expand(os.path.join(BASE_OUT,config.get("rules").get("release").get("out_dir"), "{out_name}_SingCov.txt"),out_name=PROJECT_NAME),
            expand(os.path.join(BASE_OUT,config.get("rules").get("release").get("out_dir"), "{vcf_name}_{ref_pop}_af_extrDiff.txt"), vcf_name=out_prefix,ref_pop=ref_pop)


#### Modules ####
include:
    include_prefix + "/preproc.smk"
include:
    include_prefix + "/variant_qc.smk"
include:
    include_prefix + "/pca.smk"
include:
    include_prefix + "/vcf_stats.smk"
include:
    include_prefix + "/sample_qc.smk"

if SNP_DATA != "NONE" :
    include:
        include_prefix + "/nrdr.smk"

include:
    include_prefix + "/plots.smk"

include:
    include_prefix + "/release.smk"
