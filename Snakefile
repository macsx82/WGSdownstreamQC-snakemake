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
from snakemake.exceptions import print_exception, WorkflowError
from snakemake.utils import validate, min_version

##### set minimum snakemake version #####
min_version("5.32.0")

include_prefix="rules"

#define some global variables from the config file
BASE_OUT=config["paths"]["base_out"]
MAIN_VCF_INPUT=config["paths"]["input_vcf"]
chroms=config["chrs"]
PROJECT_NAME=config["project_name"]
out_prefix=PROJECT_NAME + "_MERGED"
### path to resources needed for plots
tgRefBed = config['paths']['1000G_ref_for_king']

##### functions #####
include:
    "scripts/functions.py"

##### local rules #####
localrules: all

##### Target rules #####
rule all:
    input:
        #define target input
        expand(os.path.join(BASE_OUT,config.get("rules").get("mergeReapplyVQSR").get("out_dir"), "{vcf_name}_VQSLODrefilter.vcf.gz"), vcf_name=out_prefix),
        expand(os.path.join(config.get("paths").get("base_out"),config.get("rules").get("stats").get("out_dir"),"{vcf_name}_initial.stats"), vcf_name=out_prefix),
        #Variants qc rules
        expand(os.path.join(BASE_OUT,config.get("rules").get("cleanMissingHwe").get("out_dir"), "{vcf_name}_HWE95call.{ext}"),vcf_name=out_prefix, ext=["vcf.gz", "vcf.gz.tbi", "removed.sites","log"]),
        expand(os.path.join(BASE_OUT,config.get("rules").get("VariantsMissingRate").get("out_dir"), "{vcf_name}_missing.{ext}"), ext=["lmiss", "log"],vcf_name=out_prefix),
        expand(os.path.join(BASE_OUT,config.get("rules").get("comparePopAF").get("out_dir"), "{vcf_name}_{ext_ref}_{suffix}"),vcf_name=out_prefix, suffix=["af_extrDiff.txt","af.pdf"],ext_ref=list(config.get("rules").get("comparePopAF").get("ref_pops").keys())),

        #samples qc rules
        expand(os.path.join(BASE_OUT,config.get("rules").get("singletons").get("out_dir"), "{vcf_name}_singletons.{ext}"), ext=["singletons", "log"],vcf_name=out_prefix),
        expand(os.path.join(BASE_OUT,config.get("rules").get("coverage").get("out_dir"), "{vcf_name}_dp.{ext}"), ext=["idepth", "log"],vcf_name=out_prefix),
        expand(os.path.join(BASE_OUT,config.get("rules").get("SampleMissingRate").get("out_dir"), "{vcf_name}_missing.{ext}"), ext=["imiss", "log"],vcf_name=out_prefix),
        # expand(os.path.join(BASE_OUT,config.get("rules").get("SampleHetRate").get("out_dir"), "{vcf_name}_het.{ext}"), ext=["het", "log"],vcf_name=out_prefix),
        expand(os.path.join(BASE_OUT,config.get("rules").get("SampleHetRate").get("out_dir"), "{vcf_name}_hetRate.txt"),vcf_name=out_prefix)
#### Modules ####
include:
    include_prefix + "/preproc.smk"
include:
    include_prefix + "/variant_qc.smk"
include:
    include_prefix + "/vcf_stats.smk"
include:
    include_prefix + "/sample_qc.smk"
# include:
#     include_prefix + "/nrdr.smk"

