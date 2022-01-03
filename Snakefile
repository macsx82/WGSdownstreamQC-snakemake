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
##### functions #####
include:
    "scripts/functions.py"

##### local rules #####
localrules: all

##### Target rules #####
rule all:
	input:
		#define target input
        # expand(os.path.join(BASE_OUT,config.get("rules").get("singletons").get("out_dir"), "{vcf_name}_singletons.{ext}"), ext=["singletons", "log"])
        os.path.join(BASE_OUT,config.get("rules").get("mergeReapplyVQSR").get("out_dir"), PROJECT_NAME + "_MERGED_VQSLODrefilter.vcf.gz")

#### Modules ####
include:
    include_prefix + "/preproc.smk"
# include:
#     include_prefix + "/sample_qc.smk"
# include:
#     include_prefix + "/variant_qc.smk"
# include:
#     include_prefix + "/nrdr.smk"

