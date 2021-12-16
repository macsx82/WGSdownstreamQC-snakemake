#function useful for the pipeline
#import libraries
import pandas as pd
import pathlib
import io
import os
import re


#function to flag samples with heterozigosity rate higher or lower than 3 SD from the mean
def get_het_sample_outliers(wildcards):
	#we need to read the het table and calculate:
	# 1) het rate
	# 2) mean
	# 3) sd
	het_file = pd.read_table(, sep="\t", header=0, dtype='object')	