#function useful for the pipeline
#import libraries
import errno
import gzip 
import io
import multiprocessing
import os
import os.path
import pandas as pd
import pathlib
import psutil
import re
import sys



#define some standard functions to retrieve files more easily
def references_abs_path(ref='references'):
    references = config.get(ref)
    basepath = references['basepath']
    provider = references['provider']
    release = references['release']

    return [os.path.join(basepath, provider, release)]

def resolve_single_filepath(basepath, filename):
    return [os.path.join(basepath, filename)]

def resolve_multi_filepath(basepath, dictionary):
    for k, v in dictionary.items():
        dictionary[k] = os.path.join(basepath, v)
    return dictionary

# define a function to create and open a file in write mode
# got from https://stackoverflow.com/a/30582525
def createAndOpen(filename, mode):
	os.makedirs(os.path.dirname(filename), exist_ok=True)
	return open(filename, mode)

# small function to determine if we need to flag a sample for removal
def het_sample_remove(het_value, het_up_thr, het_down_thr):
	if (het_value > het_up_thr or het_value < het_down_thr):
		excl_flag=1
	else:
		excl_flag=0
	return excl_flag

#function to flag samples with heterozigosity rate higher or lower than 3 SD from the mean
def get_het_sample_outliers(het_table, out_file):
	#we need to read the het table and calculate:
	het_df = pd.read_table(het_table, sep="\t", header=0)
	# 1) het rate
	het_df['het_rate'] = (het_df['N_SITES'] - het_df['O(HOM)']) / het_df['N_SITES']
	# 2) mean
	het_mean = het_df['het_rate'].mean()
	# 3) sd
	het_sd = het_df['het_rate'].std()
	# 4) upper and lower threshold for samples flag
	het_up = het_mean + 4 * het_sd
	het_down = het_mean - 4 * het_sd
	#now flag samples for exclusion if their het rate is outside the defined boundaries
	het_df['het_rem']=het_df['het_rate'].apply(lambda x: het_sample_remove(x, het_up, het_down))
	het_df.to_csv(out_file,sep="\t", index=False, header=True, float_format="%.4f")