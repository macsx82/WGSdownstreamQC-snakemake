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

# small function to determine if we need to flag a record for removal
# def flag_record_remove(value, up_thr, down_thr):
def flag_record_remove(value, thresholds, comp):
	#comp values can be:
	# lt -> record flagged if the value is lower than the threshold
	# gt -> record flagged if the value is higher than the threshold
	# both -> record flagged if the value in the range defined by the thresholds
	#we know that this will be a list of one or two values
	if len(thresholds) > 1:
		up_thr=thresholds[0]
		down_thr=thresholds[1]
		if (value > up_thr or value < down_thr):
			excl_flag=1
		else:
			excl_flag=0
	else :
		thr=thresholds[0]
		if comp == 'lt':
			if (value < thr):
				excl_flag=1
			else:
				excl_flag=0
		elif comp == 'gt' :
			if (value > thr):
				excl_flag=1
			else:
				excl_flag=0
	return excl_flag

#function to calculate the het rate by site from vcftools --hardy output
def get_hetrate_variants(obs_col):
	#split genotypes
	rr=int(obs_col.split("/")[0])
	ra=int(obs_col.split("/")[1])
	aa=int(obs_col.split("/")[2])
	#count hom genotypes and total genotypes
	sum_hom=rr+aa
	sum_geno=sum_hom+ra
	#get genotype heterozygosity rate
	if sum_geno == 0:
		print(obs_col)
	het_rate=ra/sum_geno
	return het_rate

#function to flag samples with heterozigosity rate higher or lower than 5 SD from the mean
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
	het_up = het_mean + 5 * het_sd
	het_down = het_mean - 5 * het_sd
	#now flag samples for exclusion if their het rate is outside the defined boundaries
	het_df['het_rem']=het_df['het_rate'].apply(lambda x: flag_record_remove(x, [het_up, het_down], 'both'))
	het_df.to_csv(out_file,sep="\t", index=False, header=True, float_format="%.4f")

#function to flag samples with too many singletons based on the singleton distribution
def get_singleton_sample_outliers(singleton_table, out_file):
	#we need to read the het table and calculate:
	het_df = pd.read_table(het_table, sep="\t", header=0)
	# 1) het rate
	het_df['het_rate'] = (het_df['N_SITES'] - het_df['O(HOM)']) / het_df['N_SITES']
	# 2) mean
	het_mean = het_df['het_rate'].mean()
	# 3) sd
	het_sd = het_df['het_rate'].std()
	# 4) upper and lower threshold for samples flag
	het_up = het_mean + 5 * het_sd
	het_down = het_mean - 5 * het_sd
	#now flag samples for exclusion if their het rate is outside the defined boundaries
	het_df['het_rem']=het_df['het_rate'].apply(lambda x: flag_record_remove(x, het_up, het_down))
	het_df.to_csv(out_file,sep="\t", index=False, header=True, float_format="%.4f")


#function to flag variants with heterozigosity rate higher or lower than 5 SD from the mean
def get_het_hwe_variants_outliers(het_table, hwe_thr, out_file):
	#we need to read the het table and calculate:
	# het_table="/large/___SCRATCH___/burlo/cocca/WGS_JOINT_CALL/WGS_QC_pre_release/03.het/ALL_JOINT_744_callset_sorted_variants_het.hwe"
	het_df = pd.read_table(het_table, sep="\t", header=0)
	# 1) Flag sites for removal for HWE pval threshold
	# hwe_thr=0.00000001
	het_df['hwe_rem']=het_df['P_HWE'].apply(lambda x: flag_record_remove(x, [hwe_thr] , 'lt'))
	# 1) het rate
	het_df['hwe_rate']=het_df['OBS(HOM1/HET/HOM2)'].apply(lambda x: get_hetrate_variants(x))
	# 2) mean
	het_mean = het_df['het_rate'].mean()
	# 3) sd
	het_sd = het_df['het_rate'].std()
	# 4) upper and lower threshold for samples flag
	het_up = het_mean + 5 * het_sd
	het_down = het_mean - 5 * het_sd
	#now flag samples for exclusion if their het rate is outside the defined boundaries
	het_df['het_rem']=het_df['het_rate'].apply(lambda x: flag_record_remove(x, het_up, het_down))
	het_df.to_csv(out_file,sep="\t", index=False, header=True, float_format="%.4f")