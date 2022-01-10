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
import matplotlib.pyplot as plt


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
	# if sum_geno == 0:
	# 	print(obs_col)
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


#function to generate a merged file with AF to plot
def af_diff(wgs_table, ext_table, outfile, outplot):
	# define the column header
	col_head=['var_key','CHROM','POS','ID','REF','ALT','AC','AN','AF','MAF']
	# wgs_table="/large/___SCRATCH___/burlo/cocca/WGS_JOINT_CALL/WGS_QC_pre_release/20220105/01.VQSR_reapply/test_snps.tab"
	# ext_table="/storage/burlo/cocca/resources/1000GP_phase3_3202/vcf/EUR/EUR_normIndel_multiSplitted.vcf.gz.tab"
	wgs_df = pd.read_table(wgs_table,sep="\t", header=None,names=col_head)
	ext_df = pd.read_table(ext_table,sep="\t", header=None,names=col_head)
	# merge dataframes using the provided key
	merged_df = wgs_df.merge(ext_df, how='inner',on='var_key')
	#calculate the af difference
	merged_df['af_diff'] = abs(merged_df['AF_x'] - merged_df['AF_y'])
	#get the threshold value for the 99 percentile, which should contain the most differentiated snps in terms of AF
	extreme_af_diff_thr=merged_df['af_diff'].quantile(0.99)
	#now find all variants with the highest differences
	extreme_diff_variants_df=merged_df[merged_df['af_diff'] >= extreme_af_diff_thr]
	#now write the resulting table with the highest diff AF
	extreme_diff_variants_df=extreme_diff_variants_df[['var_key', 'CHROM_x', 'POS_x', 'ID_x', 'REF_x', 'ALT_x', 'AC_x', 'AN_x', 'AF_x', 'MAF_x', 'AC_y', 'AN_y', 'AF_y', 'MAF_y', 'af_diff']]
	#rename columns
	extreme_diff_variants_df.columns = ['var_key', 'CHROM', 'POS', 'ID', 'REF', 'ALT', 'AC_wgs', 'AN_wgs', 'AF_wgs', 'MAF_wgs', 'AC_ext', 'AN_ext', 'AF_ext', 'MAF_ext', 'af_diff']
	extreme_diff_variants_df.to_csv(outfile,sep="\t", index=False, header=True, float_format="%.4f")
	#plot also the data, defining the point size based on the diff value
	merged_df.plot.scatter(x='AF_x',y='AF_y',s=merged_df['af_diff'] * 100)
	plt.xlabel("WGS AF")
	plt.ylabel("EXT dataset AF")
	# plt.savefig('test.pdf')
	plt.savefig(outplot)
