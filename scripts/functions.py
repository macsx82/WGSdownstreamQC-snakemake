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
import matplotlib
import numpy as np


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
	het_df.to_csv(out_file,sep="\t", index=False, header=True, float_format="%.6f")


#function to flag variants with heterozigosity rate higher or lower than 5 SD from the mean
def get_het_hwe_variants_outliers(het_table, exc_het_thr, out_file):
	#we need to read the het table and calculate:
	# het_table="/large/___SCRATCH___/burlo/cocca/WGS_JOINT_CALL/WGS_QC_pre_release/20220105/02.sites/HetRate/WGS_ITA_PREREL_MERGED_hwe.hwe"
	het_df = pd.read_table(het_table, sep="\t", header=0)
	# 1) From VCFtools we already have info on excess of heterozigosity, in the form of a P_HET_EXCESS p value.
	#    We can use that value and its distribution, to clean the data. Since we need to set a threshold we could:
	#		a) Remove stuff in the first percentile, this way we will clean the data, removing snps putatively associated to excess of heterozygosity
	#		b) Remove all sites with a pvalue < 1e-08, using the classic genome wide significant threshold, even if it doesn't have really sense, here.
	# hwe_thr=0.00000001
	# exc_het_thr=het_df['P_HET_EXCESS'].quantile(0.001)
	# exc_het_thr=1e-8
	to_rem_df=het_df[het_df.P_HET_EXCESS < exc_het_thr]
	to_rem_df.to_csv(out_file,sep="\t", index=False, header=True)


#function to generate a merged file with AF to plot
def af_diff(wgs_table, ext_table, outfile, outplot,outplot_extreme):
	#try to fix X11 error
	matplotlib.use('Agg')
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
	#plot only data from the extreme diff dataset
	extreme_diff_variants_df.plot.scatter(x='AF_wgs',y='AF_ext',s=extreme_diff_variants_df['af_diff'] * 100)
	plt.xlabel("WGS AF")
	plt.ylabel("EXT dataset AF")
	plt.savefig(outplot_extreme)
	#plot also the data, defining the point size based on the diff value
	merged_df.plot.scatter(x='AF_x',y='AF_y',s=merged_df['af_diff'] * 100)
	plt.xlabel("WGS AF")
	plt.ylabel("EXT dataset AF")
	# plt.savefig('test.pdf')
	plt.savefig(outplot)


#function to generate a plot of N singletons vs coverage
def plot_sing_vs_cov(sing_table, cov_table, outplot,out_table):
	#try to fix X11 error
	matplotlib.use('Agg')
	# sing_table="/large/___SCRATCH___/burlo/cocca/WGS_JOINT_CALL/WGS_QC_pre_release/20220110/03.samples/singletons/WGS_ITA_PREREL_MERGED_singletons.singletons"
	# cov_table="/large/___SCRATCH___/burlo/cocca/WGS_JOINT_CALL/WGS_QC_pre_release/20220110/03.samples/coverage/WGS_ITA_PREREL_MERGED_dp.idepth"
	sing_df = pd.read_table(sing_table,sep="\t", header=0)
	cov_df = pd.read_table(cov_table,sep="\t", header=0)

	#calculate singleton number by sample. We count as singletons also doubleton sites
	sing_number=sing_df.groupby(sing_df['INDV'], as_index=False).size()
	# merge dataframes using the provided key
	merged_df = sing_number.merge(cov_df, how='inner',on='INDV')
	#rename columns
	merged_df.columns = ['INDV','SINGLETONS','N_SITES','MEAN_DEPTH']
	# tag also samples with putative singletons excess based on singletons distributions
	#mean
	sing_mean = merged_df['SINGLETONS'].mean()
	#sd
	sing_sd = merged_df['SINGLETONS'].std()
	# 4) upper and lower threshold for samples flag
	sing_up = sing_mean + 3 * sing_sd
	sing_down = sing_mean - 3 * sing_sd
	#now flag samples for exclusion if their het rate is outside the defined boundaries
	merged_df['sing_rem']=merged_df['SINGLETONS'].apply(lambda x: flag_record_remove(x, [sing_up, sing_down], 'both'))
	#add singleton_rate column (over total number sites per sample)
	merged_df['sing_rate']= merged_df.SINGLETONS/merged_df.N_SITES
	#save also the merged table
	merged_df.to_csv(out_table,sep="\t", index=False, header=True, float_format="%.8f")
	#plot the data, defining the point size based on the diff value
	merged_df.plot.scatter(x='SINGLETONS',y='MEAN_DEPTH',s=merged_df['sing_rate'] * 20000)
	plt.axvline(x=sing_up, color='r', label=' Singletons upper threshold (3SD)')
	plt.axvline(x=sing_down, color='b', label=' Singletons lower threshold (3SD)')
	plt.xlabel("SINGLETONS")
	plt.ylabel("MEAN DEPTH")
	plt.title("Singletons vs Mean Depth")
	# plt.savefig('test.pdf')
	plt.savefig(outplot)

#function to generate a plot of N singletons vs coverage
def plot_het_rate_sample(het_rate_table, outplot):
	#try to fix X11 error
	matplotlib.use('Agg')
	# het_rate_table="/large/___SCRATCH___/burlo/cocca/WGS_JOINT_CALL/WGS_QC_pre_release/20220105/03.samples/HetRate/WGS_ITA_PREREL_MERGED_hetRate.txt"
	het_rate_df = pd.read_table(het_rate_table,sep="\t", header=0)
	#get some values to plot
	het_rate_mean = het_rate_df['het_rate'].mean()
	#sd
	het_rate_sd = het_rate_df['het_rate'].std()
	# 4) upper and lower threshold for samples flag
	thr_up = het_rate_mean + 3 * het_rate_sd
	thr_down = het_rate_mean - 3 * het_rate_sd
	thr_up5 = het_rate_mean + 5 * het_rate_sd
	thr_down5 = het_rate_mean - 5 * het_rate_sd
	#get all values tagged for removal and add labels to the points
	het_rate_rem=list(het_rate_df[het_rate_df['het_rem']==1]['INDV'])
	#plot the data, defining the point size based on the diff value
	het_rate_df.plot.scatter(x='INDV',y='het_rate',s=het_rate_df['het_rate'] * 200)
	plt.axhline(y=thr_up, color='red', label='Het Rate upper threshold (3SD)')
	plt.axhline(y=thr_up5, color='green', label='Het Rate upper threshold (5SD)')
	plt.axhline(y=thr_down, color='blue', label='Het Rate lower threshold (3SD)')
	plt.axhline(y=thr_down, color='black', label='Het Rate lower threshold (5SD)')
	frame1=plt.gca()
	frame1.axes.get_xaxis().set_ticks([])
	plt.xlabel("Samples")
	plt.ylabel("Het rate")
	plt.title("Het Rate per sample")
	for s_label in het_rate_rem:
		s_name=het_rate_df[het_rate_df['INDV']==s_label]['INDV']
		s_rate=het_rate_df[het_rate_df['INDV']==s_label]['het_rate']
		plt.annotate(s_label,(s_name, s_rate))
	# plt.savefig('testHETRATE.pdf')
	plt.savefig(outplot)


def plot_het_rate_vs_coverage(het_rate_table,cov_table,manifest_table,outplot_prefix):
	#try to fix X11 error
	matplotlib.use('Agg')
	# het_rate_table="/large/___SCRATCH___/burlo/cocca/WGS_JOINT_CALL/WGS_QC_pre_release/20220105/03.samples/HetRate/WGS_ITA_PREREL_MERGED_hetRate.txt"
	het_rate_df = pd.read_table(het_rate_table,sep="\t", header=0)
	# cov_table="/large/___SCRATCH___/burlo/cocca/WGS_JOINT_CALL/WGS_QC_pre_release/20220110/03.samples/coverage/WGS_ITA_PREREL_MERGED_dp.idepth"
	cov_df = pd.read_table(cov_table,sep="\t", header=0)
	# manifest_table="/large/___HOME___/burlo/cocca/analyses/WGS_QC_pre_release/WGS_817_samples_manifest.txt"
	sex_df=pd.read_table(manifest_table,sep=" ", header=0)
	# merge dataframes using the provided key
	merged_df = het_rate_df.merge(cov_df, how='inner',on='INDV')
	#get all values tagged for removal and add labels to the points
	het_rate_rem=list(het_rate_df[het_rate_df['het_rem']==1]['INDV'])
	#define color map for sexes
	colors_map ={'Female': 'green' , 'Male' : 'orange'}
	#add sex to the dataframe and convert to string value
	merged_sex_df = merged_df.merge(sex_df, how='inner', left_on='INDV', right_on='SAMPLE_ID')
	merged_sex_df['sex'].replace({1:'Male',2:'Female'}, inplace=True)
	#group by population on the merged dataset
	merged_df_grouped= merged_sex_df.groupby('sex')
	#initialize the plot
	fig, ax = plt.subplots()
	#plot the data, defining the point size based on the het value
	for key, group in merged_df_grouped:
		group.plot(ax=ax,kind='scatter',x='het_rate',y='MEAN_DEPTH', color=colors_map[str(key)], label=key)
	# merged_df.plot.scatter(y='MEAN_DEPTH',x='het_rate',s=het_rate_df['het_rate'] * 200)
	plt.xlabel("Het rate")
	plt.ylabel("Mean depth")
	plt.title("Het Rate by mean depth per sample")
	for s_label in het_rate_rem:
		s_rate=merged_sex_df[merged_sex_df['INDV']==s_label]['het_rate']
		s_depth=merged_sex_df[merged_sex_df['INDV']==s_label]['MEAN_DEPTH']
		plt.annotate(s_label,(s_rate,s_depth))
	# plt.savefig('testHETbyDP.pdf')
	plt.savefig(outplot_prefix+"_sex.pdf")
	#need to plot also by cohort
	merged_df_grouped_cohort= merged_sex_df.groupby('COHORT')
	# get cohorts
	cohorts_list=list(merged_df_grouped_cohort.groups.keys())
	#Map cohorts to colors
	vals = np.linspace(0,1,len(cohorts_list))
	np.random.shuffle(vals)
	#use the tab20 colormap
	colors_cohort_map= dict(zip(cohorts_list,plt.cm.tab20(vals)))
	# cmap = plt.cm.colors.ListedColormap(plt.cm.tab20(vals))
	# colors_map = dict(zip(pops,colors))
	#initialize the plot
	fig, ax = plt.subplots()
	#plot the data, defining the point size based on the het value
	for key, group in merged_df_grouped_cohort:
		group.plot(ax=ax,kind='scatter',x='het_rate',y='MEAN_DEPTH', color=colors_cohort_map[str(key)], label=key)
	# merged_df.plot.scatter(y='MEAN_DEPTH',x='het_rate',s=het_rate_df['het_rate'] * 200)
	plt.xlabel("Het rate")
	plt.ylabel("Mean depth")
	plt.title("Het Rate by mean depth per sample")
	for s_label in het_rate_rem:
		s_rate=merged_sex_df[merged_sex_df['INDV']==s_label]['het_rate']
		s_depth=merged_sex_df[merged_sex_df['INDV']==s_label]['MEAN_DEPTH']
		plt.annotate(s_label,(s_rate,s_depth))
	ax.legend(loc='upper right', ncol=3, fontsize='small')
	# plt.savefig('testHETbyDPcohort.pdf')
	plt.savefig(outplot_prefix+"_cohort.pdf")