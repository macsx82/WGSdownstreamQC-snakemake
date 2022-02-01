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
import collections


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


#function to collect all NRD information to calculate the genome wide NRD by sample
def get_NRD_by_sample(all_NRD, outfile):
	# we have all files splitted by chr, we need to generate a dictionary for each sample and load the data we need
	samples_nrd_dict=collections.defaultdict(lambda: collections.defaultdict(list))
	for sample_nrd in all_NRD:
		# sample_nrd="/large/___SCRATCH___/burlo/cocca/WGS_JOINT_CALL/WGS_QC_pre_release/20220105/04.nrdr/tables/WGS_ITA_PREREL_MERGED_chr1_NRDRsamples.txt"
		current_nrd_sample=open(sample_nrd,'r')
		for sample_line in current_nrd_sample:
			#read all values by chr for the current sample
			c_sample_line=sample_line.strip().split("\t")
			#get the sample
			sample=c_sample_line[2]
			#xRR+xRA+xAA
			mismatches=sum([int(c_sample_line[7]),int(c_sample_line[8]),int(c_sample_line[9])])
			#xRR+xRA+xAA+mRA+mAA
			all_non_ref=sum([mismatches,int(c_sample_line[5]),int(c_sample_line[6])])
			samples_nrd_dict[sample]["mismatches"].append(mismatches)
			samples_nrd_dict[sample]["all_non_ref"].append(all_non_ref)
	
	nrd_dict={}
	#now calculate the wg nrd for all samples and print all in the output file
	for sample in samples_nrd_dict.keys():
		nrd_dict[sample]=(sum(samples_nrd_dict[sample]['mismatches'])/sum(samples_nrd_dict[sample]['all_non_ref']))*100
	#we want to provide a flag for removal for the highest NRD samples (+5sd), so we need a pandas dataframe!!
	nrd_df=pd.DataFrame(list(nrd_dict.items()))
	nrd_df.columns = ['SAMPLE_ID','NRD']
	#get mean and sd
	nrd_mean=nrd_df.NRD.mean()
	nrd_sd=nrd_df.NRD.std()
	# 4) upper threshold for samples flag
	nrd_up = nrd_mean + 3 * nrd_sd
	#now flag samples for exclusion if their het rate is outside the defined boundaries
	nrd_df['nrd_rem']=nrd_df['NRD'].apply(lambda x: flag_record_remove(x, [nrd_up], 'gt'))
	nrd_df.to_csv(outfile,sep="\t", index=False, header=True, float_format="%.6f")


#function to collect all NRD information by site and flag for removal sites with NRD > 3 sd
def get_NRD_by_site(all_NRD, outfile):
	# we have all files splitted by chr, we need to generate a dictionary for each sample and load the data we need
	# all_NRD = ["/large/___SCRATCH___/burlo/cocca/WGS_JOINT_CALL/WGS_QC_pre_release/20220105/04.nrdr/tables/WGS_ITA_PREREL_MERGED_chr1_NRDRsites.txt","/large/___SCRATCH___/burlo/cocca/WGS_JOINT_CALL/WGS_QC_pre_release/20220105/04.nrdr/tables/WGS_ITA_PREREL_MERGED_chr2_NRDRsites.txt","/large/___SCRATCH___/burlo/cocca/WGS_JOINT_CALL/WGS_QC_pre_release/20220105/04.nrdr/tables/WGS_ITA_PREREL_MERGED_chr3_NRDRsites.txt","/large/___SCRATCH___/burlo/cocca/WGS_JOINT_CALL/WGS_QC_pre_release/20220105/04.nrdr/tables/WGS_ITA_PREREL_MERGED_chr4_NRDRsites.txt","/large/___SCRATCH___/burlo/cocca/WGS_JOINT_CALL/WGS_QC_pre_release/20220105/04.nrdr/tables/WGS_ITA_PREREL_MERGED_chr5_NRDRsites.txt","/large/___SCRATCH___/burlo/cocca/WGS_JOINT_CALL/WGS_QC_pre_release/20220105/04.nrdr/tables/WGS_ITA_PREREL_MERGED_chr6_NRDRsites.txt","/large/___SCRATCH___/burlo/cocca/WGS_JOINT_CALL/WGS_QC_pre_release/20220105/04.nrdr/tables/WGS_ITA_PREREL_MERGED_chr7_NRDRsites.txt","/large/___SCRATCH___/burlo/cocca/WGS_JOINT_CALL/WGS_QC_pre_release/20220105/04.nrdr/tables/WGS_ITA_PREREL_MERGED_chr8_NRDRsites.txt","/large/___SCRATCH___/burlo/cocca/WGS_JOINT_CALL/WGS_QC_pre_release/20220105/04.nrdr/tables/WGS_ITA_PREREL_MERGED_chr9_NRDRsites.txt","/large/___SCRATCH___/burlo/cocca/WGS_JOINT_CALL/WGS_QC_pre_release/20220105/04.nrdr/tables/WGS_ITA_PREREL_MERGED_chr10_NRDRsites.txt","/large/___SCRATCH___/burlo/cocca/WGS_JOINT_CALL/WGS_QC_pre_release/20220105/04.nrdr/tables/WGS_ITA_PREREL_MERGED_chr11_NRDRsites.txt","/large/___SCRATCH___/burlo/cocca/WGS_JOINT_CALL/WGS_QC_pre_release/20220105/04.nrdr/tables/WGS_ITA_PREREL_MERGED_chr12_NRDRsites.txt","/large/___SCRATCH___/burlo/cocca/WGS_JOINT_CALL/WGS_QC_pre_release/20220105/04.nrdr/tables/WGS_ITA_PREREL_MERGED_chr13_NRDRsites.txt","/large/___SCRATCH___/burlo/cocca/WGS_JOINT_CALL/WGS_QC_pre_release/20220105/04.nrdr/tables/WGS_ITA_PREREL_MERGED_chr14_NRDRsites.txt","/large/___SCRATCH___/burlo/cocca/WGS_JOINT_CALL/WGS_QC_pre_release/20220105/04.nrdr/tables/WGS_ITA_PREREL_MERGED_chr15_NRDRsites.txt","/large/___SCRATCH___/burlo/cocca/WGS_JOINT_CALL/WGS_QC_pre_release/20220105/04.nrdr/tables/WGS_ITA_PREREL_MERGED_chr16_NRDRsites.txt","/large/___SCRATCH___/burlo/cocca/WGS_JOINT_CALL/WGS_QC_pre_release/20220105/04.nrdr/tables/WGS_ITA_PREREL_MERGED_chr17_NRDRsites.txt","/large/___SCRATCH___/burlo/cocca/WGS_JOINT_CALL/WGS_QC_pre_release/20220105/04.nrdr/tables/WGS_ITA_PREREL_MERGED_chr18_NRDRsites.txt","/large/___SCRATCH___/burlo/cocca/WGS_JOINT_CALL/WGS_QC_pre_release/20220105/04.nrdr/tables/WGS_ITA_PREREL_MERGED_chr19_NRDRsites.txt","/large/___SCRATCH___/burlo/cocca/WGS_JOINT_CALL/WGS_QC_pre_release/20220105/04.nrdr/tables/WGS_ITA_PREREL_MERGED_chr20_NRDRsites.txt","/large/___SCRATCH___/burlo/cocca/WGS_JOINT_CALL/WGS_QC_pre_release/20220105/04.nrdr/tables/WGS_ITA_PREREL_MERGED_chr21_NRDRsites.txt","/large/___SCRATCH___/burlo/cocca/WGS_JOINT_CALL/WGS_QC_pre_release/20220105/04.nrdr/tables/WGS_ITA_PREREL_MERGED_chr22_NRDRsites.txt"]
	# samples_nrd_dict=collections.defaultdict(lambda: collections.defaultdict(list))
	sites_nrd_df= pd.DataFrame()
	for chr_nrd in all_NRD:
		nrd_df = pd.read_table(chr_nrd,sep="\t", header=None)
		sites_nrd_df=sites_nrd_df.append(nrd_df)
	sites_nrd_df=sites_nrd_df.drop(axis='columns',labels=[0])
	sites_nrd_df.columns=['CHROM','POS','MATCHES','MISMATCHES','NRD']
	#get mean and sd
	nrd_mean=sites_nrd_df.NRD.mean()
	nrd_sd=sites_nrd_df.NRD.std()
	#upper threshold for samples flag
	nrd_up = nrd_mean + 3 * nrd_sd
	#now flag samples for exclusion if their het rate is outside the defined boundaries
	sites_nrd_df['nrd_rem']=sites_nrd_df['NRD'].apply(lambda x: flag_record_remove(x, [nrd_up], 'gt'))
	sites_nrd_df.to_csv(outfile,sep="\t", index=False, header=True, float_format="%.6f")

#function to collect and merge het rate data per sample splitted by chromosome
def collectSampleHetRate(all_het, outfile):
	# we have all files splitted by chr, we need to generate a dictionary for each sample and load the data we need
	samples_het_dict=collections.defaultdict(lambda: collections.defaultdict(list))
	for sample_het in all_het:
		# sample_nrd="/large/___SCRATCH___/burlo/cocca/WGS_JOINT_CALL/WGS_QC_pre_release/20220105/04.nrdr/tables/WGS_ITA_PREREL_MERGED_chr1_NRDRsamples.txt"
		current_het_sample=open(sample_het,'r')
		for sample_line in current_het_sample:
			if not(re.match("INDV",sample_line.strip())): 
				#read all values by chr for the current sample
				c_sample_line=sample_line.strip().split("\t")
				#get the sample
				sample=c_sample_line[0]
				#xRR+xRA+xAA
				observed_h=int(c_sample_line[1])
				expected_h=float(c_sample_line[2])
				n_sites=int(c_sample_line[3])
				#xRR+xRA+xAA+mRA+mAA
				samples_het_dict[sample]["observed_h"].append(observed_h)
				samples_het_dict[sample]["expected_h"].append(expected_h)
				samples_het_dict[sample]["n_sites"].append(n_sites)
	
	her_dict={}
	#now calculate the wg nrd for all samples and print all in the output file
	for sample in samples_het_dict.keys():
		her_dict[sample]=[sum(samples_het_dict[sample]["observed_h"]),sum(samples_het_dict[sample]["expected_h"]),sum(samples_het_dict[sample]["n_sites"])]
	#we want to provide a flag for removal for the highest NRD samples (+5sd), so we need a pandas dataframe!!
	het_df=pd.DataFrame(list(her_dict.items()))
	het_df.columns = ['INDV','O(HOM)','E(HOM)','N_SITES']
	# #get mean and sd
	# nrd_mean=nrd_df.NRD.mean()
	# nrd_sd=nrd_df.NRD.std()
	# # 4) upper threshold for samples flag
	# nrd_up = nrd_mean + 3 * nrd_sd
	# #now flag samples for exclusion if their het rate is outside the defined boundaries
	# nrd_df['nrd_rem']=nrd_df['NRD'].apply(lambda x: flag_record_remove(x, [nrd_up], 'gt'))
	het_df.to_csv(outfile,sep="\t", index=False, header=True, float_format="%.6f")