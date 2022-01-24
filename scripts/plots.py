#plot function useful for the pipeline
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
	#add labels for extreme values
	out_sing=list(merged_df[merged_df['SINGLETONS'] >= sing_up].INDV)
	for s_label in out_sing:
		s_sing=merged_df[merged_df['INDV']==s_label]['SINGLETONS']
		s_dp=merged_df[merged_df['INDV']==s_label]['MEAN_DEPTH']
		plt.annotate(s_label,(s_sing, s_dp))
	plt.axvline(x=sing_up, color='r', label=' Singletons upper threshold (3SD)')
	plt.axvline(x=sing_down, color='b', label=' Singletons lower threshold (3SD)')
	plt.xlabel("SINGLETONS")
	plt.ylabel("MEAN DEPTH")
	plt.title("Singletons vs Mean Depth")
	# plt.savefig('test.pdf')
	plt.savefig(outplot)

#function to generate a plot of N singletons vs coverage
def plot_het_rate_sample(het_rate_table, outplot_prefix):
	#try to fix X11 error
	matplotlib.use('Agg')
	# het_rate_table="/large/___SCRATCH___/burlo/cocca/WGS_JOINT_CALL/WGS_QC_pre_release/20220105/03.samples/HetRate/WGS_ITA_PREREL_MERGED_hetRate.txt"
	het_rate_df = pd.read_table(het_rate_table,sep="\t", header=0)
	# manifest_table="/large/___HOME___/burlo/cocca/analyses/WGS_QC_pre_release/WGS_817_samples_SEQ_CENTRE_manifest.txt"
	manifest_df=pd.read_table(manifest_table,sep=" ", header=0)
	#get some values to plot
	het_rate_mean = het_rate_df['het_rate'].mean()
	#sd
	het_rate_sd = het_rate_df['het_rate'].std()
	# 4) upper and lower threshold for samples flag
	thr_up1 = het_rate_mean + het_rate_sd
	thr_up = het_rate_mean + 3 * het_rate_sd
	thr_down = het_rate_mean - 3 * het_rate_sd
	thr_up5 = het_rate_mean + 5 * het_rate_sd
	thr_down5 = het_rate_mean - 5 * het_rate_sd
	#get all values tagged for removal and add labels to the points
	het_rate_rem=list(het_rate_df[het_rate_df['het_rem']==1]['INDV'])
	#merge the manifest so we can plot using also the seq centre
	merged_df=het_rate_df.merge(manifest_df,how='inner', left_on='INDV', right_on='SAMPLE_ID')
	#need to plot also by cohort
	merged_df_grouped_seq= merged_df.groupby('SEQ')
	# get seq centre
	seq_list=list(merged_df_grouped_seq.groups.keys())
	#Map cohorts to colors
	vals = np.linspace(0,1,len(seq_list))
	np.random.shuffle(vals)
	#use the tab20 colormap
	colors_seq_map= dict(zip(seq_list,plt.cm.tab20(vals)))
	#initialize the plot
	fig, ax = plt.subplots()
	#plot the data, defining the point size based on the het value
	for key, group in merged_df_grouped_seq:
		group.plot(ax=ax,kind='scatter',x='INDV',y='het_rate', color=colors_seq_map[str(key)], label=key)
	# merged_df.plot.scatter(y='MEAN_DEPTH',x='het_rate',s=het_rate_df['het_rate'] * 200)
	plt.xlabel("Samples")
	plt.ylabel("Het rate")
	plt.title("Het Rate  per sample by seq centre")
	for s_label in het_rate_rem:
		s_rate=merged_df[merged_df['INDV']==s_label]['INDV']
		s_depth=merged_df[merged_df['INDV']==s_label]['het_rate']
		plt.annotate(s_label,(s_rate,s_depth))
	ax.legend(loc='upper right', ncol=3, fontsize='small')
	# plt.savefig('testHETbySEQ.pdf')
	plt.savefig(outplot_prefix+"_SEQ.pdf")
	#plot the data, defining the point size based on the diff value
	het_rate_df.plot.scatter(x='INDV',y='het_rate',s=het_rate_df['het_rate'] * 200)
	plt.axhline(y=thr_up1, color='yellow', label='Het Rate upper threshold (1SD)')
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
	plt.savefig(outplot_prefix+".pdf")
	#we could also add the het rate density distribution



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


def plot_het_rate_vs_missing(het_rate_table,missing_table,manifest_table,outplot_prefix):
	#try to fix X11 error
	matplotlib.use('Agg')
	# het_rate_table="/large/___SCRATCH___/burlo/cocca/WGS_JOINT_CALL/WGS_QC_pre_release/20220105/03.samples/HetRate/WGS_ITA_PREREL_MERGED_hetRate.txt"
	het_rate_df = pd.read_table(het_rate_table,sep="\t", header=0)
	# missing_table="/large/___SCRATCH___/burlo/cocca/WGS_JOINT_CALL/WGS_QC_pre_release/20220105/03.samples/MissingRate/WGS_ITA_PREREL_MERGED_missing.imiss"
	missing_df = pd.read_table(missing_table,sep="\t", header=0)
	# manifest_table="/large/___HOME___/burlo/cocca/analyses/WGS_QC_pre_release/WGS_817_samples_manifest.txt"
	sex_df=pd.read_table(manifest_table,sep=" ", header=0)
	# merge dataframes using the provided key
	merged_df = het_rate_df.merge(missing_df, how='inner',on='INDV')
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
		group.plot(ax=ax,kind='scatter',x='het_rate',y='F_MISS', color=colors_map[str(key)], label=key)
	# merged_df.plot.scatter(y='MEAN_DEPTH',x='het_rate',s=het_rate_df['het_rate'] * 200)
	plt.xlabel("Het rate")
	plt.ylabel("Missing sites fraction")
	plt.title("Het Rate by Missing fraction per sample")
	for s_label in het_rate_rem:
		s_rate=merged_sex_df[merged_sex_df['INDV']==s_label]['het_rate']
		s_miss=merged_sex_df[merged_sex_df['INDV']==s_label]['F_MISS']
		plt.annotate(s_label,(s_rate,s_miss))
	# plt.savefig('testHETbyMISS.pdf')
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
		group.plot(ax=ax,kind='scatter',x='het_rate',y='F_MISS', color=colors_cohort_map[str(key)], label=key)
	# merged_df.plot.scatter(y='MEAN_DEPTH',x='het_rate',s=het_rate_df['het_rate'] * 200)
	plt.xlabel("Het rate")
	plt.ylabel("Missing sites fraction")
	plt.title("Het Rate by Missing fraction per sample")
	for s_label in het_rate_rem:
		s_rate=merged_sex_df[merged_sex_df['INDV']==s_label]['het_rate']
		s_depth=merged_sex_df[merged_sex_df['INDV']==s_label]['F_MISS']
		plt.annotate(s_label,(s_rate,s_depth))
	ax.legend(loc='upper right', ncol=3, fontsize='small')
	# plt.savefig('testHETbyMISScohort.pdf')
	plt.savefig(outplot_prefix+"_cohort.pdf")

# plot het vs singletons
def plot_het_rate_vs_singletons(het_rate_table,sing_table,manifest_table,outplot_prefix):
	#try to fix X11 error
	matplotlib.use('Agg')
	# het_rate_table="/large/___SCRATCH___/burlo/cocca/WGS_JOINT_CALL/WGS_QC_pre_release/20220105/03.samples/HetRate/WGS_ITA_PREREL_MERGED_hetRate.txt"
	het_rate_df = pd.read_table(het_rate_table,sep="\t", header=0)
	# sing_table="/large/___SCRATCH___/burlo/cocca/WGS_JOINT_CALL/WGS_QC_pre_release/20220105/03.samples/singletons/WGS_ITA_PREREL_MERGED_singletons.singletons"
	sing_df = pd.read_table(sing_table,sep="\t", header=0)
	#calculate singleton number by sample. We count as singletons also doubleton sites
	sing_number=sing_df.groupby(sing_df['INDV'], as_index=False).size()
	#rename columns
	sing_number.columns = ['INDV','SINGLETONS']
	# merge dataframes using the provided key
	merged_df = sing_number.merge(het_rate_df, how='inner',on='INDV')
	# manifest_table="/large/___HOME___/burlo/cocca/analyses/WGS_QC_pre_release/WGS_817_samples_manifest.txt"
	sex_df=pd.read_table(manifest_table,sep=" ", header=0)
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
		group.plot(ax=ax,kind='scatter',x='het_rate',y='SINGLETONS', color=colors_map[str(key)], label=key)
	plt.xlabel("Het rate")
	plt.ylabel("Singleton count")
	plt.title("Het Rate by Singleton count per sample")
	# plt.savefig('testHETbySING.pdf')
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
		group.plot(ax=ax,kind='scatter',x='het_rate',y='SINGLETONS', color=colors_cohort_map[str(key)], label=key)
	plt.xlabel("Het rate")
	plt.ylabel("Singleton count")
	plt.title("Het Rate by Singleton count per sample")
	ax.legend(loc='upper right', ncol=3, fontsize='small')
	# plt.savefig('testHETbySINGcohort.pdf')
	plt.savefig(outplot_prefix+"_cohort.pdf")