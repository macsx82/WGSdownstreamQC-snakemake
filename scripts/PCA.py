## Python port of the PCA.R script written Margherita Francescatto 2021/05/01
## adapted from king automatically created plots
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

#define a function to generate plots we need
def PCAplots(plot_dir, pca_in, proj_in, proj_ref,vcf_name):
     #try to fix X11 error
     matplotlib.use('Agg')
     #get args
     data_name = vcf_name
     # plot_dir = plot_dir
     pca_file = pca_in
     # pca_file = "/large/___SCRATCH___/burlo/cocca/WGS_JOINT_CALL/WGS_QC_pre_release/20220110/03.samples/PCA/WGS_ITA_PREREL_MERGED_cleaned.LD0.3_kingpcapc.txt"
     outplot = os.path.join(plot_dir,data_name + "_pca.png")
     outplot_pj = os.path.join(plot_dir,data_name + "_pca_projection_on_1000GP.png")
     
     ## plot simple data PCA
     pca_df = pd.read_table(pca_file,sep=" ", header=0)
     pca_df.plot.scatter(x='PC1',y='PC2')
     plt.title("Population structure")
     # plt.savefig('test.pdf')
     plt.savefig(outplot)

     ## plot PCA of data projection on 1000Genome Project reference
     pcproj_file = proj_in
     # pcproj_file = "/large/___SCRATCH___/burlo/cocca/WGS_JOINT_CALL/WGS_QC_pre_release/20220110/03.samples/PCA/WGS_ITA_PREREL_MERGED_cleaned.LD0.3_kingpcaprojpc.txt"
     ref_file = proj_ref
     # ref_file = "/large/___SCRATCH___/burlo/cocca/WGS_JOINT_CALL/WGS_QC_pre_release/20220110/03.samples/PCA/WGS_ITA_PREREL_MERGED_cleaned.LD0.3_kingpcaproj_popref.txt"
     projpc_df = pd.read_table(pcproj_file,sep=" ", header=0)
     for_col_df = pd.read_table(ref_file, sep=" ", header=0)

     merged_df = projpc_df.merge(for_col_df, how='inner',on='IID')

     ## I want to assign a color to each population, I use mapvalues for this
     pops = list((for_col_df['Population']).unique()) ## this gives me the pops in the reference file
     pops.sort()
     ## I pick 5 colors, I choose grey for EUR since I will in most cases project over this
     colors = ["orange","dodgerblue","orchid","grey","black"]
     colors_map = dict(zip(pops,colors))
     #group by population on the merged dataset
     merged_df_grouped= merged_df.groupby('Population')
     #initialize the plot
     fig, ax = plt.subplots()
     ## plot projection using the colors just picked
     for key, group in merged_df_grouped:
          group.plot(ax=ax,kind='scatter',x='PC1',y='PC2', color=colors_map[key], label=key)
     #add the study population data
     projpc_df[projpc_df['AFF']==2].plot(ax=ax,kind='scatter',x='PC1',y='PC2', color='red', label='Study pop')

     plt.title("Population structure")
     # plt.savefig('test2.pdf')
     plt.savefig(outplot_pj)
