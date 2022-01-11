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

#define a function to generate plots we need
def PCAplots(wildcards,plot_dir, pca_pref, proj_pref):
     #get args
     data_name = wildcards.vcf_name
     # plot_dir = plot_dir
     pca_file = pca_pref + "pc.txt"
     # pca_file = "/large/___SCRATCH___/burlo/cocca/WGS_JOINT_CALL/WGS_QC_pre_release/20220110/03.samples/PCA/WGS_ITA_PREREL_MERGED_cleaned.LD0.3_kingpcapc.txt"
     outplot = os.path.join(plot_dir,data_name + "_pca.png")
     
     ## plot simple data PCA
     pca_df = pd.read_table(pca_file,sep=" ", header=0)
     pca_df.plot.scatter(x='PC1',y='PC2')
     plt.title("Population structure")
     # plt.savefig('test.pdf')
     plt.savefig(outplot)

     ## plot PCA of data projection on 1000Genome Project reference
     pcproj_file = proj_pref + "pc.txt"
     # pcproj_file = "/large/___SCRATCH___/burlo/cocca/WGS_JOINT_CALL/WGS_QC_pre_release/20220110/03.samples/PCA/WGS_ITA_PREREL_MERGED_cleaned.LD0.3_kingpcaprojpc.txt"
     ref_file = proj_pref + "_popref.txt"
     # ref_file = "/large/___SCRATCH___/burlo/cocca/WGS_JOINT_CALL/WGS_QC_pre_release/20220110/03.samples/PCA/WGS_ITA_PREREL_MERGED_cleaned.LD0.3_kingpcaproj_popref.txt"

     projpc_df = pd.read_table(pcproj_file,sep=" ", header=0)
     for_col_df = pd.read_table(ref_file, sep=" ", header=0).iloc[:,2]

     ## I want to assign a color to each population, I use mapvalues for this
     from <- names(summary(as.factor(for_col))) ## this gives me the four pops
     ## I pick 5 colors, I choose grey for EUR since I will in most cases project over this
     to <- c("orange","dodgerblue","orchid","grey","black") 
     col <- mapvalues(for_col, from, to)
     ## now I plot all points using the colors I just picked
     png(paste0(plot_dir,"/",data_name,"_pca_projection_on_1000GP.png"))
     plot(projpc$PC1, projpc$PC2, type="p", xlab="PC1", ylab="PC2",
          main = paste0("Population Structure in ", data_name), col = col)
     ## I add points representing my samples in red (because they are EUR they show as red on grey)
     points(projpc$PC1[projpc$AFF==2], projpc$PC2[projpc$AFF==2], col = "red")
     ## add legend
     legend("topright", c(from, data_name),
            col=c(to,"red"),text.col = c(to, "red"), pch = 19, cex = 0.9)
     dev.off()

