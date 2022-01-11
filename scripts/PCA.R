## written Margherita Francescatto 2021/05/01
## adapted from king automatically created plots
requiredPackages = c('plyr')
for(p in requiredPackages){
  chooseCRANmirror(ind = 51)
  if(!require(p,character.only = TRUE)) install.packages(p)
  library(p,character.only = TRUE)
}

# library(plyr)

## set script to accept command line arguments and get them
args <- commandArgs(TRUE)
data_name <- args[1] 
plot_dir <- args[2] 
pca_pref <- args[3]
proj_pref <- args[4]

# ## for example, for direct script testing
# data_name = "Slo_GSAMD"
# plot_dir = "/home/marghi/projects/burlo/snparrayqc_pipeline_runs/Slo_GSAMD/00_release/plots"
# pca_pref = "/home/marghi/projects/burlo/snparrayqc_pipeline_runs/Slo_GSAMD/11_pca/Slo_GSAMD_cleaned.LD0.3_kingpca"
# proj_pref = "/home/marghi/projects/burlo/snparrayqc_pipeline_runs/Slo_GSAMD/11_pca/Slo_GSAMD_cleaned.LD0.3_kingpcaproj"

## plot simple data PCA
pc_file = paste0(pca_pref, "pc.txt")
pc <- read.table(pc_file, header=T)
png(paste0(plot_dir,"/",data_name,"_pca.png"))
plot(pc$PC1, pc$PC2, type="p", xlab="PC1", ylab="PC2", pch = 16, 
     main = paste0("Population Structure in ", data_name))
dev.off()

## plot PCA of data projection on 1000Genome Project reference
pcproj_file = paste0(proj_pref, "pc.txt")
ref_file = paste0(proj_pref, "_popref.txt")
projpc <- read.table(pcproj_file, header=T)
for_col <- read.table(ref_file, header = T)[,3]
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

