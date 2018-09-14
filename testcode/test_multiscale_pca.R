library(neuroim2)
library(fastcluster)
library(neurocluster)
mask <- read_vol("testdata/global_mask.nii")
mask2 <- read_vol("testdata/s4_episeg_1.nii")
mask3 <- read_vol("testdata/s4_episeg_2.nii")
vec <- read_vec("testdata/betas_1.nii")

mask.idx <- which(mask>0 & mask2>.1)

maskc <- mask
maskc[mask.idx] <- 1
maskc[-mask.idx] <- 0

clus <- commute_cluster(vec, mask=maskc, K=350, weight_mode="heat", alpha=.5)
cres <- meta_clust(clus, cuts = c(16,32,64,128,256))