#! /usr/bin/env Rscript

library(dclust)
library(neuroim2)
library(neuroca)
library(dendextend)

init_future()

mask <- read_vol("testcode/sub-1003-MNI152NLin2009cAsym_global_mask.nii")
mask.idx <- which(mask>0)
vec <- read_vec("testcode/sub-1003_ses-1_task-speech_run-01_bold_space-MNI152NLin2009cAsym_preproc.nii.gz", mask=mask)

X <- series(vec, which(mask>0))

cds <- index_to_coord(mask, mask.idx)
hclus <- dclust(cds)

res <- neuroca:::hpca(t(X), hclus, cuts=c(2,4,8,16,32), est_method="smooth",
               comp_method="log",
               cds=cds, spat_smooth=c(16,8,4,2,2))
               
