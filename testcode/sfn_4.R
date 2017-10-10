library(tibble)
library(dplyr)
library(neuroim)
library(assertthat)
library(devtools)
library(energy)
library(pls)
library(parallel)

devtools::load_all()

PATH <- "~/Dropbox/Brad/neuropls/video_marie/SFN/"
setwd(PATH)




sids <- scan("sids.txt", "")


clusvol <- loadVolume("nw_marie_video_1000.nii")

clus.idx <- which(clusvol > 0)

emat <- as.matrix(loadVector("encode_tstat.nii.gz"))[clus.idx,]
rmat <- as.matrix(loadVector("recall_tstat.nii.gz"))[clus.idx,]

emat <- scale(emat, center=TRUE, scale=TRUE)
rmat <- scale(rmat, center=TRUE, scale=TRUE)


red_emat <- t(splitReduce(emat, factor(clusvol[clus.idx])))
red_rmat <- t(splitReduce(rmat, factor(clusvol[clus.idx])))

pls.1 <- plsr(red_rmat ~ red_emat, ncomp=10, method="kernelpls")

df1 <- data.frame(sid=rep(sids, each=11))

perfmat <- mclapply(1:length(sids), function(i) {
    print(i)
    heldout <- which(sids[i] == df1$sid)
    keep <- which(sids[i] != df1$sid)
    
    Xtest <- red_rmat[heldout,]
    pls.x <- plsr(red_rmat[keep,] ~ red_emat[keep,], ncomp=10, method="kernelpls")
    ##pls.y <- spls(red_emat[keep,], red_rmat[keep,], K=10, kappa=.5, eta=.5)
    
    ret <- lapply(1:10, function(nc) {
      sum((red_rmat[heldout,] - predict(pls.x, newdata=red_emat[heldout,], ncomp=nc)[,,1])^2)
    })
    
    unlist(ret)
}, mc.cores=8)



