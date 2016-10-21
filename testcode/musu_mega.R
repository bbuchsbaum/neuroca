library(neuroim)

videos <- scan("~/Dropbox/Brad/neuropls/megalocalizer/clus1000/videos.txt", "")
Xlist <- readRDS("~/Dropbox/Brad/neuropls/megalocalizer/clus1000/all_subs_clus_1000.rds")
sXlist <- lapply(Xlist, function(x) x/svd(x)$d)

Y <- factor(videos)
sids <- names(Xlist)


#Xcat <- do.call(rbind, Xlist)
#Yb <- rep(Y, length(Xlist))
#sid <- factor(rep(sids, each=length(videos)))
#bres <- bada(Yb, Xcat, center=TRUE, strata=sid)

for (i in 1:length(Xlist)) {
  
  lds <- loadings(result,i)
  vol <- loadVolume(paste0("~/Dropbox/Brad/neuropls/megalocalizer/clus1000/nw_", 
                           sids[i], "_s3_clus_betas_1000.nii"))
  
  for (j in 1:20) {
    out <- fill(vol, cbind(1:1000, lds[,j]))
    writeVolume(out, paste0(sids[i], "_lds_", j, ".nii"))
  }                
}
