library(neuroim)

path <- "/Users/bbuchsbaum/rstudio/neuropls/test_data/musutest/"
mask <- loadVolume(paste0(path, "group_mask.nii"))
mask.idx <- which(mask > 0)

sids <- c(1001,1002,1003,1007,1008,1009,1010)
vn <- paste0(sids, "_nw_videos_runavg.nii.gz")
vnames <- paste0(path, vn)

dn <- paste0(sids, "_design.txt")
dnames <- paste0(path, dn)
des <- do.call(rbind, lapply(dnames, read.table, header=TRUE))

Y <- paste0(des$video, ":", des$condition)
  
Xlist <- lapply(vnames, function(fname) {
  x <- scale(as.matrix(loadVector(fname))[mask.idx,], center=TRUE, scale=FALSE)
  t(x)
})

