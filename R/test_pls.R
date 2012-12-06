

testpath <- "/Users/brad/rstudio/neuropls/test_data/"
library(neuroim)
df1 <- read.table(paste(testpath, "/1014/block_design.txt",sep=""), header=TRUE)

fnames <- paste(testpath, "/1014/", df1$encodeFiles, sep="")
mask <- loadVolume(paste(testpath, "/1014/global_mask.nii",sep=""))
mask.idx <- which(mask > 0)

X <- t(do.call(cbind, lapply(fnames, function(fname) loadVolume(fname)[mask.idx])))
Video <- df1$Video