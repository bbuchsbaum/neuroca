
library(neuroim)

testpath <- "/Users/brad/rstudio/neuropls/test_data/betas/"

df1 <- read.table(paste(testpath, "/3006/mvpa_design.txt",sep=""), header=TRUE)
df1 <- subset(df1, Condition == "Encod")


mask <- loadVolume(paste(testpath, "/3006/global_mask.nii",sep=""))
mask.idx <- which(mask > 0)

X <- do.call(cbind, lapply(df1$offset, function(off) {
  print(off)
  loadVolume(as.character(paste0(testpath, "/3006/", paste0(df1$image[1]))), off)[mask.idx]
}))

X <- t(X)

Xs <- scale(X, scale=FALSE)
Xss <- scale(X)
Y <- df1$Video


fit <- pls.meancen(Y, Xss, ncomp=7, cv=FALSE, boot=FALSE, svd.method="fast")

thresh <- c(.000000001, .0000001, .0001, .001, .05, 1)

acc <- do.call(cbind, lapply(thresh, function(thresh) {
  print(thresh)
  feat <- AnovaFeatureSel(thresh=thresh, criterion="pval")
  res = cv(fit, 5, 10, featSel=feat)
  apply(res, 2, function(ypred) {
    sum(ypred ==Y)/length(Y)
  })
}))

