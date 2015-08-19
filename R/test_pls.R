
library(neuroim)

testpath <- "/Users/brad/rstudio/neuropls/test_data/betas/"
mask <- loadVolume(paste(testpath, "/3006/global_mask.nii",sep=""))
mask.idx <- which(mask > 0)

df1 <- read.table(paste(testpath, "/3006/mvpa_design.txt",sep=""), header=TRUE)

loadData <- function(cond) {
  df1 <- subset(df1, Condition == cond)
  X <- do.call(cbind, lapply(df1$offset, function(off) {
    print(off)
    loadVolume(as.character(paste0(testpath, "/3006/", paste0(df1$image[1]))), off)[mask.idx]
  }))
  X <- t(X)


  Y <- df1$Video
  list(X=X, Y=Y)
}

Encod <- loadData("Encod")
Recall <- loadData("Recall")

df2 <- subset(df1, Condition == "Encod")
for (i in 1:7) {
  Xtrain <- Encod$X[df2$run != i,]
  Ytrain <- Encod$Y[df2$run != i]
  
  xtest <- Encod$X[df2$run == i,]
  ytest <- Encod$Y[df2$run == i] 
  
  sda.1 <- sda(Xtrain, Ytrain)
  P <- predict(sda.1, xtest)
}


fit <- pls.meancen(Encod$Y, Encod$X, ncomp=10, cv=FALSE, boot=FALSE, svd.method="fast")
pred <- predict(fit, Recall$X, 10)

thresh <- c(.05, .1, .2, .5, .9)

acc <- do.call(cbind, lapply(thresh, function(thresh) {
  print(thresh)
  feat <- AnovaFeatureSel(thresh=.2)
  res = cv(fit, 5, 10, featSel=feat, metric="distanceRank")
  apply(res, 2, function(ypred) {
    sum(ypred ==Y)/length(Y)
  })
}))


