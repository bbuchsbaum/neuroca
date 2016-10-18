library(neuroim)


mask <- loadVolume("~/Dropbox/Brad/neuropls/weeks/group_mask.nii")
mask.idx <- which(mask > 0)

bvec <- loadVector("~/Dropbox/Brad/neuropls/weeks/s5_pcue_all.nii.gz")
mat <- as.matrix(bvec)[mask.idx,]
time <- seq(0,24,by=2)


N <- 39
dframe <- data.frame(time=rep(time, N), sid=factor(rep(1:N, each=13)), 
                     group=c(rep("old", each=19*13), rep("young", each=20*13)))

Xlist <- lapply(1:N, function(i) {
  idx <- which(dframe$sid == i)
  t(mat[,idx])
})

Y <- ordered(as.numeric(dframe$time))
Y <- Y[1:13]

mres <- musubada(Y, Xlist,ncomp=12, normalization="None", svd.method="fast")

Xcat <- do.call(rbind, Xlist)

bres <- bada(Y=ordered(dframe$time), Xcat, ncomp=11, strata=dframe$sid)
p <- permutation(bres, ncomp=10, nperm=50)

res <- lapply(1:50, function(i) {
  print(i)
  k1 <- sort(sample(1:39, 19))
  k2 <- sort(seq(1,39)[-k1])
  mres1 <- bada(ordered(dframe$time[dframe$sid %in% k1]), do.call(rbind, Xlist[k1]),ncomp=12, svd.method="fast", strata=factor(dframe$sid[dframe$sid %in% k1]))
  mres2 <- bada(ordered(dframe$time[dframe$sid %in% k2]), do.call(rbind, Xlist[k2]),ncomp=12, svd.method="fast", strata=factor(dframe$sid[dframe$sid %in% k2]))
  sapply(1:12, function(j) {
    s1 <- mres1$scores[,1:j, drop=FALSE]
    s2 <- mres2$scores
  
    D1 <- dist(s1)
    D2 <- dist(s2)
  
    cor(as.vector(D1), as.vector(D2))
  }) 
})


m2 <- mres$permute_refit(Xlist[1:39])

tscore <- function(k) {
  lmat <- do.call(cbind, lapply(1:39, function(i) {
    loadings(mres, i)[,k]
  }))
  
  M <- rowMeans(lmat)
  se <- apply(lmat,1,sd)/sqrt(ncol(lmat))
  
  t=M/se
  
}

tscore2 <- function(k) {
  
  idx1 <- 1:19
  idx2 <- 20:39
  
  lmat <- do.call(cbind, lapply(1:39, function(i) {
    loadings(mres, i)[,k]
  }))
  
  M1 <- rowMeans(lmat[, idx1])
  M2 <- rowMeans(lmat[, idx2])
  
  se <- apply(lmat,1,sd)/sqrt(ncol(lmat)-1)
  (M1 - M2)/se
  
}

tval1 <- tscore(1)
tval2 <- tscore(2)
tval3 <- tscore(3)

bv1 <- BrainVolume(tval1, space(mask), indices=mask.idx)
bv2 <- BrainVolume(tval2, space(mask), indices=mask.idx)
bv3 <- BrainVolume(tval3, space(mask), indices=mask.idx)

writeVolume(bv1, "pc1.nii")
writeVolume(bv2, "pc2.nii")
writeVolume(bv3, "pc3.nii")

bv1 <- BrainVolume(tscore2(1), space(mask), indices=mask.idx)
bv2 <- BrainVolume(tscore2(2), space(mask), indices=mask.idx)
bv3 <- BrainVolume(tscore2(3), space(mask), indices=mask.idx)

writeVolume(bv1, "diff_pc1.nii")
writeVolume(bv2, "diff_pc2.nii")
writeVolume(bv3, "diff_pc3.nii")

Xavg <- Reduce("+", Xlist)/39
pres <- prcomp(Xavg)

library(matrixStats)

do_perm <- function(N=20, nmean=0, altmean=2, prop=.5) {
  #print("running")
  #vals=rgamma(N,3,2) + runif(1) * sign(rnorm(1))
  vals = c(rnorm((1-prop)*N, mean=nmean), rnorm(prop*N, mean=altmean))
  tval <- t.test(vals)$statistic
  ranmat = replicate(1000, { vals * sign(rnorm(length(vals))) } )
  
  M <- colMeans(ranmat)
  sds <- colSds(ranmat)
  ptvals <- M/ (sds/ (sqrt(N)))
  
  c(tval=tval, min=min(ptvals), max=max(ptvals), sd=sd(vals), sdran=sd(ptvals), npos=sum(ptvals > 1.72), nneg=sum(ptvals < -1.72))
}
  
  
ret1 <- t(replicate(1000, do_perm()))
ret2 <- t(replicate(1000, do_perm(nmean=1, altmean=1)))
ret3 <- t(replicate(1000, do_perm(nmean=0, altmean=3, prop=.3)))
ret4 <- t(replicate(1000, do_perm(nmean=0, altmean=3, prop=.7)))

grid <- expand.grid(prop = seq(0, 1, by=.1), altmean=c(1,3))

res <- lapply(1:nrow(grid), function(i) {
  print(i)
  t(replicate(1000, do_perm(20, nmean=0, altmean=grid$altmean[i], prop=grid$prop[i])))
})
  
cvals <- unlist(lapply(res, function(x) {
  cor(x[,1], x[,3])
}))







