library(neuroim)
library(foreach)
library(FNN)
library(Matrix)

dset <- readRDS("/Users/brad/code/neuropls/test_data/video/video_testset.rds")
grid <- indexToGrid(dset$mask, dset$mask.idx)

design <- dset$deslist[[1]]
enc_idx <- which(design$Condition == "Encod")

X <- dset$bvecs[[1]]@data[enc_idx,]

X <- t(scale(t(X)))

labels <- design$Video[enc_idx]

A <- spatial_adjacency(grid)
L <- spatial_laplacian(grid)
diag(A) <- 1

X <- scale(X, center=TRUE, scale=FALSE)
Lnorm <- L/RSpectra::eigs(L, k=1)$values[1]
Anorm <- A/RSpectra::eigs(A, k=1)$values[1]

Ta <- temporal_adjacency(1:nrow(X), window=4, weight_mode="heat", sigma=2)
diag(Ta) <- rowSums(Ta)
D <- Diagonal(x=rowSums(Ta))
Tl <- D - Ta
#diag(Tl) <- rowSums(Tl)

Tlnorm <- Tl/RSpectra::eigs(Tl, k=1)$values[1]
Tanorm <- Ta/RSpectra::eigs(Ta, k=1)$values[1]

gres1 = genpca(X, Anorm, Tlnorm, ncomp=10)
gres2 = genpca(X, Anorm, Tanorm, ncomp=10)


Q <- construct_weight_matrix(X, neighbor_mode="supervised", weight_mode="normalized", sigma=.4, k=5, labels=labels)
diag(Q) <- rowSums(Q)

Qnorm <- Q/RSpectra::eigs(Q, k=1)$values[1]

pgrid = expand.grid(sigma=seq(.1,.3, by=.05), k=seq(3, 20, by=4))

res <- lapply(1:nrow(pgrid), function(i) {
  print(i)
  sig <- pgrid$sigma[i]
  k <- pgrid$k[i]
  Q <- construct_weight_matrix(X, neighbor_mode="supervised", weight_mode="normalized", sigma=sig, k=k, labels=labels)
  diag(Q) <- rowSums(Q)
  Qnorm <- Q/RSpectra::eigs(Q, k=1)$values[1]
  gres1 = sGPCA::gpca(X, Qnorm, Anorm, K=10)
  fstats <- sapply(1:10, function(j) {
    summary(aov(gres1$U[,j] ~ labels))[[1]]$F[1]
  })
  
  print(fstats)
  fstats
  
})

gres1 = sGPCA::gpca(X, Qnorm, Anorm, K=6)
gres2 = sGPCA::gpca(X, diag(nrow(X)), Anorm, K=6)
gres3 = sGPCA::gpca(X, Qnorm, Lnorm, K=6)
gres4 <- svd(X)



savepc <- function(vals, oname) {
  bv <- BrainVolume(scale(vals), space(dset$mask), indices=dset$mask.idx)
  writeVolume(bv, oname)
}

for (i in 1:5) {
  savepc(gres1$v[,i], paste0("/Users/brad/code/neuropls/test_data/video/super_gpc_laptime_", i, ".nii"))
  savepc(gres2$v[,i], paste0("/Users/brad/code/neuropls/test_data/video/super_gpc_smoothtime_", i, ".nii"))
}

for (i in 1:5) {
  savepc(gres1$V[,1], "~/Dropbox/Brad/neuropls/video_marie/1001/gpc1.nii")
}
  

savepc(gres$V[,1], "~/Dropbox/Brad/neuropls/video_marie/1001/gpc1.nii")
savepc(gres$V[,2], "~/Dropbox/Brad/neuropls/video_marie/1001/gpc2.nii")
savepc(gres$V[,3], "~/Dropbox/Brad/neuropls/video_marie/1001/gpc3.nii")
savepc(gres$V[,4], "~/Dropbox/Brad/neuropls/video_marie/1001/gpc4.nii")
savepc(gres$V[,5], "~/Dropbox/Brad/neuropls/video_marie/1001/gpc5.nii")
savepc(gres$V[,6], "~/Dropbox/Brad/neuropls/video_marie/1001/gpc6.nii")
savepc(gres$V[,7], "~/Dropbox/Brad/neuropls/video_marie/1001/gpc7.nii")
savepc(gres$V[,8], "~/Dropbox/Brad/neuropls/video_marie/1001/gpc8.nii")


sres <- svd(Xs)
savepc(sres$v[,1], "~/Dropbox/Brad/neuropls/video_marie/1001/svdpc1.nii")
savepc(sres$v[,2], "~/Dropbox/Brad/neuropls/video_marie/1001/svdpc2.nii")
savepc(sres$v[,3], "~/Dropbox/Brad/neuropls/video_marie/1001/svdpc3.nii")
savepc(sres$v[,4], "~/Dropbox/Brad/neuropls/video_marie/1001/svdpc4.nii")
savepc(sres$v[,5], "~/Dropbox/Brad/neuropls/video_marie/1001/svdpc5.nii")
savepc(sres$v[,6], "~/Dropbox/Brad/neuropls/video_marie/1001/svdpc6.nii")
savepc(sres$v[,7], "~/Dropbox/Brad/neuropls/video_marie/1001/svdpc7.nii")
savepc(sres$v[,8], "~/Dropbox/Brad/neuropls/video_marie/1001/svdpc8.nii")

Q <- Exp.cov(grid, theta=5)
gres = sGPCA::gpca(Xs, diag(nrow(Xs)), Lnorm, K=8)

  
  