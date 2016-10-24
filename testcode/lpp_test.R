library(neuroim)
library(foreach)
library(FNN)

mask <- loadVolume("~/Dropbox/Brad/neuropls/video_marie/1001/mean1.nii")
gmask <- loadVolume("~/Dropbox/Brad/neuropls/video_marie/s5_young_gray.nii")
des <- read.table("~/Dropbox/Brad/neuropls/video_marie/1001/trial_data_all_blockavg.txt", header=TRUE)
des1 <- subset(des, run ==1 & condition == "Encod")
idx <- which(des$run == 1 & des$condition == "Encod")

mask.idx <- which(mask !=0 & gmask > .2)
mask2 <- LogicalBrainVolume(rep(TRUE, length(mask.idx)), space(mask), indices=mask.idx)

fnames <- list.files("~/Dropbox/Brad/neuropls/video_marie/1001/", "nwtrial.*nii.gz", full.names=TRUE)

mlist <- lapply(fnames, function(fn) {
  as.matrix(loadVector(fn))[mask.idx,idx]
})

h=.73
X <- Reduce("+", mlist)/length(mlist)
Xs <- scale(t(X), scale=FALSE, center=TRUE)
Xvec <- SparseBrainVector(X, addDim(space(mask), ncol(X)), mask=mask2)
grid <- indexToGrid(mask2, mask.idx)

K <- 50
res <- lapply(1:nrow(grid), function(i) {
  print(i)
  cen <- grid[i,,drop=FALSE]
  nn <- get.knnx(grid, cen, 27)
  keep <- which(nn$nn.dist < 3 & nn$nn.dist > 0)
  
  ind <- nn$nn.index[keep]
  ovox <- grid[ind,]
  
 
  #cvals <- as.vector(cor(series(Xvec, cen), m))
  #S <- exp(cvals/(h^2))
  
  D <- nn$nn.dist[keep]
  wts <- exp(-D/(1.5^2))
  #ord <- order(S, decreasing=TRUE)
  #cbind(i=mask.idx[i], j=mask.idx[ind[ord[1:K]]], S[ord[1:K]])
  cbind(i=mask.idx[i], j=mask.idx[keep], wts)
  
})
res <- res[-7151]
smat <- do.call(rbind, res)
adj <- sparseMatrix(i=lookup(Xvec, smat[,1]), j=lookup(Xvec, smat[,2]), x=smat[,3], dims=c(length(mask.idx), length(mask.idx)))
adj <- (adj + t(adj))/2

g <- graph_from_adjacency_matrix(adj, mode="undirected", weighted=TRUE)
#D <- rowSums(adj)  
#D <- sparseMatrix(i=seq(1,length(mask.idx)), j=seq(1,length(mask.idx)), x=D)
#L <- D - adj

#Lnorm <- D_inv %*% adj %*% D_inv
Lnorm <- laplacian_matrix(g, normalized=FALSE)
Lnorm <- Lnorm/eigs(Lnorm, k=1)$values[1]
Snorm <- adj/eigs(adj, k=1)$values[1]

gres1 = sGPCA::gpca(Xs, diag(nrow(Xs)), Lnorm, K=8)
gres2 = sGPCA::gpca(Xs, diag(nrow(Xs)), Snorm, K=8)
sres <- svd(Xs)
savepc <- function(vals, oname) {
  bv <- BrainVolume(scale(vals), space(mask), indices=mask.idx)
  writeVolume(bv, oname)
}

for (i in 1:8) {
  savepc(gres1$V[,1], paste0("~/Dropbox/Brad/neuropls/video_marie/1001/gpc_laplacian_", i, ".nii"))
  savepc(gres2$V[,1], paste0("~/Dropbox/Brad/neuropls/video_marie/1001/gpc_adj_", i, ".nii"))
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

  
  