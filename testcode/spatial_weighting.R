library(Matrix)
dat <- readRDS("testcode/normals_betas_14.rds")

perc_dat <- lapply(dat, function(x) {
  keep <- which(x$design$Condition == "Encod")
  m <- residualize(~ factor(Run), x$mat, x$design)
  
  des <- x$design[keep,]
  m <- m[keep,]
  
  cscore <- sda::catscore(as.matrix(m), des$Video, diagonal=TRUE)
  fscore <- apply(cscore, 1, function(x) max(abs(x)))
  
  if (nrow(m) != 77) {

    m <- rbind(m, m[56:66,])
    tmp <- subset(des, Run==6)
    tmp$Run <- as.integer(tmp$Run)
    tmp$Run <- 7
    des$Run <- as.integer(des$Run)
    des <- rbind(des, tmp)
    des$Run <- factor(des$Run)
  }
  
  m[which(m == 0)] <- rnorm(length(which(m==0))) * .1
  m <- t(scale(t(m)))

  list(sid=x$sid, design=des, mat=m, idx=x$idx,cds=x$cds, fscore=fscore)
})

fmat <- do.call(cbind, lapply(perc_dat, function(x) x$fscore))
avg_fscore <- rowMeans(fmat)
coords <- dat[[1]]$cds
indices <- rep(1:nrow(coords), length(dat))

Xlist <- lapply(perc_dat, function(x) x$mat)


Xbar <- block_matrix(lapply(Xlist, function(x) {
  gm <- group_means(perc_dat[[1]]$design$Video, x)
  gm[gm==0] <- rnorm(sum(gm==0))*.1
  gm
}))
  
wg <- unlist(lapply(perc_dat, function(x) x$fscore/sum(x$fscore)))
wg <- wg/max(wg)
cds2 <- do.call(rbind, lapply(1:length(dat), function(i) cbind(coords, i*10)))

folds <- lapply(perc_dat, function(x) x$design$Run)



do_wmubada <- function(sigma, bw, bdthresh=1) {
  Swithin <- neighborweights::spatial_smoother(cds2,sigma=sigma,nnk=27,stochastic = TRUE)
  Sbetween <- neighborweights::spatial_adjacency(as.matrix(indices),weight_mode="binary", 
                                                 normalize=FALSE,dthresh=bdthresh, include_diagonal=FALSE)/length(perc_dat)
  diag(Sbetween) <- 0

  Wg <- Diagonal(x=sqrt(wg))

  S <- Wg %*% Swithin %*% Wg
  S <- S/(RSpectra::eigs_sym(S, k=1, which="LA")$values[1])
  S <- S + Sbetween*bw
  S <- S/(RSpectra::eigs_sym(S, k=1, which="LA")$values[1])
  
  mu.1 <- mubada(Y, Xlist, normalization="custom", A=S, ncomp=10)
  p <- performance(mu.1, metric="ACC", type="prob")
  
  data.frame(p=unlist(p), sid=1:length(p), bw=bw, sigma=sigma)
}

grid <- expand.grid(bw=seq(0,.1, by=.01), sigma=c(2,5,8))

res <- do.call(rbind, lapply(1:nrow(grid), function(i) {
  do_wmubada(sigma=grid$sigma[i], bw=grid$bw[i])
}))

#mfa.1 <- mfa(Xbar, normalization="custom", A=S, ncomp=6, center=TRUE, scale=TRUE)

Y <- rep(perc_dat[[1]]$design$Video, 14)
mu.1 <- mubada(Y, Xlist, normalization="custom", A=S, ncomp=10)
mu.2 <- mubada(Y, Xlist, normalization="MFA", ncomp=10, center=TRUE)

folds <- lapply(perc_dat, function(x) x$design$Run)

performance(mu.1, metric="ACC", type="prob")
performance(mu.2, metric="ACC", type="cosine", folds=folds)
