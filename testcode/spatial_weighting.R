library(Matrix)
dat <- readRDS("testcode/normals_betas_14.rds")

perc_dat <- lapply(dat, function(x) {
  keep <- which(x$design$Condition == "Encod")
  m <- residualize(~ factor(Run), x$mat, x$design)
  
  des <- x$design[keep,]
  m <- m[keep,]
  
  cscore <- sda::catscore(as.matrix(m), des$Video, diagonal=FALSE)
  fscore <- apply(cscore, 1, function(x) sqrt(sum(x^2)))
  
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

wg <- unlist(lapply(perc_dat, function(x) x$fscore/sum(x$fscore)))
wg <- wg/max(wg)

folds <- lapply(perc_dat, function(x) x$design$Run)


Y <- rep(perc_dat[[1]]$design$Video, 14)

do_wmubada <- function(sigma, sigmab, bw, wg) {
  SA <- spatial_constraints(coords,nblocks=length(perc_dat), feature_scores=wg, sigma=sigma, sigma_between=sigmab, shrinkage_factor=bw)
  
  mu.1 <- mubada(Y, Xlist, normalization="custom", A=SA, ncomp=10)
  p <- performance(mu.1, metric="ACC", type="prob")

  df1 <- data.frame(p=unlist(p), sid=1:length(p), bw=bw, sigmab=sigmab, sigma=sigma)
  print(df1)
  df1
}

res <- do.call(rbind, lapply(1:nrow(grid), function(i) {
  do_wmubada(sigma=grid$sigma[i], sigmab=grid$sigmab[i], bw=grid$bw[i], wg=rep(1, length(wg)))
}))



library(dplyr)
library(ggplot2)
res2 <- res %>% group_by(bw, sigmab, sigma) %>% summarize(p=mean(p), pmed=median(p))
res3 <- res %>% group_by(sigma) %>% summarize(p=mean(p), pmed=median(p))
res4 <- res %>% group_by(sigmab) %>% summarize(p=mean(p), pmed=median(p))
res5 <- res %>% group_by(bw) %>% summarize(p=mean(p), pmed=median(p))
res6 <- res %>% group_by(bw,sigmab) %>% summarize(p=mean(p), pmed=median(p))
res7 <- res %>% group_by(bw,sigma) %>% summarize(p=mean(p), pmed=median(p))


qplot(sigma, p, colour=factor(sigmab), shape=factor(bw), data=res2)
qplot(sigma, p, data=res3) + geom_line()
qplot(sigmab, p, data=res4) + geom_line()
qplot(bw, p, data=res5) + geom_line()
qplot(bw, p, colour=factor(sigmab), data=res6) + geom_line()
qplot(bw, p, colour=factor(sigma), data=res7) + geom_line()


#mfa.1 <- mfa(Xbar, normalization="custom", A=S, ncomp=6, center=TRUE, scale=TRUE)


mu.1 <- mubada(Y, Xlist, normalization="custom", A=S, ncomp=10)
mu.2 <- mubada(Y, Xlist, normalization="MFA", ncomp=10, center=TRUE)

folds <- lapply(perc_dat, function(x) x$design$Run)

performance(mu.1, metric="ACC", type="prob")
performance(mu.2, metric="ACC", type="cosine", folds=folds)
