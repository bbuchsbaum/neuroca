library(neighborweights)
library(dplyr)
library(Matrix)

base_path <- "/Users/bbuchsbaum/analysis/airplane_paper/gpca"

dat <- readRDS(paste0(base_path, "/", "betas_loc_mvideo.rds"))


des_all <- do.call(rbind, lapply(dat, function(x) {
  x$des
}))

des_enc <- des_all %>% filter(Condition == "Encod")
des_rec <- des_all %>% filter(Condition == "Recall")

Xenc <- lapply(dat, function(x) {
  idx <- which(x$des$Condition == "Encod")
  m <- x$mat[idx,]
  m <- scale(m, center=TRUE, scale=FALSE)
  m2 <- t(apply(m, 1, function(v) {
    v0 <- v-mean(v)
    v0/sqrt(sum(v0^2))
  }))
  
  m2/svd(m2)$d[1]
})

Xenc_resid <- lapply(dat, function(x) {
  idx <- which(x$des$Condition == "Encod")
  m <- x$mat[idx,]
  m2 <- t(apply(m, 1, function(v) {
    v0 <- v-mean(v)
    v0/sqrt(sum(v0^2))
  }))
  
  form = ~ Video + factor(run)
  neuroca:::residualize(form, m2, x$des[idx,])
})

Xrec <- lapply(dat, function(x) {
  idx <- which(x$des$Condition == "Recall")
  x$mat[idx,]
})

Srepnum <- graph_weights(matrix(des_enc$repnum), k=5*length(dat)*11, 
                         neighbor_mode="knn", weight_mode="heat", sigma=5)
Srepnum <- neighborweights::make_doubly_stochastic(Srepnum, iter=10)

Svideo <- binary_label_matrix(des_enc$Video, des_enc$Video)
SSubject <- binary_label_matrix(des_enc$sid, des_enc$sid)
Sboth <- Srepnum * Svideo

Sboth <- (Sboth + t(Sboth))/2
pp=RSpectra::eigs(Sboth, k=1)$values[1]
Sboth <- Sboth/pp

sids <- sapply(dat, function(x) x$des$sid[1])
old <- substr(sids, 1,1) == "2"
young <- substr(sids, 1,1) != "2"

Xencall <- do.call(Matrix::bdiag, Xenc_resid[young])
destmp = des_enc %>% 
fit <- genpca(Xencall,M=Sboth[1:nrow(Xencall),1:nrow(Xencall)], preproc=pass(), ncomp=11)


doit <- function(keep) {
  Xencall <- do.call(Matrix::bdiag, Xenc[keep])
  ids <- which(des_enc$sid %in% sids[keep])
  destmp <- des_enc[ids,]

  Smod <- (Svideo * SSubject) + Svideo/20
  Smod <- Smod/RSpectra::eigs_sym(Smod,k=1)$values[1]
  fit <- genpca(Xencall,M=Smod[ids,ids], preproc=pass(), ncomp=4)
  list(fit=fit, destmp=destmp)
}

fit2 <- doit(young | old)
qplot(fit2$fit$scores[,1], fit2$fit$scores[,2], colour=Video, data=fit2$destmp, geom=c("blank")) +
  geom_text(aes(x=fit2$fit$scores[,1], y=fit2$fit$scores[,2], label=sid),data=fit2$destmp)



