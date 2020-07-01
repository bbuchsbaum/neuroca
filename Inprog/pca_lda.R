
#' @keywords internal
pooled_scatter <- function(X, Y) {
  ina <- as.integer(Y)
  s <- crossprod(X)
  ni <- sqrt(tabulate(ina))
  mi <- rowsum(X, ina)/ni
  k <- length(ni)
  denom <- dim(X)[1] - k
  for (i in 1:k) s <- s - tcrossprod(mi[i, ])
  s
}

#' @keywords internal
total_scatter <- function(X, mu) {
  p <- ncol(X)
  St <- array(0, c(p,p))
  for (i in 1:nrow(X)) {
    delta <- X[i,] - mu
    St <- St + outer(delta,delta)
  }
  
  St
}

#' @keywords internal
between_class_scatter <- function(X, Y, mu) {
  p <- ncol(X)
  levs <- levels(Y)
  
  gmeans <- group_means(Y,X)
  gmeans <- sweep(gmeans, 2, mu, "-")
  
  n <- tabulate(Y)
  
  res <- lapply(seq_along(levs), function(i) {
    n[i] * tcrossprod(gmeans[i,], gmeans[i,])
  })
  
  Reduce("+", res)
  
}

#' @keywords internal
within_class_scatter <- function(X, Y) {
  pooled_scatter(X,Y)
}


## A Unified Framework for
## Subspace Face Recognition
## Xiaogang Wang, Student Member, IEEE, and
## Xiaoou Tang, Senior Member, IEEE
#' @keywords internal
pca_lda <- function(Y, X, preproc=center(), dp=min(dim(X)), di, dl=dp-1) {
  Y <- as.factor(Y)
  
  procres <- prep(preproc, X)
  Xp <- procres$Xp
  
  pca_basis <- pca(Xp, ncomp=dp)
  
  proj_dp <- loadings(pca_basis)
  
  Xpca <- scores(pca_basis)
  mu <- colMeans(Xpca)
  
  Sw <- within_class_scatter(Xpca, Y)
  Sb <- between_class_scatter(Xpca, Y, mu)
  
  gmeans <- group_means(Y, Xpca)
  
  E_i <- RSpectra::eigs_sym(Sw, k=dl)
  
  proj_di <- E_i$vectors %*% diag(1/sqrt(E_i$values))
  
  gmeans_proj <- (gmeans %*% proj_di)
  
  E_l <- pca(gmeans_proj, ncomp=dl)
  
  proj_dl <- loadings(E_l)
  
  proj_final <- proj_dp %*% proj_di %*% proj_dl
  
  projector(procres, ncomp=ncol(proj_final), v=proj_final, classes="pca_lda")
}