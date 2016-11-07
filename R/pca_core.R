
#' @export
#' @param X
#' @param ncomp
#' @param center
#' @param scale
#' @param svd.method
pca_core <- function(X, ncomp=min(dim(X)), center=TRUE, scale=FALSE, svd.method="base") {
  
  if (center || scale) {
    X <- scale(X, center=center, scale=scale)
  }
  
  
  svdres <- svd.wrapper(X, ncomp, method=svd.method)
  
  #scores <- t(t(svdres$u) %*% diag(svdres$d, nrow=svdres$ncomp, ncol=svdres$ncomp))
  scores <- t(t(as.matrix(svdres$u)) * svdres$d)
  ret <- list(v=svdres$v, u=svdres$u, d=svdres$d, scores=scores, ncomp=ncomp, svd.method=svd.method)
  
  class(ret) <- c("pca")
  ret
}

#' @export
loadings.pca <- function(x) {
  x$v
}

scores.pca <- function(x) {
  x$scores
}

