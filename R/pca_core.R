
#' @export
#' @param X
#' @param ncomp
#' @param center
#' @param scale
#' @param svd.method
pca_core <- function(X, ncomp=2, center=FALSE, scale=FALSE, svd.method="base") {
  
  if (center || scale) {
    X <- scale(X, center=center, scale=scale)
  }
  
  
  svdres <- svd.wrapper(X, ncomp, method=svd.method)
  
  scores <- svdres$v %*% diag(svdres$d, nrow=svdres$ncomp, ncol=svdres$ncomp)
  
  ret <- list(v=svdres$v, u=svdres$u, d=svdres$d, scores=scores, ncomp=ncomp, svd.method=svd.method)
  
  class(ret) <- c("pca_result")
  ret
}

#' @export
loadings.pca_result <- function(x) {
  x$v
}