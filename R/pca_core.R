

#' @param X
#' @param ncomp
#' @param center
#' @param scale
#' @param svd.method
#' @export
pca_core <- function(X, ncomp=min(dim(X)), center=TRUE, scale=FALSE, svd.method="base") {
  
  X <- pre_processor(X, center,scale)
  
  svdres <- svd_wrapper(X, ncomp, method=svd.method)
  
  scores <- t(t(as.matrix(svdres$u)) * svdres$d)
  
  ret <- list(v=svdres$v, u=svdres$u, d=svdres$d, scores=scores, ncomp=ncomp, svd.method=svd.method, pre_process=attr(X, "pre_process"))
  
  class(ret) <- c("pca", "projector")
  ret
}

#' @export
loadings.pca <- function(x) {
  x$v
}

#' @export
scores.pca <- function(x) {
  x$scores
}

#' @export
project.pca <- function(obj, newX) {
  obj$pre_process(newX) %*% obj$v
}


hierarchical_pca <- function(X, comps, hierarchy, svd.method="base") {
}