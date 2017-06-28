

#' @param X
#' @param ncomp
#' @param center
#' @param scale
#' @param svd.method
#' @export
pca_core <- function(X, ncomp=min(dim(X)), center=TRUE, scale=FALSE, svd.method="fast") {
  
  X <- pre_processor(X, center=center,scale=scale)
  
  svdres <- svd_wrapper(X, ncomp, method=svd.method)
  
  scores <- t(t(as.matrix(svdres$u)) * svdres$d)
  
  ret <- list(v=svdres$v, u=svdres$u, d=svdres$d, scores=scores, ncomp=ncomp, svd.method=svd.method, pre_process=attr(X, "pre_process"), reverse_pre_process=attr(X, "reverse"))
  
  class(ret) <- c("pca", "projector", "list")
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

project_cols <- function(obj, newX) {
  newX %*% obj$v
}

#' @export
reconstruct.pca <- function(x, ncomp=x$ncomp) {
  x$reverse_pre_process(x$scores[,1:ncomp,drop=FALSE] %*% t(x$v[,1:ncomp,drop=FALSE]))
}


truncate.pca <- function(obj, ncomp) {
  if (ncomp >= obj$ncomp) {
    warning("number of components to keep is greater than or equal to rank of pca fit, returning original model")
    ret <- obj
  } else {
    ret <- list(v=obj$v[,1:ncomp], u=obj$u[,1:ncomp], d=obj$d[1:ncomp], scores=obj$scores[,1:ncomp], ncomp=ncomp, svd.method=obj$svd.method, pre_process=obj$pre_process)
    class(ret) <- c("pca", "projector", "list")
  }
  
  ret
}


reduce_rank.matrix <- function(x, k=min(dim(x)), center=TRUE, scale=FALSE, reducer=pca_core, ...) {
  res <- reducer(x, k, center=center, scale=scale, ...)
  
  #ret <- list(
}


# 
# hierarchical_pca <- function(X, comps, hierarchy, svd.method="base") {
# }