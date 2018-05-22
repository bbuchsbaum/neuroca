

#' block_pca
#' 
#' @param X
#' @param groups
#' @param method
#' @param ncomp
#' @param min_k
#' @param max_k
#' @param center
#' @param scale
#' @importFrom parallel mclapply
block_pca <- function(X, est_method=c("gcv", "shrink", "fixed"), 
                   ncomp=2, min_k=1, max_k=3, 
                   center=TRUE, scale=FALSE, ncores=max(parallel::detectCores()/2,1), ...) {
  
  
  assert_that(inherits(X, "block_matrix"))
  
  est_method <- match.arg(est_method)
  
  ngroups <- nblocks(X)
  
  if (length(ncomp) == 1) {
    ncomp <- rep(ncomp, ngroups)
  } else {
    assertthat::assert_that(length(ncomp) == ngroups)
  }

  
  fits <- if (est_method == "gcv") {
    fits <- parallel::mclapply(1:nblocks(X), function(i) {
      xb <- get_block(X, i)
      est <- fast_estim_ncomp(xb, ncp.min=min_k, ncp.max=max_k)
      pca(xb, ncomp=max(min_k, est$bestcomp), center=center, scale=scale)
    }, mc.cores=ncores)
  } else if (est_method == "fixed") {
    ## method is fixed
    parallel::mclapply(1:nblocks(X), function(i) {
      xb <- get_block(X, i)
      pca(xb, ncomp=ncomp[i], center=center, scale=scale)
    },mc.cores=ncores)
  } else {
    parallel::mclapply(1:nblocks(X), function(i) {
      xb <- get_block(X, i)
      shrink_pca(xb, center=center, scale=scale, ...)
    },mc.cores=ncores)
  }
  
  bm <- block_projector(fits)

}
  