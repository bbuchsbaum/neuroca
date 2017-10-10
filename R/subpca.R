

subpca <- function(X, groups, method=c("adaptive", "fixed"), ncomp=2, min_k=1, max_k=3, center=TRUE, scale=FALSE, ...) {
  assertthat(length(groups) == ncol(X))
  
  method <- match.arg(method)
  
  ngroups <- length(unique(groups))
  
  if (length(ncomp) == 1) {
    ncomp <- rep(ncomp, ngroups)
  } else {
    assertthat::assert_that(length(ncomp) == ngroups)
  }
  
  group_indices <- split(1:length(groups), factor(groups))
  
  X <- pre_processor(X, scale, center)
  Xblock <- block_matrix(lapply(group_indices, function(ind) {
    X[,ind,drop=FALSE]
  }))
  
  if (method == "adaptive") {
    fits <- lapply(1:nblocks(Xblock), function(i) {
      xb <- get_block(Xblock, i)
      est <- fast_estim_ncomp(xb, ncp.min=min_k, ncp.max=max_k)
      pca(xb, ncomp=max(min_k, est$bestcomp), center=FALSE, scale=FALSE)
    })
  } else {
    ## method is fixed
    fits <- lapply(1:nblocks(Xblock), function(i) {
      xb <- get_block(Xblock, i)
      pca(xb, ncomp=ncomp[i], center=FALSE, scale=FALSE)
    })
  }
  
}
  