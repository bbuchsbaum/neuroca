

subpca <- function(X, groups, method=c("gcv", "fixed"), 
                   ncomp=2, min_k=1, max_k=3, 
                   center=TRUE, scale=FALSE, ...) {
  
  
  assertthat(length(groups) == ncol(X))
  
  method <- match.arg(method)
  
  ngroups <- length(unique(groups))
  
  if (length(ncomp) == 1) {
    ncomp <- rep(ncomp, ngroups)
  } else {
    assertthat::assert_that(length(ncomp) == ngroups)
  }
  
  group_indices <- split(1:length(groups), factor(groups))
  
  Xblock <- block_matrix(lapply(group_indices, function(ind) {
    X[,ind,drop=FALSE]
  }))
  
  fits <- if (method == "gcv") {
    fits <- lapply(1:nblocks(Xblock), function(i) {
      xb <- get_block(Xblock, i)
      est <- fast_estim_ncomp(xb, ncp.min=min_k, ncp.max=max_k)
      pca(xb, ncomp=max(min_k, est$bestcomp), center=center, scale=scale)
    })
  } else if (method == "fixed") {
    ## method is fixed
    lapply(1:nblocks(Xblock), function(i) {
      xb <- get_block(Xblock, i)
      pca(xb, ncomp=ncomp[i], center=center, scale=scale)
    })
  } else {
    lapply(1:nblocks(Xblock), function(i) {
      xb <- get_block(Xblock, i)
      denoiseR::adashrink(xb)
      #pca(xb, ncomp=ncomp[i], center=center, scale=scale)
    })
  }
  
  components <- lapply(fits, function(x) scores(x))
  
  bm <- block_matrix_list(components)
  
  projector <- function(x,i) {
    project(fits[[i]], x)
  }
  
  global_projector <- function(x) {
    out <- lapply(1:ngroups, function(i) {
      xi <- x[, group_indices[[i]]]
      projector(xi,i)
    })
    
    ## need to reorder
    block_matrix_list(out)
  }
  
  ret <- list(x=bm, block_projector=projector, global_projector=global_projector)
  class(ret) <- c("subpca", "projector", "list")
  ret
  
}
  