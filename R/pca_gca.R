#' @importFrom RGCCA rgcca
pca_gca <- function(X, ncomp=rep(2, length(X)), preproc=center(), cor_min=.7) {
  
  
  assertthat::assert_that(inherits(X, "block_matrix"), msg="X must be a 'block_matrix'")
  assert_that(all(ncomp > 1), "all entries in `ncomp`` must be greater than 1")
  
  svdlist <- lapply(1:nblocks(X), function(i) {
    Xi <- get_block(X,i)
    svd_i <- svd_wrapper(Xi, ncomp=ncomp[i], method="fast")
  })
  
  sclist <- lapply(svdlist, "[[", "u")
  ccorres <- RGCCA::rgcca(sclist, C=1-diag(length(svdlist)), tau = rep(0, length(svdlist)),
                                ncomp = rep(min(ncomp), length(svdlist)), verbose = FALSE)
  com_comp <- sum(sqrt(ccorres$AVE$AVE_inner) >= cor_min)
  
}
