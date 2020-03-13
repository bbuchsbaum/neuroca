#' @export
#' 
#' @example 
#' 
#' X1 <- matrix(rnorm(20*10), 20, 10)
#' X2 <- matrix(rnorm(20*20), 20, 20)
#' X3 <- matrix(rnorm(20*30), 20, 30)
#' 
#' pc1 <- pca(X1, ncomp=10)
#' pc2 <- pca(X2, ncomp=20)
#' pc3 <- pca(X3, ncomp=20)
#' 
#' fits <- list(pc1,pc2,pc3)
#' pfit <- metapca(fits, ncomp=15)
#' 
metapca <- function(fits, ncomp=2, preproc=pass(), fit_weights=rep(1,length(fits))) {
  assert_that(all(sapply(fits,function(f) inherits(f, "projector"))))
  
  X <- block_matrix(lapply(fits, function(x) project(x)))
  
  if (all(fit_weights[1] == fit_weights[2:length(fit_weights)])) {
    pres <- pca(X, ncomp=ncomp, preproc=preproc)
  } else {
    wts <- rep(fit_weights/sum(fit_weights), block_lengths(X))
    pres <- genpca(unclass(X), A=wts, ncomp=ncomp, preproc=preproc)
  }
  
  ret <- bi_projector(
    preproc=pass(),
    ncomp=length(pres$d),
    v=pres$v, 
    u=pres$u, 
    d=pres$d, 
    scores=pres$scores,
    metafit=pres,
    fits=fits,
    block_indices=block_index_list(X),
    classes=c("metapca", "pca"))
  
}

collapse.metapca <- function(x, newdata, block_ind) {
  
}

reconstruct.metapca <- function(x, newdata, comp=1:length(x$d), 
                                block_ind=seq(1,length(x$fits)), 
                                reverse_pre_process=TRUE) {
  
  
  ## reonstruct here means reconstruct the original data
  ret <- lapply(block_ind, function(i) {
    recon <- reconstruct(x$metafit, comp=comp, colind=x$block_indices[[i]], reverse_pre_process=TRUE)
    as.matrix(reconstruct(x$fits[[i]], newdata=recon, reverse_pre_process=reverse_pre_process))
  })
  
  block_matrix(ret)
  
}