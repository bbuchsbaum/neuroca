

#' svd.wrapper
#' 
#' @param X the \code{matrix}
#' @param ncomp number of components to estimate
#' @param method the svd method to use. One of: 'base', 'fast', 'irlba', 'propack'
#' @export
svd.wrapper <- function(X, ncomp=min(dim(X)), method=c("base", "fast", "irlba", "propack")) {
  assert_that(method %in% c("base", "fast", "irlba", "propack"))
  
  res <- switch(method[1],
                base=svd(X),
                fast=corpcor:::fast.svd(X),
                propack=svd::propack.svd(X, neig=ncomp),
                irlba=irlba:::irlba(X, nu=min(ncomp, min(dim(XC)) -3), nv=min(ncomp, min(dim(XC)) -3)))
  
  
  res$d <- res$d[1:length(res$d)]
  res$u <- res$u[,1:length(res$d), drop=FALSE]
  res$v <- res$v[,1:length(res$d), drop=FALSE]
  res$ncomp <- length(res$d)
  res
}
