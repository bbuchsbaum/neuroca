

#' svd_wrapper
#' 
#' @param X the \code{matrix}
#' @param ncomp number of components to estimate
#' @param method the svd method to use. One of: 'base', 'fast', 'irlba', 'propack'
#' @export
svd_wrapper <- function(X, ncomp=min(dim(X)), method=c("base", "fast", "irlba", "propack", "rsvd", "svds"), ...) {
  method <- match.arg(method)

  res <- switch(method,
                base=svd(X,...),
                fast=corpcor:::fast.svd(X,...),
                rsvd=rsvd::rsvd(X, k=ncomp, q=2, ...),
                svds=RSpectra::svds(X,k=ncomp),
                propack=svd::propack.svd(X, neig=ncomp,...),
                irlba=irlba:::irlba(X, nu=min(ncomp, min(dim(X)) -3), nv=min(ncomp, min(dim(X)) -3)), ...)
  
 
  
  ncomp <- min(ncomp,length(res$d))
  res$d <- res$d[1:ncomp]
  res$u <- res$u[,1:ncomp, drop=FALSE]
  res$v <- res$v[,1:ncomp, drop=FALSE]
  res$ncomp <- length(res$d)
  class(res) <- c("svd", "projector", "list")
  res
}



project.svd <- function(obj, newX) {
  newX %*% obj$v
}

# compose <- function(obj1, obj2) {
#   f <- function(newX) {
#     project(obj2, project(obj1, newX))
#   }
#   
#   ret <- list(f=f)
#   class(ret) <- c("composition", "list")
# }
# 
# project.compose <- function(obj, newX) {
#   obj$f(newX)
# }


