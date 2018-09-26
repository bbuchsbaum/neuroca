
#' @importFrom multiway sca
#' @export
sca <- function(X, ncomp=2, center=TRUE, scale=FALSE,
                type=c("sca-p","sca-pf2","sca-ind","sca-ecp"), ...) {
  
  assertthat::assert_that(inherits(X, "block_matrix"))
  type <- match.arg(type)
  
  preproc <- pre_processor(X, 
                     center=center, 
                     scale=scale)
  
  Xr <- pre_process(preproc, X)
  
  Xl <- lapply(as.list(Xr), t)
  
  bind <- block_index_list(X)
  
  sca_fit <- multiway::sca(Xl, nfac=ncomp, type=type, ...)
  
  Dmat <- do.call(rbind, sca_fit$D)
  vscale <- apply(Dmat, 2, function(x) sqrt(sum(x^2)))
  
  
  v <- sweep(Dmat, 2, vscale, "/")
  u <- sca_fit$B / sqrt(nrow(X))
  d <- vscale * sqrt(nrow(X))
  
  fit <- pseudo_svd(u, v, d, rnames=row.names(X))
  
  ret <- list(
    X=X,
    Xr=Xr,
    sca_fit=sca_fit,
    preproc=preproc,
    fit=fit,
    center=center,
    scale=scale,
    ncomp=length(d),
    block_indices=bind,
    nvars=ncol(X),
    ntables=length(block_lengths(X)))
  
  class(ret) <- c("sca", "multiblock", "list")
  ret  
}


#' @export
singular_values.sca <- function(x) x$fit$d

#' @export
scores.sca <- function(x) scores(x$fit) 

#' @export
loadings.sca <- function(x) loadings(x$fit) 

#' @export 
block_index_list.sca <- function(x) x$block_indices

#' @export
ncomp.sca <- function(x) x$ncomp

ncol.sca <- function(x) ncol(x$fit)

nrow.sca <- function(x) nrow(x$fit)

#' @export
reconstruct.sca <- function(x, ncomp=x$ncomp) {
  ret <- block_matrix(lapply(x$D, function(d) {
    t(tcrossprod(d, x$B))
  }))
}


