
#' @importFrom multiway sca
#' @export
#' @examples 
#' 
#' X <- block_matrix(replicate(10, { matrix(rnorm(10*10), 10, 10) }, simplify=FALSE))
#' res <- sca(X, ncomp=5, type="sca-pf2")
#' lds <- loadings(res)
#' bind <- block_index_list(res)
#' blds <- lapply(seq_along(bind), function(i) lds[bind[[i]],])
#' stopifnot(ncol(scores(res)) == 3)
sca <- function(X, ncomp=2, preproc=center(),
                type=c("sca-p","sca-pf2","sca-ind","sca-ecp"), ...) {
  
  assertthat::assert_that(inherits(X, "block_matrix"))
  type <- match.arg(type)
  
  
  ## pre-process the variables.
  procres <- prep(preproc, X)
  Xp <- procres$init(X)

  Xl <- lapply(as.list(Xp), t)
  
  bind <- block_index_list(X)
  
  sca_fit <- multiway::sca(Xl, nfac=ncomp, type=type, ...)
  
  Dmat <- do.call(rbind, sca_fit$D)
  vscale <- apply(Dmat, 2, function(x) sqrt(sum(x^2)))
  
  
  v <- sweep(Dmat, 2, vscale, "/")
  u <- sca_fit$B / sqrt(nrow(X))
  d <- vscale * sqrt(nrow(X))
  
  fit <- pseudo_svd(u, v, d, rnames=row.names(X))
  
  
  ret <- list(
    X=Xp,
    preproc=procres,
    sca_fit=sca_fit,
    fit=fit,
    ncomp=fit$ncomp,
    block_indices=bind,
    nvars=ncol(X),
    ntables=length(block_lengths(X)))
  
  class(ret) <- c("sca", "multiblock", "bi-projector", "list")
 
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


