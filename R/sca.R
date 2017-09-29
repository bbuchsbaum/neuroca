
#' @importFrom multiway sca
sca <- function(X, ncomp=2, center=TRUE, scale=FALSE, rank_k=NULL,
                type=c("sca-p","sca-pf2","sca-ind","sca-ecp"), ...) {
  
  assertthat::assert_that(inherits(X, "block_matrix"))
  type <- match.arg(type)
  
  X <- pre_processor(X, 
                     center=center, 
                     scale=scale)
  
  Xr <- if (!is.null(rank_k)) {
    is_reduced <- TRUE
    reducer <- reduce_rank(X, rank_k)
    reducer$x
  } else {
    is_reduced=FALSE
    reducer <- NULL
    X
  }
  
  Xl <- lapply(as.list(Xr), t)
  
  sca_fit <- multiway::sca(Xl, nfac=ncomp, type=type, ...)
  
  reprocess <- function(newdat, table_index) {
    prep <- attr(X, "pre_process")
    newdat <- prep(newdat, bind[[table_index]])
    
    if (is_reduced) {
      newdat <- project(reducer, newdat, table_index)
    }
    
    newdat
  }
  
  ret <- list(
    sca_fit=sca_fit,
    reprocess=reprocess,
    center=center,
    scale=scale,
    rank_k=rank_k,
    is_reduced=is_reduced)
  
  class(ret) <- c("sca", "list")
  ret  
}