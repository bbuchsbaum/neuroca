
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



project.sca <- function(x, newdata, ncomp=x$ncomp, pre_process=TRUE, 
                        table_index=1:x$ntables) {
  if (length(table_index) > 1) {
    assert_that(is.list(newdata) && length(newdata) == length(table_index))
  }
  
  if (is.vector(newdata)) {
    assert_that(length(table_index) == 1)
    newdata <- list(matrix(newdata, 1, length(newdata)))
  }
  
  if (is.matrix(newdata)) {
    assert_that(length(table_index) == 1)
    newdata <- list(newdata)
  }
  
  ## project new data-point(s)
  res <- lapply(1:length(table_index), function(i) {
    tbind <- table_index[i]
    xnewdat <- x$reprocess(newdata[[i]], tbind)
    x$ntables *  t(newdata[[i]]) %*% x$D[[tbind]] 
      
    
  })
  
  names(res) <- paste0("table_", table_index)
  res
} 

#' @export
reconstruct.sca <- function(x, ncomp=x$ncomp) {
  ret <- block_matrix(lapply(x$D, function(d) {
    t(tcrossprod(d, x$B))
  }))
}

#' @export
contributions.sca <- function(x, type=c("column", "row", "table")) {
  contr <- contributions(x$pca_fit)
  type <- match.arg(type)
  
  if (type == "table") {
    out <- do.call(cbind, lapply(1:x$ncomp, function(i) {
      sapply(1:x$ntables, function(j) {
        ind <- x$block_indices[[j]]
        sum(contr[i,ind])
      })
    }))
    
  } else if (type == "row") {
    contributions(x$pca_fit, type="row")
  } else {
    contributions(x$pca_fit, type="column")
  }
  
}

