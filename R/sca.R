
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
  bind <- block_index_list(X)
  
  sca_fit <- multiway::sca(Xl, nfac=ncomp, type=type, ...)
  
  Dmat <- do.call(rbind, sca_fit$D)
  
  vscale <- apply(Dmat, 2, function(x) sqrt(sum(x^2)))
  
  
  v <- sweep(Dmat, 2, vscale, "/")
  u <- sca_fit$B / sqrt(nrow(X))
  d <- vscale * sqrt(nrow(X))
  
  reprocess <- function(newdat, table_index) {
    prep <- attr(X, "pre_process")
    newdat <- prep(newdat, bind[[table_index]])
    
    if (is_reduced) {
      newdat <- project(reducer, newdat, table_index)
    }
    
    newdat
  }
  
  fit <- pseudo_pca(u, v, d)
  
  ret <- list(
    fit=sca_fit,
    projector=fit,
    reprocess=reprocess,
    center=center,
    scale=scale,
    ncomp=length(d),
    rank_k=rank_k,
    block_indices=bind,
    ntables=length(block_lengths(X)),
    is_reduced=is_reduced)
  
  class(ret) <- c("sca", "multiblock", "list")
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
    ind <- x$block_indices[[tbind]]
    x$ntables *  project(x$fit, xnewdat, 1:ncomp) 
      
  })
  
  names(res) <- paste0("table_", table_index)
  res
} 

#' @export
singular_values.sca <- function(x) x$fit$d


#' @export
scores.sca <- function(x) {
  x$fit$scores
}

#' @export
loadings.sca <- function(x) {
  x$fit$v
}

#' @export 
block_index_list.sca <- function(x) x$block_indices


#' @export
ncomp.sca <- function(x) length(x$d)

#' @export
reconstruct.sca <- function(x, ncomp=x$ncomp) {
  ret <- block_matrix(lapply(x$D, function(d) {
    t(tcrossprod(d, x$B))
  }))
}

#' @export
contributions.sca <- function(x, type=c("column", "row", "table")) {
  contr <- contributions(x$fit)
  type <- match.arg(type)
  
  if (type == "table") {
    out <- do.call(cbind, lapply(1:x$ncomp, function(i) {
      sapply(1:x$ntables, function(j) {
        ind <- x$block_indices[[j]]
        sum(contr[i,ind])
      })
    }))
    
  } else if (type == "row") {
    contributions(x$fit, type="row")
  } else {
    contributions(x$fit, type="column")
  }
  
}

predict.sca <- function(x, newdata, ncomp=x$ncomp, table_index=1:x$ntables, pre_process=TRUE) {
  assert_that(is.matrix(newdata))
  assert_that(length(table_index) == 1 || length(table_index) == x$ntables)
  
  fscores <- if (length(table_index) == x$ntables) {
    assert_that(ncol(newdata) == ncol(x$X))
    Reduce("+", lapply(table_index, function(i) {
      ind <- x$block_indices[[i]]
      
      Xp <- if (pre_process) {
        x$reprocess(newdata[, ind], i)
      } else {
        newdata[, ind]
      }
      
      project(x, Xp, ncomp=ncomp, table_index=i) 
    }))
  } else if (length(table_index) == 1) {
    
    ind <- x$block_indices[[table_index]]
    assert_that(length(ind) == ncol(newdata))
    
    fscores <- project(x, newdata, ncomp=ncomp, table_index=table_index) 
    
  }
  
}

