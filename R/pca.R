

#' shrink_pca
#' 
#' adaptive shrinakge pca from the \code{denoiseR} package
#' 
#'   
#' @param X
#' @param center
#' @param scale
#' @importFrom denoiseR adashrink
#' @export
shrink_pca <- function(X, center=TRUE, scale=FALSE,  method = c("GSURE", "QUT", "SURE"), ...) {
  preproc <- pre_processor(X, center, scale)
  Xp <- pre_process(preproc, X)
  
  res <- denoiseR::adashrink(Xp, method=method, center=FALSE, ...)
  
  keep <- res$singval > 1e-06
  
  if (sum(keep) == 0) {
    warning("all singular values are zero, computing rank-1 svd")
    res <- RSpectra::svds(Xp, k=1)
    v <- res$v
    u <- res$u
    d <- res$d
  } else {
    v=res$low.rank$v[,keep,drop=FALSE]
    u=res$low.rank$u[,keep,drop=FALSE]
    d=res$low.rank$d[keep]
  }
  
  ret <- list(v=v, 
              u=u,
              d=d,
              scores=t(t(as.matrix(u)) * d),
              ncomp=length(d), 
              preproc=preproc)
  
  
  class(ret) <- c("shrink_pca", "pca", "projector", "list")
  ret
}
                       


#' pseudo_pca
#' 
#' @param u the row weights
#' @param v the column weights
#' @param d the singular values
#' @param rnames row names of observations
pseudo_pca <- function(u, v, d, rnames=NULL) {
  
  scores <- t(t(as.matrix(u)) * d)
  
  if (!is.null(rnames)) {
    row.names(scores) <- rnames
  }
  
  assert_that(length(d) == nrow(u))
  
  ret <- list(v=v, 
              u=u, 
              d=d, 
              scores=scores, 
              ncomp=length(d), 
              preproc=pre_processor(matrix(), center=FALSE, scale=FALSE))
  
  
  class(ret) <- c("pseudo_pca", "pca", "projector", "list")
  ret
}


#' @param X
#' @param ncomp
#' @param center
#' @param scale
#' @export
pca <- function(X, ncomp=min(dim(X)), center=TRUE, scale=FALSE, ...) {
  
  preproc <- pre_processor(X, center=center,scale=scale)
  Xp <- pre_process(preproc)
  
  svdres <- svd_wrapper(Xp, ncomp, ...)
  
  scores <- t(t(as.matrix(svdres$u)) * svdres$d)
  
  if (!is.null(row.names(scores))) {
    row.names(scores) <- row.names(Xp)[seq_along(svdres$d)]
  }
  
  ret <- list(v=svdres$v, 
              u=svdres$u, 
              d=svdres$d, 
              scores=scores, 
              ncomp=length(svdres$d), 
              preproc=preproc)

  
  class(ret) <- c("pca", "projector", "list")
  ret
}


#' @export
reprocess.pca <- function(x, newdata, subind=NULL) {
  if (is.null(subind)) {
    assert_that(ncol(newdata) == ncol(x))
    pre_process(x$preproc, newdata)
  } else {
    assert_that(length(subind) == ncol(newdata), 
                msg=paste("length of subind not equal to number of columns of newdata", length(subind), "!=", ncol(newdata)))
    pre_process(x$preproc, newdata, subind)
  }
  
}

#' @export
singular_values.pca <- function(x) {
  x$d
}

#' @export
loadings.pca <- function(x) {
  x$v
}


#' @export
projection_fun.pca <- function(x, comp=1:ncomp(x), pre_process=TRUE) {
  .comp <- comp
  .pre_process <- pre_process
  
  function(newdata, subind=NULL) {
    .subind <- subind
    project(x, newdata, comp=.comp, pre_process=.pre_process, subind=.subind)
  }
}


#' @export
project.pca <- function(x, newdata, comp=1:ncomp(x), pre_process=TRUE, subind=NULL) {
  
  ## if no newdata, then simply return the factor scores
  if (missing(newdata)) {
    return(scores(x)[,comp, drop=FALSE])
  }
  
  if (is.vector(newdata)) {
    newdata <- matrix(newdata,nrow=1)
  }
  
  if (is.null(subind)) {
    reprocess(x, newdata) %*% loadings(x)[,comp, drop=FALSE]
  } else {
    ## subind must be in the input space
    assertthat::assert_that(max(subind) <= ncol(newdata))
    reprocess(x, newdata, subind=subind) %*% (loadings(x)[subind, comp, drop=FALSE])
  }
}

#' @export
residuals.pca <- function(x, ncomp=1, xorig) {
  recon <- reconstruct(x,comp=1:ncomp)
  #recon <- scores(x)[,1:ncomp,drop=FALSE] %*% t(loadings(x)[,1:ncomp,drop=FALSE])
  orig <- pre_process(x$preproc, xorig)
  orig - recon
}


#' @export
reconstruct.pca <- function(x, newdata=NULL, comp=1:x$ncomp) {
  if (!is.null(newdata)) {
    assert_that(ncol(newdata) == length(comp) && nrow(newdata) == nrow(scores(x)))
  } else {
    newdata <- scores(x)[,comp, drop=FALSE]
  }
  reverse_pre_process(x$preproc, newdata %*% t(loadings(x)[,comp,drop=FALSE]))
}


#' @export
ncol.pca <- function(x) {
  nrow(x$v)
}

#' @export
nrow.pca <- function(x) {
  nrow(x$u)
}

#' @export
ncomp.pca <- function(x) {
  x$ncomp
}

#' @export
scores.pca <- function(x) {
  x$scores
}


#' @export
truncate.pca <- function(obj, ncomp) {
  if (ncomp >= obj$ncomp) {
    warning("number of components to keep is greater than or equal to rank of pca fit, returning original model")
    ret <- obj
  } else {
    ret <- list(v=obj$v[,1:ncomp,drop=FALSE], 
                u=obj$u[,1:ncomp,drop=FALSE], 
                d=obj$d[1:ncomp], 
                scores=obj$scores[,1:ncomp,drop=FALSE], 
                ncomp=ncomp, 
                preproc=obj$preproc)
    
    class(ret) <- c("pca", "projector", "list")
  }
  
  ret
}


#' @export
contributions.pca <- function(x, type=c("column", "row")) {
  type <- match.arg(type)
  if (type == "column") {
    loadings(x)^2
  } else {
    x$u^2
  }
}

#' @export
reduce_rank.matrix <- function(x, k=min(dim(x)), center=TRUE, scale=FALSE, 
                               reducer=pca, ...) {
  res <- reducer(x, k, center=center, scale=scale, ...)
}

