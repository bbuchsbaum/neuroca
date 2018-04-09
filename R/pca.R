

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
  xdim <- dim(X)
  proj_fun <- projection_fun(X)
  
  Xp <- pre_processor(project(X), center,scale)
  res <- adashrink(Xp, method=method, center=FALSE, ...)
  
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
  
  ret <- list(proj_fun=proj_fun,
              xdim=xdim,
              v=v, 
              u=u,
              d=d,
              scores=t(t(as.matrix(u)) * d),
              ncomp=length(d), 
              pre_process=attr(Xp, "pre_process"), 
              reverse_pre_process=attr(Xp, "reverse_pre_process"))
  
  
  class(ret) <- c("shrink_pca", "pca", "projector", "list")
  ret
}
                       


#' pseudo_pca
#' 
#' @export
#' @param u the row weights
#' @param v the column weights
#' @param d the singular values
#' @param pre_process the function to pre_process the data
#' @param reverse_pre_process the function to reverse pre_process the data
#' @param row names of observations
pseudo_pca <- function(u, v, d, 
                       pre_process=function(x, subind) x, 
                       reverse_pre_process=function(x, subind) x, rnames=NULL) {
  
  scores <- t(t(as.matrix(u)) * d)
  
  if (!is.null(rnames)) {
    row.names(scores) <- rnames
  }
  
  assert_that(length(d) == nrow(u))
  
  ret <- list(xdim=c(length(d), nrow(v)),
              proj_fun=identity,
              v=v, 
              u=u, 
              d=d, 
              scores=scores, ncomp=length(d), 
              pre_process=pre_process, 
              reverse_pre_process=reverse_pre_process)
  
  
  class(ret) <- c("pseudo_pca", "pca", "projector", "list")
  ret
}


#' @param X
#' @param ncomp
#' @param center
#' @param scale
#' @export
pca <- function(X, ncomp=min(dim(X)), center=TRUE, scale=FALSE, ...) {
  xdim <- dim(X)
  proj_fun <- projection_fun(X)
  
  Xp <- pre_processor(project(X), center=center,scale=scale)
  
  svdres <- svd_wrapper(Xp, ncomp, ...)
  
  scores <- t(t(as.matrix(svdres$u)) * svdres$d)
  
  if (!is.null(row.names(scores))) {
    row.names(scores) <- row.names(Xp)[seq_along(svdres$d)]
  }
  
  ret <- list(xdim=xdim,
              proj_fun=proj_fun,
              v=svdres$v, u=svdres$u, d=svdres$d, 
              scores=scores, ncomp=length(svdres$d), 
              pre_process=attr(Xp, "pre_process"), 
              reverse_pre_process=attr(Xp, "reverse"))

  
  class(ret) <- c("pca", "projector", "list")
  ret
}


#' @export
reprocess.pca <- function(x, newdata, subind=1:ncol(newdata)) {
  x$pre_process(x$proj_fun(newdata[,subind,drop=FALSE]))
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
projection_fun.pca <- function(x, comp=1:ncomp(x), pre_process=TRUE, subind=NULL, ...) {
  function(newdata, .comp=comp, .pre_process=pre_process, .subind=subind) {
    project(x, newdata, comp=.comp, pre_process=.pre_process, subind=.subind)
  }
}

#' @export
project.pca <- function(x, newdata, comp=1:ncomp(x), pre_process=TRUE, subind=NULL) {
  if (missing(newdata)) {
    return(scores(x)[,comp])
  }
  
  if (is.vector(newdata)) {
    newdata <- matrix(newdata,nrow=1)
  }
  
  if (is.null(subind)) {
    ## subind must be in the input space
    if (pre_process) {
      reprocess(x, newdata) %*% loadings(x)[,comp]
    } else {
      x$proj_fun(newdata) %*% loadings(x)[,comp]
    }
  } else {
    assertthat::assert_that(length(subind) < ncol(newdata))
    
    Xsup <- if (pre_process) {
      reprocess(x, newdata, subind=subind)
    } else {
      x$proj_fun(newdata, subind=subind)
    }
    
    Xsup %*% (loadings(x)[, comp])
  }
}


#' @export
residuals.pca <- function(x, ncomp=1, xorig) {
  recon <- scores(x)[,1:ncomp,drop=FALSE] %*% t(loadings(x)[,1:ncomp,drop=FALSE])
  newdat <- x$pre_process(xorig)
  newdat - recon
}


#' @export
reconstruct.pca <- function(x, comp=1:x$ncomp) {
  x$reverse_pre_process(scores(x)[,comp,drop=FALSE] %*% t(loadings(x)[,comp,drop=FALSE]))
}


#' @export
ncol.pca <- function(x) {
  x$xdim[2]
}

#' @export
nrow.pca <- function(x) {
  x$xdim[1]
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
    ret <- list(v=obj$v[,1:ncomp], u=obj$u[,1:ncomp], d=obj$d[1:ncomp], 
                scores=obj$scores[,1:ncomp], ncomp=ncomp, 
                pre_process=obj$pre_process,
                reverse_pre_process=obj$pre_process)
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





# 
# hierarchical_pca <- function(X, comps, hierarchy, svd.method="base") {
# }