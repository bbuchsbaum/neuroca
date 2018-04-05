

#' shrink_pca
#' 
#' adaptive shrinakge pca from the \code{denoiseR} package
#' 
#'   
#' @param X
#' @param center
#' @param scale
#' @export
shrink_pca <- function(X, center=TRUE, scale=FALSE,  method = c("GSURE", "QUT", "SURE"), ...) {
  X <- pre_processor(X, center,scale)
  res <- adashrink(X, method=method, center=FALSE, ...)
  
  keep <- res$singval > 1e-06
  if (sum(keep) == 0) {
    warning("all singular values are zero, computing rank-1 svd")
    res <- RSpectra::svds(X, k=1)
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
              pre_process=attr(X, "pre_process"), 
              reverse_pre_process=attr(X, "reverse_pre_process"))
  
  
  class(ret) <- c("shrink_pca", "pseudo_pca", "projector", "list")
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
  
  ret <- list(v=v, u=u, d=d, 
              scores=scores, ncomp=length(d), 
              pre_process=pre_process, 
              reverse_pre_process=reverse_pre_process)
  
  
  class(ret) <- c("pseudo_pca", "projector", "list")
  ret
}

#' @export
truncate.pseudo_pca <- function(obj, ncomp) {
  if (ncomp >= obj$ncomp) {
    warning("number of components to keep is greater than or equal to rank of pca fit, returning original model")
    ret <- obj
  } else {
    ret <- list(v=obj$v[,1:ncomp], u=obj$u[,1:ncomp], d=obj$d[1:ncomp], 
                scores=obj$scores[,1:ncomp], ncomp=ncomp, 
                pre_process=obj$pre_process)
    class(ret) <- c("pseudo_pca", "projector", "list")
  }
  
  ret
}
  
  
  
#' @param X
#' @param ncomp
#' @param center
#' @param scale
#' @export
pca <- function(X, ncomp=min(dim(X)), center=TRUE, scale=FALSE, ...) {
  
  X <- pre_processor(X, center=center,scale=scale)
  
  svdres <- svd_wrapper(X, ncomp, ...)
  
  scores <- t(t(as.matrix(svdres$u)) * svdres$d)
  row.names(scores) <- row.names(X)
  
  ret <- list(v=svdres$v, u=svdres$u, d=svdres$d, 
              scores=scores, ncomp=ncomp, 
              pre_process=attr(X, "pre_process"), reverse_pre_process=attr(X, "reverse"))

  
  class(ret) <- c("pca", "projector", "list")
  ret
}


#' @export
reprocess.projector <- function(x, newdata, subind=NULL) {
  x$pre_process(newdata, subind)
}

#' @export
singular_values.pseduo_pca <- function(x) {
  x$d
}


#' @export
singular_values.pca <- function(x) {
  x$d
}

#' @export
loadings.projector <- function(x) {
  x$v
}

#' @export
ncomp.projector <- function(x) {
  length(x$d)
}

#' @export
scores.projector <- function(x) {
  x$scores
}

#' @export
project.projector <- function(x, newdata, comp=1:ncomp(x), pre_process=TRUE, subind=NULL) {
  if (is.vector(newdata)) {
    newdata <- matrix(newdata,nrow=1)
  }
  if (is.null(subind)) {
    if (pre_process) {
      reprocess(x, newdata) %*% x$v[,comp]
    } else {
      newdata %*% x$v[,comp]
    }
  } else {
    assertthat::assert_that(length(subind) == ncol(newdata))
    Xsup <- if (pre_process) {
      #browser()
      reprocess(x, newdata, subind)
    } else {
      newdata
    }
    Xsup %*% x$v[subind, comp]
  }
}

#' @export
predict.projector <- function(x, newdata, ncomp=ncomp(x)) {
  project(x, newdata, comp=1:ncomp)
}

#' @export
residuals.projector <- function(x, ncomp=1, xorig) {
  recon <- x$scores[,1:ncomp,drop=FALSE] %*% t(x$v[,1:ncomp,drop=FALSE])
  newdat <- x$pre_process(xorig)
  newdat - recon
}


#' @export
reconstruct.projector <- function(x, comp=1:x$ncomp) {
  x$reverse_pre_process(x$scores[,comp,drop=FALSE] %*% t(x$v[,comp,drop=FALSE]))
}



#' @export
truncate.pca <- function(obj, ncomp) {
  if (ncomp >= obj$ncomp) {
    warning("number of components to keep is greater than or equal to rank of pca fit, returning original model")
    ret <- obj
  } else {
    ret <- list(v=obj$v[,1:ncomp], u=obj$u[,1:ncomp], d=obj$d[1:ncomp], 
                scores=obj$scores[,1:ncomp], ncomp=ncomp, 
                pre_process=obj$pre_process)
    class(ret) <- c("pca", "projector", "list")
  }
  
  ret
}


#' @export
contributions.projector <- function(x, type=c("column", "row")) {
  type <- match.arg(type)
  if (type == "column") {
    x$v^2
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