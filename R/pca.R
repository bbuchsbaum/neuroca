

## TODO projector needs to be better defined. May need a class called "dimred": "projector" (X -> D), "dimred" (orthogonal or non-orthogonal), "pca" (orthogonal)

#' construct a `projector` instance
#' @export
#' @param preproc
#' @param ncomp
#' @param v
#' @param classes
projector <- function(preproc, ncomp, v, classes, ...) {
  out <- list(
    preproc=preproc,
    ncomp=ncomp,
    v=v,
    ...)
  
  class(out) <- c(classes, "projector")
  out
}


#' construct a `projector` instance
#' @export
#' @inheritParams projector
#' @param u
#' @param d
#' @param scores
bi_projector <- function(preproc, ncomp, v, u, d, scores, classes, ...) {
  out <- list(
    preproc=preproc,
    ncomp=ncomp,
    v=v,
    u=u,
    d=d,
    scores=scores,
    ...)
  
  class(out) <- c(classes, "bi_projector", "projector")
  out
}
    



#' shrink_pca
#' 
#' adaptive shrinkage pca from the \code{denoiseR} package
#' 
#'   
#' @param X
#' @param center
#' @param scale
#' @importFrom denoiseR adashrink
#' @export
shrink_pca <- function(X, preproc=center(), method = c("GSURE", "QUT", "SURE"), ...) {
  assert_that(is.matrix(X) || inherits(X, "Matrix"))
  
  procres <- prep(preproc, X)
  Xp <- procres$Xp
  
  res <- denoiseR::adashrink(Xp, method=method, center=FALSE, ...)
  
  keep <- res$singval > 1e-06
  
  if (sum(keep) == 0) {
    keep <- 1
  } 
  
  v=res$low.rank$v[,keep,drop=FALSE]
  u=res$low.rank$u[,keep,drop=FALSE]
  d=res$low.rank$d[keep]

  ret <- bi_projector(
              preproc=procres,
              ncomp=length(d),
              v=v, 
              u=u,
              d=d,
              scores=t(t(as.matrix(u)) * d),
              classes=c("shrink_pca", "pca"))

}
                       

#' pseudo_svd
#' 
#' @param u the row weights
#' @param v the column weights
#' @param d the singular values
#' @param rnames row names of observations
#' @export
pseudo_svd <- function(u, v, d, rnames=NULL) {
  
  scores <- t(t(as.matrix(u)) * d)
  
  if (!is.null(rnames)) {
    row.names(scores) <- rnames
  }
  
  assert_that(length(d) == ncol(u))
  assert_that(ncol(v) == ncol(u))
  

  ret <- bi_projector(
              preproc=pre_processor(matrix(), center=FALSE, scale=FALSE),
              ncomp=length(d),
              v=v, 
              u=u, 
              d=d, 
              scores=scores,
              classes="pseudo_svd")
  

  ret
}

#' principal components analysis
#' 
#' @param X
#' @param ncomp
#' @param preproc
#' @export
#' 
#' @examples 
#' 
#' X <- matrix(rnorm(10*20), 10, 20)
#' res <- pca(X, ncomp=3, preproc=standardize())
pca <- function(X, ncomp=min(dim(X)), preproc=center(), ...) {
  assert_that(is.matrix(X) || inherits(X, "Matrix"))
  
  procres <- prep(preproc, X)
  
  Xp <- procres$Xp
  
  svdres <- svd_wrapper(Xp, ncomp, ...)
  
  scores <- t(t(as.matrix(svdres$u)) * svdres$d)
  
  if (!is.null(row.names(scores))) {
    row.names(scores) <- row.names(Xp)[seq_along(svdres$d)]
  }
  
  ret <- bi_projector(
              preproc=procres,
              length(svdres$d),
              v=svdres$v, 
              u=svdres$u, 
              d=svdres$d, 
              scores=scores,
              classes=c("pca"))
  ret
}

#' @export
refit.pca <- function(x, X) {
  pca(X, ncomp=x$ncomp, preproc=x$preproc$preproc)
}

#' @export
permutation.pca <- function(x, X, nperm=100) {
  
  evals <- x$d^2
  Fa <- sapply(1:length(evals), function(i) evals[i]/sum(evals[i:length(evals)]))
  
  F1_perm <- sapply(1:nperm, function(i) {
    Xperm <- apply(X, 2, function(x) sample(x))
    fit <- refit(x, Xperm)
    evals <- fit$d^2
    F1_perm <- evals[1]/sum(evals[1:length(evals)])
  })
  
  if (x$ncomp > 1) {
    for (i in 2:x$ncomp) {
      Ea <- residuals(x, X, i)
      
    }
  }
}
 
#' @export   
residuals.pca <- function(x, X, ncomp) {
  recon <- reconstruct(x,ncomp)
  X - recon
}


#' @export
rotate.pca <- function(x, rot) {
  u_rot <- x$u %*% rot
  v_rot <- x$v %*% rot
  sc_rot <- scores(x) %*% rot
  
  ret <- bi_projector(
    preproc=x$preproc,
    ncomp=length(x$d),
    v=v_rot, 
    u=u_rot, 
    d=apply(sc_rot, 2, function(x) sqrt(sum(x^2))),
    scores=sc_rot, 
    classes=c("rotated_pca", "pca"))
  
  ret
  
}


#' @export
reprocess.bi_projector <- function(x, newdata, colind=NULL) {
  if (is.null(colind)) {
    assert_that(ncol(newdata) == nrow(loadings(x)))
    x$preproc$transform(newdata)
  } else {
    assert_that(length(colind) == ncol(newdata), 
                msg=paste("length of colind not equal to number of columns of newdata", length(colind), "!=", ncol(newdata)))
    x$preproc$transform(newdata, colind)
  }
  
}

#' @export
singular_values.bi_projector <- function(x) {
  x$d
}

#' @export
loadings.bi_projector <- function(x) {
  x$v
}


#' @export
projection_fun.bi_projector <- function(x, comp=1:ncomp(x), pre_process=TRUE) {
  .comp <- comp
  .pre_process <- pre_process
  
  function(newdata, colind=NULL) {
    .colind <- colind
    project(x, newdata, comp=.comp, pre_process=.pre_process, colind=.colind)
  }
}

#' @export
project_cols.bi_projector <- function(x, newdata, comp=1:ncomp(x)) {

  ## if no newdata, then simply return the loadings
  if (missing(newdata)) {
    return(loadings(x)[,comp, drop=FALSE])
  } else {
    if (is.vector(newdata)) {
      newdata <- as.matrix(newdata, ncol=1)
    }
    assert_that(nrow(newdata) == nrow(x$u))
    t(newdata) %*% (x$u[,comp,drop=FALSE]) %*% diag(1/x$d[comp], nrow=length(comp), ncol=length(comp))
  }
}

#' @export
project.projector <- function(x, newdata, comp=1:ncomp(x), colind=NULL) {
  ## if no newdata, then simply return the factor scores
  if (missing(newdata)) {
    return(scores(x)[,comp, drop=FALSE])
  }
  
  if (is.vector(newdata)) {
    newdata <- matrix(newdata,nrow=1)
  }
  
  if (is.null(colind)) {
    reprocess(x, newdata) %*% loadings(x)[,comp, drop=FALSE]
  } else {
    ## colind must be in the input space
    assertthat::assert_that(max(colind) <= ncol(newdata))
    reprocess(x, newdata, colind=colind) %*% (loadings(x)[colind, comp, drop=FALSE])
  }
}

#' @export
residuals.bi_projector <- function(x, ncomp=1, xorig) {
  recon <- reconstruct(x,comp=1:ncomp)
  xorig - recon
}


#' @export
reconstruct.bi_projector <- function(x, newdata=NULL, comp=1:x$ncomp, colind=NULL, rowind=NULL, reverse_pre_process=TRUE) {
  if (!is.null(newdata)) {
    assert_that(ncol(newdata) == length(comp) && nrow(newdata) == nrow(scores(x)))
  } else {
    newdata <- scores(x)[,comp, drop=FALSE]
  }
  
  if (is.null(rowind)) {
    rowind <- 1:nrow(newdata)
  } else {
    assert_that(all(rowind > 0) && max(rowind) < nrow(newdata))
  }
  
  if (is.null(colind)) {
    if (reverse_pre_process) {
      x$preproc$reverse_transform(newdata[rowind,,drop=FALSE] %*% t(loadings(x)[,comp,drop=FALSE]))
    } else {
      newdata[rowind,,drop=FALSE] %*% t(loadings(x)[,comp,drop=FALSE])
    }
  } else {
    if (reverse_pre_process) {
      x$preproc$reverse_transform(newdata[rowind,,drop=FALSE] %*% t(loadings(x)[,comp,drop=FALSE])[,colind], 
                        colind=colind)
    } else {
      newdata[rowind,,drop=FALSE] %*% t(loadings(x)[,comp,drop=FALSE])[,colind]
    }
  }
}


#' @export
ncol.projector <- function(x) {
  nrow(x$v)
}

#' @export
nrow.bi_projector <- function(x) {
  nrow(x$u)
}

#' @export
ncomp.projector <- function(x) {
  x$ncomp
}

#' @export
scores.bi_projector <- function(x) {
  x$scores
}


#' @export
dim.projector <- function(x) {
  c(nrow(x$u), nrow(x$v))
}



#' @export
truncate.bi_projector <- function(obj, ncomp) {
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
    
    class(ret) <- c("bi_projector", "list")
  }
  
  ret
}


#' @export
contributions.projector <- function(x) {
  loadings(x)^2
}


#' @export
contributions.bi_projector <- function(x, type=c("column", "row")) {
  type <- match.arg(type)
  if (type == "column") {
    loadings(x)^2
  } else {
    x$u^2
  }
}




#' @export
print.bi_projector <- function(object) {
  showk <- 1:min(object$ncomp, 5)
  cat(class(object)[1], "\n")
  cat("  number of components: ", object$ncomp, "\n")
  cat("  number of variables: ", nrow(object$v), "\n")
  cat("  number of observations: ", nrow(object$u), "\n")
  cat("  % variance explained (top 5): ", ((object$d[showk]^2)/sum(object$d^2)) * 100, "\n")
}





