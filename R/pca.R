
## TODO projector needs to be better defined. May need a class called "dimred": "projector" (X -> D), "dimred" (orthogonal or non-orthogonal), "pca" (orthogonal)


projector <- function(preproc, ncomp, v, classes, ...) {
  out <- list(
    preproc=preproc,
    ncomp=ncomp,
    v=v,
    ...)
  
  class(out) <- c(classes, "projector")
  out
}

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
    


kmeans_projector <- function(X, ncomp=max(as.integer(nrow(X)/2),2),center=TRUE, scale=FALSE,  ...) {

  preproc <- pre_processor(X, center, scale)
  Xp <- pre_process(preproc, X)
  
  kres <- kmeans(t(Xp), centers=ncomp,...)
  mat <- do.call(rbind, purrr::imap(split(1:ncol(Xp), kres$cluster), function(.x, .y) {
  cbind(as.numeric(.x), as.numeric(.y), 1/length(.x))

  }))
  
  v <- sparseMatrix(i=mat[,1], j=mat[,2], x=mat[,3])
  
  cen <- t(kres$centers)
  u <- apply(cen, 2, function(vals) normalize(as.matrix(vals)))
  
  d <- apply(cen,2,function(x) sqrt(sum(x^2)))
  dord <- order(d, decreasing=TRUE)
  d <- d[dord]
  
  clus <- sapply(kres$cluster, function(i) { dord[i] })
  ret <- bi_projector(
                preproc=preproc,
                ncomp=ncomp,
                v=v[,dord],
                u=u[,dord], 
                d=d,
                scores=cen[,dord],
                clusters=clus,
                withinss=kres$withinss[dord],
                totss=kres$totss,
                kres$betweenss,
                classes=c("kmeans_projector"))
  ret
  
}

#' @export
reconstruct.kmeans_projector <- function(x, newdata=NULL, comp=1:x$ncomp, subind=NULL) {
  if (!is.null(newdata)) {
    assert_that(ncol(newdata) == length(comp) && nrow(newdata) == nrow(scores(x)))
  } else {
    newdata <- x$scores[,comp, drop=FALSE]
  }
  
  ##recon2 = predict(rst) %*% ginv(rst$rotation) + matrix(1,5,1) %*% rst$center
  lds <- x$v[,comp,drop=FALSE]
  
  if (is.null(subind)) {
    reverse_pre_process(x$preproc, newdata %*% t(lds))
  } else {
    reverse_pre_process(x$preproc, newdata %*% t(lds)[,subind], subind=subind)
  }
}



#' nneg_pca
#' 
#' non-negative pca
#' 
#'   
#' @inheritParams pca
#' @importFrom nsprcomp nsprcomp
#' @export
nneg_pca <- function(X, ncomp=min(dim(X)), center=TRUE, scale=FALSE,  ...) {
  preproc <- pre_processor(X, center, scale)
  Xp <- pre_process(preproc, X)
  
  ret <- nsprcomp::nsprcomp(Xp, center=FALSE, scale=FALSE, nneg=TRUE, ...)
  keep <- ret$sdev > 1e-06
  
  v=ret$rotation[,keep,drop=FALSE]
  u=ret$x[,keep,drop=FALSE]
  d=ret$sdev[keep] * 1/sqrt(max(1,nrow(X)-1))
  
  ret <- bi_projector(
              preproc=preproc,
              ncomp = length(d),
              v=v, 
              u=u,
              d=d,
              scores=ret$x,
              classes=c("nneg_pca", "pca"))
  
  ret
  
}

#' @export
reconstruct.nneg_pca <- function(x, newdata=NULL, comp=1:x$ncomp, subind=NULL) {
  if (!is.null(newdata)) {
    assert_that(ncol(newdata) == length(comp) && nrow(newdata) == nrow(scores(x)))
  } else {
    newdata <- scores(x)[,comp, drop=FALSE]
  }
  
  ##recon2 = predict(rst) %*% ginv(rst$rotation) + matrix(1,5,1) %*% rst$center
  
  piv <- corpcor::pseudoinverse(loadings(x)[,comp,drop=FALSE])
  if (is.null(subind)) {
    reverse_pre_process(x$preproc, newdata %*% piv)
  } else {
    reverse_pre_process(x$preproc, newdata %*% piv[,subind], subind=subind)
  }
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
shrink_pca <- function(X, center=TRUE, scale=FALSE,  method = c("GSURE", "QUT", "SURE"), ...) {
  preproc <- pre_processor(X, center, scale)
  Xp <- pre_process(preproc, X)
  
  res <- denoiseR::adashrink(Xp, method=method, center=FALSE, ...)
  
  keep <- res$singval > 1e-06
  
  if (sum(keep) == 0) {
    keep <- 1
  } 
  
  
  v=res$low.rank$v[,keep,drop=FALSE]
  u=res$low.rank$u[,keep,drop=FALSE]
  d=res$low.rank$d[keep]

  ret <- bi_projector(
              preproc=preproc,
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


#' @param X
#' @param ncomp
#' @param center
#' @param scale
#' @export
pca <- function(X, ncomp=min(dim(X)), center=TRUE, scale=FALSE, ...) {
  assert_that(is.matrix(X) || inherits(X, "Matrix"))
  
  preproc <- pre_processor(X, center=center,scale=scale)
  Xp <- pre_process(preproc,X)
  
  svdres <- svd_wrapper(Xp, ncomp, ...)
  
  scores <- t(t(as.matrix(svdres$u)) * svdres$d)
  
  if (!is.null(row.names(scores))) {
    row.names(scores) <- row.names(Xp)[seq_along(svdres$d)]
  }
  
  ret <- bi_projector(
              preproc,
              length(svdres$d),
              v=svdres$v, 
              u=svdres$u, 
              d=svdres$d, 
              scores=scores,
              classes=c("pca"))

  ret
}


#' @export

reprocess.bi_projector <- function(x, newdata, subind=NULL) {
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
  
  function(newdata, subind=NULL) {
    .subind <- subind
    project(x, newdata, comp=.comp, pre_process=.pre_process, subind=.subind)
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
project.projector <- function(x, newdata, comp=1:ncomp(x), pre_process=TRUE, subind=NULL) {
  
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

residuals.bi_projector <- function(x, ncomp=1, xorig) {
  recon <- reconstruct(x,comp=1:ncomp)
  #recon <- scores(x)[,1:ncomp,drop=FALSE] %*% t(loadings(x)[,1:ncomp,drop=FALSE])
  orig <- pre_process(x$preproc, xorig)
  orig - recon
}


#' @export
reconstruct.bi_projector <- function(x, newdata=NULL, comp=1:x$ncomp, subind=NULL) {
  if (!is.null(newdata)) {
    assert_that(ncol(newdata) == length(comp) && nrow(newdata) == nrow(scores(x)))
  } else {
    newdata <- scores(x)[,comp, drop=FALSE]
  }
  
  if (is.null(subind)) {
    reverse_pre_process(x$preproc, newdata %*% t(loadings(x)[,comp,drop=FALSE]))
  } else {
    reverse_pre_process(x$preproc, newdata %*% t(loadings(x)[,comp,drop=FALSE])[,subind], subind=subind)
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
    
    class(ret) <- c("pca", "projector", "list")
  }
  
  ret
}


#' @export
contributions.projector <- function(x, type=c("column", "row")) {
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

print.bi_projector <- function(object) {
  showk <- 1:min(object$ncomp, 5)
  cat("pca", "\n")
  cat("  number of components: ", object$ncomp, "\n")
  cat("  number of variables: ", nrow(object$v), "\n")
  cat("  number of observations: ", nrow(object$u), "\n")
  cat("  center: ", object$preproc$center, "\n")
  cat("  scale: ", object$preproc$scale, "\n")
  cat("  % variance explained (top 5): ", ((object$d[showk]^2)/sum(object$d^2)) * 100, "\n")
}





