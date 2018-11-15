
#' @export
reconstruct.kmeans_projector <- function(x, newdata=NULL, comp=1:x$ncomp, colind=NULL) {
  if (!is.null(newdata)) {
    assert_that(ncol(newdata) == length(comp) && nrow(newdata) == nrow(scores(x)))
  } else {
    newdata <- x$scores[,comp, drop=FALSE]
  }
  
  ##recon2 = predict(rst) %*% ginv(rst$rotation) + matrix(1,5,1) %*% rst$center
  lds <- x$v[,comp,drop=FALSE]
  
  if (is.null(colind)) {
    reverse_pre_process(x$preproc, newdata %*% t(lds))
  } else {
    reverse_pre_process(x$preproc, newdata %*% t(lds)[,colind], colind=colind)
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
reconstruct.nneg_pca <- function(x, newdata=NULL, comp=1:x$ncomp, colind=NULL) {
  if (!is.null(newdata)) {
    assert_that(ncol(newdata) == length(comp) && nrow(newdata) == nrow(scores(x)))
  } else {
    newdata <- scores(x)[,comp, drop=FALSE]
  }
  
  ##recon2 = predict(rst) %*% ginv(rst$rotation) + matrix(1,5,1) %*% rst$center
  
  piv <- corpcor::pseudoinverse(loadings(x)[,comp,drop=FALSE])
  if (is.null(colind)) {
    reverse_pre_process(x$preproc, newdata %*% piv)
  } else {
    reverse_pre_process(x$preproc, newdata %*% piv[,colind], colind=colind)
  }
}




#' @keywords internal
nipals <- function(X, center=TRUE, scale=FALSE, ncomp=min(dim(X)), thresh=1e-5) {
  
  iterate <- function(E) {
    t <- E[,sample(1:ncol(X),1), drop=FALSE]
    crit <- .Machine$integer.max
    
    while (crit > thresh) {
      told <- t
      p <- crossprod(E,t)/sum(t^2)
      p <- p * (sum(p^2))^-.5
      t <- (E %*% p)/sum(p^2)
      crit <- abs(sum(t^2) - sum(told^2))
    }
    
    list(p=p, t=t)
  }
  
  last <- NULL
  p <- matrix(0, nrow(X), ncomp)
  t <- matrix(0, ncol(X), ncomp)
  
  for (i in 1:ncomp) {
    if (i == 1) {
      E <- X
    } else {
      E <- E - (last$t %*% t(last$p))
    }
    
    last <- iterate(E)
    t[,i] <- last$t
    p[,i] <- last$p
    
  }
  
  list(t=t, p=p)
}