
#' @importFrom RGCCA rgcca
canautocor <- function(X, preproc=pass(), ncomp=2, tau=c(.2,.2), npcs=NULL) {
  
  Xorig <- X
  
  if (!is.null(npcs)) {
    pcres <- pca(Xorig, ncomp=npcs, preproc=pass())
    X <- scores(pcres)
  } else {
    pcres <- NULL
    X <- Xorig
  }
  
  X1 <- X[1:nrow(X)-1,]
  X2 <- X[2:nrow(X),]
  
  fit <- RGCCA::rgcca(A= list(X1, X2), ncomp=rep(ncomp,2),
        C = matrix(c(0, 1, 1, 0), 2, 2),
        tau = tau)
  
  scores <- X %*% fit$a[[1]]
  
  lds <- if (!is.null(pcres)) {
    loadings(pcres) %*% fit$a[[1]]
  } else {
    fit$a[[1]]
  }
  
  ret <- list(X=Xorig, scores=scores, loadings=lds, fit=fit, ncomp=ncomp, tau=tau, pcres=pcres)
  class(ret) <- "canautocor"
  ret
}

scores.canautocor <- function(x) {
  x$scores
}

loadings.canautocor <- function(x) {
  x$loadings
}
  