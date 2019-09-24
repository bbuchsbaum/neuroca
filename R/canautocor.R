
#' @importFrom RGCCA rgcca
canautocor <- function(X, preproc=center(), ncomp=2, tau=c(.2,.2)) {
  X1 <- X[1:nrow(X)-1,]
  X2 <- X[2:nrow(X),]
  
  fit <- rgcca(A= list(X1, X2), ncomp=rep(ncomp,2),
        C = matrix(c(0, 1, 1, 0), 2, 2),
        tau = tau)
  
  scores <- X %*% fit$a[[1]]
  ret <- list(X=X, scores=scores, fit=fit, ncomp=ncomp, tau=tau)
  class(ret) <- "canautocor"
  ret
}

scores.canautocor <- function(x) {
  x$X %*% x$fit$a[[1]]
}


loadings.canautocor <- function(x) {
  x$fit$a[[2]]
}
  