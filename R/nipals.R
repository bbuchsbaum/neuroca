

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