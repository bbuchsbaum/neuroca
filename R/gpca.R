
.msqrt <- function(a) {
  a.eig <- eigen(a)
  a.sqrt <- a.eig$vectors %*% diag(sqrt(a.eig$values)) %*% solve(a.eig$vectors)
}

#' genpca
#' 
#' @param X the data matrix
#' @param A the column constraints. Can be a \code{vector}, \code{matrix}, or sparse matrix.
#' @param M the row constraints. Can be a \code{vector}, \code{matrix}, or sparse matrix.
#' @param ncomp the number of components to estimate
#' @param center
#' @param scale
#' @importFrom assertthat assert_that
#' @export
genpca <- function(X, A, M, ncomp=min(dim(X)), center=FALSE, scale=FALSE) {
  
  if (is.vector(A)) {
    assert_that(length(A) == ncol(X))
    A <- sparseMatrix(i=1:length(A), j=1:length(A),x=A)
  }
  
  if (is.vector(M)) {
    assert_that(length(M) == nrow(X))
    M <- sparseMatrix(i=1:length(M), j=1:length(M),x=M)
  }
  
  assert_that(nrow(A) == ncol(A))
  assert_that(nrow(M) == ncol(M))
  assert_that(nrow(A) == ncol(X))
  assert_that(nrow(M) == nrow(X))
  assert_that(ncomp > 0)
  
  ncomp <- min(min(dim(X)), ncomp)
  
  n = nrow(X)
  p = ncol(X)
  
  if(n < p){
    ret = gmdLA(t(X), A,M, ncomp,p,n)
    svdfit = list(u=ret$v, v=ret$u,d=ret$d, cumv=ret$cumv,propv=ret$propv)
  }else{
    svdfit = gmdLA(X, M,A,ncomp,n,p)
  }
  
  scores <- t(t(as.matrix(svdfit$u)) * svdfit$d)
  
  ret <- list(v=svdfit$v, u=svdfit$u, d=svdfit$d, scores=scores, ncomp=ncomp)

  class(ret) <- c("gpca")
  ret
  
}

gmdLA <- function(X,Q,R,k,n,p){

  ##computation

  decomp <- eigen(R)
  keep <- which(abs(decomp$values) > 1e-7)
  
  decomp$vectors <- decomp$vectors[,keep]
  decomp$values <- decomp$values[keep]
  
  Rtilde = decomp$vectors %*% diag(sqrt(decomp$values)) %*% t(decomp$vectors)
  
  inv.values = 1/sqrt(decomp$values)
  Rtilde.inv = decomp$vectors %*% diag(inv.values) %*% t(decomp$vectors)
  
  #inmat =  t(X) %*% Q %*% X
  
  print(class(X))
  print(class(Q))
  inmat <- Matrix::crossprod(X, Q) %*% X

  RtinRt = Rtilde %*% inmat %*% Rtilde
  xtilde.decomp = eigen(RtinRt)
  
  keep <- which(abs(xtilde.decomp$values) > 1e-7)
  k <- length(keep)
  xtilde.decomp$vectors <- xtilde.decomp$vectors[,keep]
  xtilde.decomp$values <- xtilde.decomp$values[keep]
  
  vgmd <- Rtilde.inv %*% xtilde.decomp$vectors
  dgmd <- sqrt(xtilde.decomp$values[1:k])

  ugmd = matrix(nrow = n,ncol = k)
  cumv = rep(0,k)


  propv = dgmd^2/sum(diag(as.matrix(inmat %*% R)))

  normalizing.number = 1

  XR = X %*% R
  RnR = R %*% inmat %*% R
  
  k <- length(keep)

  for( i in 1:k) {
    print(i)
    normalizing.number = sqrt(vgmd[,i] %*% RnR %*% vgmd[,i])
    ugmd[,i] = as.matrix(XR %*% vgmd[,i])/as.double(normalizing.number)
    cumv[i] = sum(propv[1:i])
  }

  list(u=ugmd[,1:k],
       v=vgmd[,1:k],
       d=dgmd[1:k],
       cumv=cumv,
       propv=propv)
  
}
