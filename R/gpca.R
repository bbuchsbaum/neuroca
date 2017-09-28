
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
#' @param center whether to center the columns
#' @param scale whether to standardize the columns
#' @importFrom assertthat assert_that
#' @export
genpca <- function(X, A=rep(1:ncol(X)), M=rep(1,nrow(X)), ncomp=min(dim(X)), center=FALSE, scale=FALSE) {
  
  #if (!is.vector(A) && isDiagonal(A)) {
  #  A <- diag(A) 
  #}
  
  #if (!is.vector(M) && isDiagonal(M)) {
  #  M <- diag(M) 
  #}
  
  if (is.vector(A)) {
    assert_that(length(A) == ncol(X))
    A <- sparseMatrix(i=1:length(A), j=1:length(A),x=A)
  } else {
    assert_that(nrow(A) == ncol(A))
    assert_that(nrow(A) == ncol(X))
  }
  
  if (is.vector(M)) {
    assert_that(length(M) == nrow(X))
    M <- sparseMatrix(i=1:length(M), j=1:length(M),x=M)
  } else {
    assert_that(nrow(M) == ncol(M))
    assert_that(nrow(M) == nrow(X))
  }
  
  assert_that(ncomp > 0)
  ncomp <- min(min(dim(X)), ncomp)
  
  n = nrow(X)
  p = ncol(X)
  
  X <- pre_processor(X, center=center,scale=scale)
  
  if(n < p){
    ret = gmdLA(t(X), A,M, ncomp,p,n)
    svdfit = list(u=ret$v, v=ret$u,d=ret$d, cumv=ret$cumv,propv=ret$propv)
  }else{
    svdfit = gmdLA(X, M,A,ncomp,n,p)
  }
  
  scores <- t(t(as.matrix(svdfit$u)) * svdfit$d)
  #norm_loadings <- t(t(as.matrix(svdfit$v)) * svdfit$d)
  
  ret <- list(v=svdfit$v, 
              u=svdfit$u, 
              d=svdfit$d, 
              eig=svdfit$d^2,
              scores=scores, 
              ncomp=length(svdfit$d), 
              A=A,
              M=M,
              pre_process=attr(X, "pre_process"), 
              reverse_pre_process=attr(X, "reverse"))

  class(ret) <- c("genpca", "list")
  ret
  
}

project_xav <- function(X, A, V) {
  if (is.vector(A)) {
    t(t(X) * A) %*% V
  } else {
    X %*% A %*% V
  }
}
  
  
  
predict.genpca <- function(x, newdata, ncomp=x$ncomp, pre_process=TRUE) {
  Xsup <- if (pre_process) {
    x$pre_process(newdata)
  } else {
    newdata
  }
  
  project_xav(Xsup, x$A, x$v[,1:ncomp,drop=FALSE])
  
}

contributions.genpca <- function(x) {
  coord.var <- t(t(as.matrix(x$v)) * x$d)
  t(t(coord.var^2) %*% x$A)
}

project.genpca <- function(x, newdata, ncomp=x$ncomp, pre_process=TRUE, subind=NULL) {
  if (is.null(subind)) {
    predict(x, newdata, ncomp, pre_process)
  } else {
    assertthat::assert_that(length(subind) == ncol(newdata))
    Xsup <- if (pre_process) {
      x$pre_process(newdata, subind)
    } else {
      newdata
    }
    project_xav(Xsup, x$A[subind,subind], x$v[,1:ncomp,drop=FALSE])
  }
}

# mmult <- function(X, q) {
#   if (is.vector(q)) {
#     t(t(X) * q)
#   } else {
#     X %*% q
#   }
# }
# 
# mmult2 <- function(q, X) {
#   if (is.vector(q)) {
#     t(q * t(X))
#   } else {
#     q %*% X
#   }
#   
# }
# 
# cprod <- function(X, q) {
#   if (is.vector(q)) {
#     t(X * q)
#   } else {
#     Matrix::crossprod(X, q)
#   }
# }

gmdLA <- function(X, Q, R, k, n, p) {
  ##computation
  
  
  decomp <- eigen(R)
  keep <- which(abs(decomp$values) > 1e-7)
  
  decomp$vectors <- decomp$vectors[, keep]
  decomp$values <- decomp$values[keep]
  
  Rtilde <-
    decomp$vectors %*% diag(sqrt(decomp$values)) %*% t(decomp$vectors)
  
  inv.values = 1 / sqrt(decomp$values)
  Rtilde.inv = decomp$vectors %*% diag(inv.values) %*% t(decomp$vectors)
  inmat <- crossprod(X, Q) %*% X
  
  RtinRt <- Rtilde %*% inmat %*% Rtilde
  
  XR <- X %*% R
  RnR <- R %*% inmat %*% R
  

  xtilde.decomp <- eigen(RtinRt)
  keep <- which(abs(xtilde.decomp$values) > 1e-7)
  k <- length(keep)
  xtilde.decomp$vectors <- xtilde.decomp$vectors[, keep]
  xtilde.decomp$values <- xtilde.decomp$values[keep]
  
  Rtilde.inv %*% xtilde.decomp$vectors
  
  vgmd <- Rtilde.inv %*% xtilde.decomp$vectors
  dgmd <- sqrt(xtilde.decomp$values[1:k])
  ugmd <- matrix(nrow = n, ncol = k)
  cumv <- rep(0, k)
  propv <- dgmd ^ 2 / sum(diag(as.matrix(inmat %*% R)))
  normalizing.number <- 1
  
  for (i in 1:k) {
    normalizing.number = sqrt(vgmd[, i] %*% RnR %*% vgmd[, i])
    ugmd[, i] = as.matrix(XR %*% vgmd[, i]) / as.double(normalizing.number)
    cumv[i] = sum(propv[1:i])
  }
  
  list(
    u = ugmd[, 1:k, drop=FALSE],
    v = vgmd[, 1:k, drop=FALSE],
    d = dgmd[1:k],
    k = k,
    cumv = cumv,
    propv = propv
  )
  
}
