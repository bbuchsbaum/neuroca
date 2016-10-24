
.msqrt <- function(a) {
  a.eig <- eigen(a)
  a.sqrt <- a.eig$vectors %*% diag(sqrt(a.eig$values)) %*% solve(a.eig$vectors)
}

#' @param X the data matrix
#' @param A the column constraints. Can be a \code{vector}, \code{matrix}, or sparse matrix.
#' @param M the row constraints. Can be a \code{vector}, \code{matrix}, or sparse matrix.
#' @param ncomp the number of components to estimate
#' @param center
#' @param scale
#' @importFrom assertthat assert_that
#' @importFrom sGPCA gpca
#' @export
generalized_pca <- function(X, A, M, ncomp=min(dim(X)), center=FALSE, scale=FALSE) {
  
  if (is.vector(A)) {
    assert_that(length(A) == ncol(X))
    A <- sparseMatrix(i=1:length(A), j=1:length(A),x=A)
  }
  
  if (is.vector(M)) {
    assert_that(length(M) == nrow(X))
    M <- sparseMatrix(i=1:length(M), j=1:length(M),x=M)
  }
  
  gpfit <- sGPCA::gpca(X, M, A, K=10)
  scores <- t(t(as.matrix(gpfit$U)) * gpfit$D)
  
  ret <- list(v=gpfit$v, u=gpfit$u, d=gpfit$d, scores=scores, ncomp=ncomp)
  
  class(ret) <- c("pca_result")
  ret
  
}

# gmdLA <- function(X,Q,R,k,n,p){
#   
#   ##computation
#   
#   decomp = eigen(R)
#   Rtilde = decomp$vectors %*% diag(sqrt(decomp$values)) %*% t(decomp$vectors)
#   inv.values = decomp$values
#   inv.values[which(decomp$values!=0)] = 1/sqrt(decomp$values[which(decomp$values!=0)])
#   Rtilde.inv = decomp$vectors %*% diag(inv.values) %*% t(decomp$vectors)
#   inmat =  t(X) %*% Q %*% X
#   
#   RtinRt = Rtilde %*% inmat %*% Rtilde
#   xtilde.decomp = eigen(RtinRt)
#   
#   vgmd = Rtilde.inv %*% xtilde.decomp$vectors
#   dgmd = sqrt(xtilde.decomp$values[1:k])
#   
#   ugmd = matrix(nrow = n,ncol = k)
#   cumv = rep(0,k)
#   
#   
#   propv = dgmd^2/sum(diag(as.matrix(inmat %*% R)))
#   
#   normalizing.number = 1
#   
#   XR = X %*% R
#   RnR = R %*% inmat %*% R
#   
#   for( i in 1:k){
#     normalizing.number = sqrt(vgmd[,i] %*% RnR %*% vgmd[,i])
#     ugmd[,i] = as.matrix(XR %*% vgmd[,i])/as.double(normalizing.number)
#     cumv[i] = sum(propv[1:i]) 
#   }
#   
#   return.object = list("u"=ugmd,"v"=vgmd[,1:k],"d"=dgmd,"cumv"=cumv,"propv"=propv)
#   return(return.object)
#   
# }
