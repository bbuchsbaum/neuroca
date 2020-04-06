
#' trace_ratio optimization
#' 
#' @param A the numerator matrix
#' @param B the denominator matrix
#' 
#' @keywords 
#' 
#' 2012 Ngo et al 
#' @export
trace_ratio <- function(A, B, ncomp, eps=1e-6, maxiter=100) {
  ## in the the language
  n = nrow(A)
  p = ncomp
  
  ## prepare the initializer
  Vold = qr.Q(qr(matrix(rnorm(n*p),ncol=p)))
  rhoold = 0
  for (i in 1:maxiter){
    Vnew   = RSpectra::eigs(A-rhoold*B,p,which="LR")$vectors
    rhonew = sum(diag(t(Vnew)%*%A%*%Vnew))/sum(diag(t(Vnew)%*%B%*%Vnew))
    
    rhoinc = abs(rhonew-rhoold)
    Vold   = Vnew
    rhoold = rhonew
    
    if (rhoinc < eps){
      break
    }
  }
  
  V <- Vold
  
  values <- sapply(1:ncomp, function(i) {
    (V[,i] %*% A %*% V[,i])/ (V[,i] %*% A %*% V[,i])
  })
  
  #values = diag(t(V)%*%A%*%V)/diag(t(V)%*%B%*%V)
  list(vectors=V,values=values)
}