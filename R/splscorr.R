normv <- function(x) sqrt(sum(x**2))

## TODO generic ALS method?
svd_iter <- function(X, ncomp=1,svd.method="irlba", cardinality=rep(ncol(X), ncomp), tol=1e-05) {
  
  Xcur <- X
  for (i in 1:ncomp) {
    svdres <- svd.wrapper(Xcur, ncomp=1, method=svd.method)
    u.new <- svdres$u[,1]
    v.new <- svdres$v[,1]
    iter <- 0
    v.stab <- FALSE
    u.stab <- FALSE
    
    nx <- ncol(X) - cardinality[i]
    
    while ((v.stab == FALSE) || (u.stab == FALSE) && iter < iter.max) {
      iter <- iter + 1
      print(paste("iter", iter))
      u.old <- u.new
      v.temp <- t(Xcur) %*% u.old
      v.old <- v.new
    
      if (cardinality[i] < ncol(X)) {
        
        thresh <- abs(sort(abs(v.temp), decreasing=TRUE)[cardinality[i]])
        print("thresh is", thresh)
        v.new <- ifelse(abs(v.temp) > thresh, (abs(v.temp) - thresh) * sign(v.temp), 0)
        v.new <- v.new/normv(v.new)
      }
      
      u.new <- (Xcur %*% v.new)[,1]
      u.new <- u.new/normv(u.new)
      
      print(paste("u", crossprod(u.new-u.old)))
      print(paste("v ", crossprod(v.new-v.old)))
      if(crossprod(u.new-u.old)<tol){ u.stab=TRUE }
      if(crossprod(v.new-v.old)<tol){ v.stab=TRUE }
      
    }
    
    ## deflation
    Xcur <- Xcur - svdres$d[1] * svdres$u[,1] %*% t(svdres$v[,1])
  }
  
}


#' @export
splscorr_da <- function(Y, X, ncomp=2, center=TRUE, scale=FALSE) {
  assert_that(is.factor(Y))
  
  if (any(table(Y) > 1)) {
    ## compute barycenters
    XB <- group_means(Y, X)
    XBc <- scale(XB, center=center, scale=scale)
  } else {
    XB <- X
    XBc <- scale(XB, center=center, scale=scale)
  }
  
  #svdres <- svd.wrapper(XBc, ncomp, svd.method)
  svdres <- svd.wrapper(t(XBc), ncomp, svd.method)
  
  scores <- svdres$v %*% diag(svdres$d, nrow=svdres$ncomp, ncol=svdres$ncomp)
  row.names(scores) <- levels(Y)
  
  ret <- list(Y=Y,X=X,ncomp=svdres$ncomp, condMeans=XBc, center=center, scale=scale, pre_process=apply_scaling(XBc), 
              svd.method=svd.method, scores=scores, v=svdres$v, u=svdres$u, d=svdres$d)
  
  class(ret) <- "plscorr_result_da"
  ret
  
}