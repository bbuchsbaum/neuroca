spatial_mfa <- function(X, ncomp=2, center=TRUE, scale=FALSE, 
                        normalization=c("MFA", "RV", "None", "RV-MFA")) {
  
  blocklens <- block_lengths(X)
  assert_that(all(blocklens == blocklens[1]))
  assert_that(nrow(coord_mat) == blocklens[1])
  
  normalization <- match.arg(normalization)
  
  ## pre-process the projected variables.
  preproc <- pre_processor(X, center=center,scale=scale)
  Xr <- pre_process(preproc)
  
  ## normalize the matrices 
  #alpha <- normalization_factors(Xr, type=normalization)

  bind <- block_index_list(X)
  
  colind <- unlist(lapply(bind, function(i) seq_along(i)))
  Smat <- neighborweights::spatial_laplacian(as.matrix(colind), dthresh=1, weight_mode="binary", normalized=FALSE)
  
  
}



empca <- function(x, k=4, alpha=1, Smat, thresh=1e-4) {
  U <- matrix(rnorm(nrow(x)*k), nrow(x), k)
  A <- t(matrix(rnorm(ncol(x)*k), k, ncol(x)))
  
  crit = 1
  
  while(crit > thresh) {
    print(crit)
    Uinv <- corpcor::pseudoinverse(U)
    
    A <- t(Uinv %*% x)
    Aprime <- -Smat %*% A
    
    #A <- A + alpha*Aprime
    A <- A + (alpha * Aprime)
    #A <- Aprime
    Ui <- x %*% t(corpcor::pseudoinverse(as.matrix(A)))
    #Ui <- sweep(U, 2, FUN="/", STATS=apply(Ui, 2, function(x) sqrt(sum(x ^2))))
    
    d2 <- sum(abs(Ui - U))
    crit <- d2
    
    U <- Ui
    
    
  }
  
  U <- svd(U)$u
  A <- t(U) %*% x
  
  
  
}
   

        
    
    
    
  
  
}