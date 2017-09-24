
laplacian_scores <- function(X, W) {
  nsamples <- nrow(X)
  nfeatures <- ncol(X)

  if (nrow(W) != nsamples) {
    stop("nrow(W) must match nrow(X)")
  }
  
  D = Diagonal(x=rowSums(W))
  L = W
  
  tmp1 <- diag(D) %*% X
  #tmp1 = t(D) %*% X
  
  #allone <- rep(1, nsamples)
  #Xr <- sweep(X, 2, as.vector((t(X) %*% D %*% allone)/sum(diag(D))), FUN="-")
  
  DPrime <- colSums(t(t(X) %*% D) * X) - tmp1 * tmp1/sum(diag(D))
  LPrime <- colSums(t(t(X) %*% L) * X) - tmp1 * tmp1/sum(diag(D))
  #DPrime = sum((X'*D)'.*X)-tmp1.*tmp1/sum(diag(D));
  #LPrime = sum((X'*L)'.*X)-tmp1.*tmp1/sum(diag(D));
                            
  DPrime[which(DPrime) < 1e-12] = 10000;

  Y = t(LPrime/DPrime)
  #Y = full(Y);
}

