
#' @keywords internal
unit_norm <- function(X) {
  Xs <- scale(X, center=TRUE, scale=TRUE)
  div <- sqrt(nrow(X) - 1)
  Xs/div
}


relief_scores <- function(X, labels, k=10) {
  labels <- as.factor(labels)
  
  X <- unit_norm(X)
  neighborweights:::label_matrix2(labels, labels)
  
  nclasses <- length(levels(labels))
  
  D <- Diagonal(nrow(X)) 
  knabes <- neighborweights::similarity_matrix(X, k=k, neighbor_mode="knn", weight_mode="binary", sigma=1) 
  lambat <- neighborweights:::label_matrix2(labels,labels)
  
  S <- lambat
  S[Matrix::which(knabes==1 & lambat == 1, arr.ind=TRUE)] <- 1/k
  S[Matrix::which(knabes==1 & lambat == 0, arr.ind=TRUE)] <- -1/((nclasses-1)*k)
  S[Matrix::which(knabes==0 & lambat == 1, arr.ind=TRUE)] <- 0
  diag(S) <- rep(1, nrow(S))
  
  D <- rowSums(S)
  L <- Diagonal(x=D) - S

  Dtilde <- Diagonal(x= D^(-1/2))
  Lnorm <- Dtilde %*% L %*% Dtilde
  X %*% Lnorm %*% X
  
}

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

