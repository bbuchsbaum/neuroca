

#' locality preserving projections
#' 
#' @param X the data matrix
#' @param W the affinity matrix
#' @param ndim number of dimensions
#' @importFrom RSpectra eigs
#' @export
embed_lpp <- function(X, W, ndim=2) {
  
  D <- Diagonal(x=rowSums(W))
  L <- D - W
  
  Zl <- t(X) %*% L %*% X
  Zr <- t(X) %*% D %*% X
  
  decomp <- RSpectra::eigs(solve(Zr,Zl), k=ndim, which="SM")
  
}


cut_pca <- function(U,D,V, ratio=.95) {
  stopifnot(ratio > 0 && ratio < 1)
  cutoff <- sum(D) * ratio
  k <- which(cumsum(D) < cutoff)
  k <- k[length(k)]
  
  U <- U[, 1:k]
  V <- V[,1:k]
  D <- D[1:k]
  
  list(d=D, u=U, v=V)
  
}


lge <- function(X, W, D=NULL, ndim=3, regularize=FALSE, alpha=.1, pca_ratio=1) {
  do_pca <- !regularize
  
  nsamples <- nrow(X)
  nfeatures <- ncol(X)
  
  if (nrow(W) != nsamples) stop('W and X mismatch!');
 
  if (!is.null(D) && (nrow(D) != nsamples)) {
    stop('D and data mismatch!')
  }
  
  #   bChol = 0;
  #   skip bChol
  
  if (do_pca) {   
    
    svd_full <- corpcor::fast.svd(X)
    svd_red  <- cut_pca(svd_full$u,svd_full$d,svd_full$v,pca_ratio)
    eigvalue_PCA <- diag(svd_red$d)
    
    if (!is.null(D)) {
      X <- svd_red$u %*% Diagonal(x=svd_red$d)
      eigvector_PCA <- svd_red$v
      DPrime <- t(X) %*% D %*% X
      ##DPrime = max(DPrime,DPrime');
    } else {
      X <- svd_red$u
      eigvector_PCA = svd_red$v %*% Matrix::Diagonal(x=eigvalue_PCA.^-1)
    }
  } else {
    if (!is.null(D)) {
      DPrime <- t(X) %*% D %*% X
    } else {
      DPrime = t(X) %*% X
    }
    
     DPrime <- DPrime + Diagonal(n=nrow(DPrime), x=alpha)
    }
    
    WPrime = t(X)  %*% W %*% X
    #WPrime = max(WPrime,WPrime');
    
    mdim = ncol(WPrime)
    ndim <- min(mdim, ndim)
    
    edecomp <- if (do_pca && is.null(D)) {
      RSpectra::eigs(WPrime,ndim,'LM')
    } else {
      #browser()
      RSpectra::eigs(solve(DPrime, WPrime),ndim,'LM')
      #eigen(solve(DPrime, WPrime))
    }
    
    eigvector <- edecomp$vectors[,1:ndim]
    
      
    if (do_pca) {
       eigvector = eigvector_PCA %*% eigvector
    }
    
    for (i in 1:ncol(eigvector)) {
      eigvector[,i] = eigvector[,i]/norm(eigvector[,i], "2")
    }
    
    list(vectors=eigvector, values=edecomp$values[1:ndim])
  
}

lsda <- function(X, labels, k=5, sigma=.7, beta=.1, regularize=FALSE, alpha=.1) {
  labels <- as.factor(labels)
  nlabel <- length(levels(labels))
  
  Ww <- label_matrix(labels, labels, type="s")
  Wb <- label_matrix(labels, labels, type="d")
  
  if (k > 0) {
    W <- construct_weight_matrix(X, neighbor_mode = "knn", weight_mode="normalized", k=k, sigma=sigma, labels=labels)
    Ww <- Ww * W
    Wb <- Wb * W
  }
  
  Db <- rowSums(Wb,2)
  Wb <- -Wb
  
  Wb <- Wb + Matrix::Diagonal(x=Db)
  
  D <- rowSums(Ww)
  
  if (regularize) {
    alpha <- alpha * sum(D)/length(D)
  }
  
  Wcombined = (beta/(1-beta))*Wb + Ww 
  D <- Diagonal(x=D)
  
  v=lge(X, Wcombined, D, regularize=regularize, alpha=alpha, pca_ratio=pca_ratio)
  

}

repmat <- function(a,n,m) kronecker(matrix(1,n,m),a) 

