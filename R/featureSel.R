matrixAnova <- function(Y, mat) {
  K <- length(levels(Y))
  N <- length(Y)
  
  Vbetween <- colSums(do.call(rbind, lapply(levels(Y), function(lev) {
    idx <- which(Y == lev)
    G <- colMeans(mat[idx,])
    n <- length(idx)
    (G^2 * n)/(K-1)  
  })))
  
  Vwithin <- colSums(do.call(rbind, lapply(levels(Y), function(lev) {
    idx <- which(Y == lev)
    G <- colMeans(mat[idx,])
    (sweep(mat[idx,], 2, G)^2)/(N-K)
  })))
  
  F <- Vbetween/Vwithin
  pval <- pf(Vbetween/Vwithin, K-1, N-K, lower.tail=FALSE)   
  list(F=F, pval=pval)
}

SDAFeatureSel <- function(thresh=.5, criterion="fdr") {
  function(Y, mat) {
    if (!is.factor(Y)) {
      stop("argument is not a factor: SDA feature selection only performed on factors")
    }
    ranks <- sda.ranking(mat, Y)
    keep <- ranks[, "lfdr"] < thresh
    
    
    ret <- logical(ncol(mat))
    ret[ranks[keep, "idx"]] <- TRUE
    ret
  }
}

AnovaFeatureSel <- function(thresh=.05, criterion=c("pval", "F")) {
  function(Y, mat) {
    print(thresh)
    if (!is.factor(Y)) {
      stop("argument is not a factor: Anova feature selection only performed on factors")
    }
    scores <- matrixAnova(Y, mat)
    if (criterion == "pval") {
      scores$pval < thresh
    } else if (criterion == "F") {
      scores$F > thresh
    } else {
      stop(paste("unrecognized criterion: ", criterion))
    }
  }
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

