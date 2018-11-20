library(fmrireg)

N <- 27
fac <- factor(sapply(1:9, function(i) sample(letters[1:3], 3)))
onsets <- cumsum(rnorm(N, m=3.9,sd=1.7))
df1 <- data.frame(fac=fac, onsets=onsets,run=1)


sf <- sampling_frame(round(max(onsets))+5, TR=1)
ev <- event_model(onsets ~ hrf(fac, basis="gaussian"), block= ~ run, 
                  sampling_frame=sf, data=df1)



dmat <- as.matrix(design_matrix(ev))
X <- matrix(rnorm(nrow(dmat)*100), nrow(dmat), 100)


weighted_group_means <- function(X, F) {
  ret <- do.call(rbind, lapply(1:nrow(F), function(i) {
    w <- F[i,]
    matrixStats::colWeightedMeans(X, w/sum(w))
  }))
  
  row.names(ret) <- row.names(F)
  ret
}


soft_lda <- function(C, X, preproc=pass(), dp=min(dim(X)), di=dp, dl=ncol(C)-1) {
  assert_that(nrow(C) == nrow(X))
  if (is.null(colnames(C))) {
    colnames(C) <- paste0("c_", 1:ncol(C))
  }
  
  assert_that(all(C >= 0), msg="all weights in 'C' must be positive")

  E <- diag(rowSums(C))
  G <- diag(colSums(C))
  F <- t(as.matrix(C))
  #e <- rep(1, nrow(dmat))

  FtGF <- (t(F) %*% diag(1/diag(G)) %*% F)
  
  sw_scatter <- function(X) {
    Xt <- t(X)
    Xt %*% (E - FtGF) %*% t(Xt)
  }
  
  sb_scatter <- function(X) {
    Xt <- t(X)
    e <- matrix(e)
    #num <- E %*% e %*% t(e) %*% E
    num <- tcrossprod(diag(E), rep(1,ncol(E))) %*% E
    #denom <- t(e) %*% E %*% e
    denom <- sum(diag(E))
    M <- FtGF - num/denom[,]
    Xt %*% M %*% t(Xt)
  }
  
  
  gmeans <- weighted_group_means(X, F)
  mu <- colMeans(gmeans)
  Sw <- sw_scatter(X, F, gmeans)
  Sb <- sb_scatter(X, F, gmeans,mu)
  
  E_i <- RSpectra::eigs_sym(Sw, k=di)
  
  proj_di <- E_i$vectors %*% diag(1/sqrt(E_i$values))
  
  gmeans_proj <- (gmeans %*% proj_di)
  
  E_l <- pca(gmeans_proj, ncomp=dl)
  
  proj_dl <- loadings(E_l)
  
  proj_final <- proj_dp %*% proj_di %*% proj_dl
  
  
  #St <- st_scatter(X, F, mu)
  
}
  
  
  
  



  




sw_scatter_naive <- function(X, F, gmeans) {
  Reduce("+", lapply(1:nrow(gmeans), function(i) {
    Reduce("+", lapply(1:nrow(X), function(j) {
      d <- X[j,] - gmeans[i,]
      F[i,j] * outer(d,d)
    }))
  }))
}

sb_scatter_naive <- function(X, F, gmeans, mu) {
  Reduce("+", lapply(1:nrow(gmeans), function(i) {
    d <- gmeans[i,] - mu
    Reduce("+", lapply(1:nrow(X), function(j) {
      F[i,j] * outer(d,d)
    }))
  }))
}

st_scatter_naive <- function(X, F, mu) {
  Reduce("+", lapply(1:nrow(F), function(i) {
    Reduce("+", lapply(1:nrow(X), function(j) {
      d <- X[j,] - mu
      F[i,j] * outer(d,d)
    }))
  }))
}


