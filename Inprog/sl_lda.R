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

##A general soft label based Linear Discriminant Analysis for
##semi-supervised dimensionality reduction
## https://paperpile.com/app/p/03eb8cb3-2326-0cc5-a440-02c6aa543bf3 

## see also, maybe similar:
## A weighted linear discriminant analysis framework for multi-label
## feature extraction

soft_lda <- function(C, X, preproc=pass(), dp=min(dim(X)), di=dp, dl=ncol(C)-1) {
  assert_that(nrow(C) == nrow(X))
  
  if (is.null(colnames(C))) {
    colnames(C) <- paste0("c_", 1:ncol(C))
  }
  
  assert_that(all(C >= 0), msg="all weights in 'C' must be positive")
  
  ## pre-process X 
  procres <- prep(preproc, X)
  Xp <- procres$Xp
  
  ## reduce Xp into dp dimensions
  pca_red <- pca(Xp, ncomp=dp, preproc=center())
  
  ## get projection -- assumes Xp is centered...
  proj_dp <- loadings(pca_red)
  
  ## get pca basis
  Xpca <- scores(pca_red)
  
  ## row sums of soft-label matrix
  E <- diag(rowSums(C))
  
  ## column sums of soft-label matrix
  G <- diag(colSums(C))
  
  ## transposed weight matrix: each row is a weight vector for a class
  F <- t(as.matrix(C))
 
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
    M <- FtGF - num/denom
    Xt %*% M %*% t(Xt)
  }
  
  
  gmeans <- weighted_group_means(Xpca, F)
  mu <- colMeans(gmeans)
  Sw <- sw_scatter(Xpca)
  Sb <- sb_scatter(Xpca)
  
  E_i <- RSpectra::eigs_sym(Sw, k=di)
  
  proj_di <- E_i$vectors %*% diag(1/sqrt(E_i$values))
  
  gmeans_proj <- (gmeans %*% proj_di)
  
  E_l <- pca(gmeans_proj, ncomp=dl)
  
  proj_dl <- loadings(E_l)
  
  proj_final <- proj_dp %*% proj_di %*% proj_dl
  
  projector(procres, ncomp=ncol(proj_final), v=proj_final, classes="sl_lda")
  
  #St <- st_scatter(X, F, mu)
  
}
  
  
  

# 
# 
# 
# sw_scatter_naive <- function(X, F, gmeans) {
#   Reduce("+", lapply(1:nrow(gmeans), function(i) {
#     Reduce("+", lapply(1:nrow(X), function(j) {
#       d <- X[j,] - gmeans[i,]
#       F[i,j] * outer(d,d)
#     }))
#   }))
# }
# 
# sb_scatter_naive <- function(X, F, gmeans, mu) {
#   Reduce("+", lapply(1:nrow(gmeans), function(i) {
#     d <- gmeans[i,] - mu
#     Reduce("+", lapply(1:nrow(X), function(j) {
#       F[i,j] * outer(d,d)
#     }))
#   }))
# }
# 
# st_scatter_naive <- function(X, F, mu) {
#   Reduce("+", lapply(1:nrow(F), function(i) {
#     Reduce("+", lapply(1:nrow(X), function(j) {
#       d <- X[j,] - mu
#       F[i,j] * outer(d,d)
#     }))
#   }))
# }


