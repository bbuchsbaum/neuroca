

trace_ratio <- function(A, B, ncomp) {
  Y <- iris$Species
  A <- between_class_scatter(X,Y, colMeans(X))
  B <- within_class_scatter(X, Y)
  
  ncomp = 2
  V <- matrix(rnorm(ncol(A)*ncomp), ncol(A), ncomp)
  V.qr <- qr(V)
  V <- qr.Q(V.qr)
  num <- (t(V) %*% A %*% V)
  denom <- (t(V) %*% B %*% V)
  rho <- sum(diag(num))/sum(diag(denom))
  print(rho)
  for (i in 1:10) {
    G <- A - rho*B
    V <- eigen(G)$vectors[,1:ncomp]
    num <- (t(V) %*% A %*% V)
    denom <- (t(V) %*% B %*% V)
    rho <- sum(diag(num))/sum(diag(denom))
    print(rho)
  }
}

ulda <- function(Y, X, preproc) {
  Y <- as.factor(Y)
  
  procres <- prep(preproc)
  Xp <- procres$init(X)
  
  levs <- levels(Y)
  nc <- length(levs)
  
  cprobs <- table(Y)/length(Y)
  cmeans <- group_means(Y,Xp)
  gmean <- colMeans(Xp)
  
  Hb <- do.call(cbind, lapply(1:nc, function(i) {
    sqrt(cprobs[i]) * (cmeans[i,] - gmean)
  }))
  
  Ht <- t(sweep(Xp, 2, gmean, "-"))
  t <- min(dim(Ht))
  
  svd_ht <- svd(Ht)
  B <- diag(1/svd_ht$d[1:t]) %*% t(svd_ht$u[,1:t]) %*% Hb
  
  svd_B <- svd(B)
  q <- min(dim(B))
  vecs <- svd_ht$u[,1:t,drop=FALSE] %*% diag(1/svd_ht$d[1:t,drop=FALSE]) %*% svd_B$u
  
  projector(procres, ncomp=t, v=vecs, class=c("ulda"))
}



lda_outer <- function (X) {
  p = ncol(X)
  output = array(0, c(p, p))
  for (i in 1:nrow(X)) {
    output = output + outer(X[i, ], X[i, ])
  }
  return(output)
}

## A Multiple Maximum Scatter Difference
## Discriminant Criterion for Facial Feature Extraction
## Fengxi Song, David Zhang, Senior Member, IEEE, Dayong Mei, and Zhongwei Guo
## https://paperpile.com/app/p/c1daa0f0-754f-08db-81e1-bd5484e4c675
mmsd <- function(X, Y, ncomp=min(dim(X)), scale=FALSE, c=1) {
  Y <- as.factor(Y)
  preproc <- pre_processor(X, center=FALSE, scale)
  Xp <- pre_process(preproc, X)
  
  pcred <- pca(Xp, center=TRUE)
  xscores <- scores(pcred)
  
  levs <- levels(Y)
  ncount <- table(Y)
  
  p <- ncol(xscores)
  Sw <- array(0, c(p, p))
  for (i in 1:length(levs)) {
    idxnow <- which(Y == levs[i])
    gmean <- colMeans(xscores[idxnow,])
    xs <- sweep(xscores[idxnow,], 2, gmean, "-")
    Sw <- Sw + lda_outer(xs)
  }
  
  
  Sb = array(0, c(p, p))
  m <- colMeans(xscores)
  for (i in 1:length(levs)) {
    idxnow = which(Y == levs[i])
    Nk = length(idxnow)
    mdiff = (colMeans(xscores[idxnow, ]) - m)
    Sb = Sb + Nk * outer(mdiff, mdiff)
  }
  
  SwSb <- Sb - c*Sw
  
  P <- loadings(pcred)
  decomp <- RSpectra::eigs_sym(P %*% SwSb %*% t(P), k=ncomp)
  
  posidx <- which(decomp$values > 1e-6)
  ncomp <- length(posidx)
  
  v <- t(P) %*% decomp$vectors[,posidx]
  scores <- Xp %*% v
  u <- apply(scores, 2, function(x) x/sqrt(sum(x^2)))
  
  fit <- bi_projector(preproc=preproc,
                      ncomp=ncomp,
                      u=u,
                      v=v,
                      d=apply(scores, 2, function(x) sum(x^2)),
                      scores=scores,
                      classes="mmsd",
                      c=c)
  fit$Y <- Y
  fit
  
}

#' mmpca
#' 
#' @inheritParams bada
#' @param knn the number of nearest neighbors used to form margin
#' @export
#' @inheritsParams pca
mmpca <- function(X, Y, ncomp=min(dim(X)), preproc=standardize(), knn=1, sigma=.7) {
  assert_that(nrow(X) == length(Y))
  Y <- as.factor(Y)
  
  preproc <- prep(preprox)
  Xp <- preproc$init(X)
  
  nabes <- neighborweights::graph_weights(Xp, "knearest_misses", k=knn, labels=Y, type="asym", sigma=sigma)
  nnbc <- apply(nabes, 1, function(x) which(x > 0))
  
  if (is.vector(nnbc)) {
    nnbc <- matrix(nnbc, ncol=nrow(X))
  }
  
  levs <- levels(Y)
  
  Xt <- do.call(rbind, lapply(1:nrow(Xp), function(i) {
    if (knn > 1) {
      Xp[i,] - colMeans(Xp[nnbc[,i],])
    } else {
      Xp[i,] - Xp[nnbc[,i],]
    }
  }))
  
  fit <- pca(Xt, center=FALSE, ncomp=ncomp)
  fit$Y <- Y
  class(fit) <- c("mmpca", class(fit))
  fit
  
}


# testrun <- function(knn=1, iter=25, ncomp=2) {
#   
#     tridx <- sort(unlist(lapply(levels(Y), function(lev) {
#       idx <- which(Y == lev)
#       tmp <- which(diff(idx) == 1)
#       idx[c(tmp, tmp+1)]
#     })))
#   
#     testidx <- seq(1,360)[-tridx]
#     fit <- mmpca(X[tridx,], Y[tridx], knn=knn, ncomp=ncomp)
#     #fit <- pca(X[tridx,], ncomp=ncomp)
#     fscores <- project(fit, X[testidx,])
#   
#     sfit <- scores(fit)
#     row.names(sfit) <- as.character(Y[tridx])
#     pred <- scorepred(fscores, sfit, type="class")
#     sum(pred == as.character(Y[testidx]))/length(testidx)
# }
# 
# grid <- expand.grid(ncomp=1:4, knn=1:24)
# grid$perf <- unlist(lapply(1:nrow(grid), function(i) {
#   testrun(grid$knn[i], iter=0, grid$ncomp[i])
# }))




