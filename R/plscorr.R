
#group_means1 <- function(Y, X) {
# G <- model.matrix(~ Y - 1)
# colnames(G) <- levels(Y)
# GW <- G/colSums(G)
# R <- t(crossprod(X, GW))
 #centroid <-  rowMeans(R)
 #Rcent <- sweep(R, 1, centroid)
 #list(Rcent=Rcent, Ymat=G, centroid=centroid)
#}\\


scores <- function(x) UseMethod("scores")
loadings <- function(x) UseMethod("loadings")



svd.wrapper <- function(XC, ncomp=min(dim(XC)), method=c("base", "fast", "irbla")) {
  res <- switch(method[1],
                base=svd(XC),
                fast.svd=corpcor:::fast.svd(XC),
                irlba=irlba:::irlba(XC, nu=min(ncomp, min(dim(XC)) -3), nv=min(ncomp, min(dim(XC)) -3)))
  
 
  res$d <- res$d[1:ncomp]
  res$u <- res$u[,1:ncomp]
  res$v <- res$v[,1:ncomp]
  res$ncomp <- ncomp
  res
}

group_means <- function(Y, X) {
  Rs <- rowsum(X,Y)
  yt <- table(Y)
  sweep(Rs, 1, yt, "/")
}


plscorr <- function(Y, X, ncomp=2, center=TRUE, scale=FALSE, svd.method="fast.svd") {
  if (is.factor(Y)) {
    plscorr_f(Y, X, ncomp, svd.method)
  }
}

scores.pls_result_f <- function(x) {
  x$scores
}

predict.plscorr_result_f <- function(x, newdata, type=c("class", "prob", "scores"), ncomp=NULL) {
  if (is.vector(newdata) && (length(newdata) == ncol(x$condMeans))) {
    x <- matrix(newdata, nrow=1, ncol=ncol(x$condMeans))
  }
  
  assert_that(is.matrix(newdata))
  assert_that(ncol(newdata) == ncol(x$condMeans))
  assert_that(type %in% c("class", "prob", "scores"))
  
  if (is.null(ncomp)) {
    ncomp <- x$ncomp
  }
  

  xc <- x$pre_process(newdata)
  fscores <- xc %*% x$v[,1:ncomp,drop=FALSE]
  
  if (type == "scores") {
    fscores
  } else if (type =="class") {
    D <- rdist(fscores[,1:ncomp,drop=FALSE], x$scores[,1:ncomp,drop=FALSE])
    D2 <- D^2
    min.d <- apply(D2, 1, which.min)
    levels(x$Y)[min.d]
  } else {
    ## type is 'prob'
    probs <- tcrossprod(fscores[,1:ncomp,drop=FALSE], x$scores[,1:ncomp,drop=FALSE])
    maxid <- apply(probs,1,which.max)
    maxp <- probs[cbind(1:nrow(probs), maxid)]
    probs <- exp(probs - maxp)
    probs <- zapsmall(probs/rowSums(probs))
  }
  
}


apply_scaling <- function(Xc) {
  force(Xc)
  function(M) {
    center <- !is.null(attr(Xc, "scaled:center"))
    scale <- !is.null(attr(Xc, "scaled:scale"))
    
    if (!center && !scale) {
      M
    } else if (center && !scale) {
      sweep(M, 2, attr(Xc, "scaled:center"))
    } else if (scale && !center) {
      sweep(M, 2, attr(Xc, "scaled:scale"), "/")
    } else {
      M <- sweep(M, 2, attr(Xc, "scaled:center"))
      sweep(M, 2, attr(Xc, "scaled:scale"), "/")
    }
  }
}

plscorr_f <- function(Y, X, ncomp=2, center=TRUE, scale=FALSE, svd.method="base") {
  assert_that(is.factor(Y))
  Xc <- scale(X, center=center, scale=scale)
  XB <- group_means(Y, Xc)
  svdres <- svd.wrapper(XB, ncomp, svd.method)
  
 

  scores <- svdres$u %*% diag(svdres$d, nrow=ncomp, ncol=ncomp)
  row.names(scores) <- levels(Y)

  ret <- list(Y=Y,ncomp=ncomp, condMeans=XB, center=center, scale=scale, pre_process=apply_scaling(Xc), 
              svd.method=svd.method, scores=scores, v=svdres$v, u=svdres$u, d=svdres$d)
  
  class(ret) <- "plscorr_result_f"
  ret
  
}

