
#group_means1 <- function(Y, X) {
# G <- model.matrix(~ Y - 1)
# colnames(G) <- levels(Y)
# GW <- G/colSums(G)
# R <- t(crossprod(X, GW))
 #centroid <-  rowMeans(R)
 #Rcent <- sweep(R, 1, centroid)
 #list(Rcent=Rcent, Ymat=G, centroid=centroid)
#}\\

combinedACC <- function(Pred, Obs) {
  levs <- levels(as.factor(Obs))
  maxind <- apply(Pred, 1, which.max)
  pclass <- levs[maxind]
  sum(pclass == Obs)/length(pclass)
  
}

combinedAUC <- function(Pred, Obs) {
  mean(sapply(1:ncol(Pred), function(i) {
    lev <- levels(Obs)[i]
    pos <- Obs == lev
    pclass <- Pred[,i]
    pother <- rowMeans(Pred[,-i,drop=FALSE])
    Metrics::auc(as.numeric(pos), pclass - pother)-.5
  }))
}


#' @export
scores <- function(x) UseMethod("scores")

#' @export
loadings <- function(x) UseMethod("loadings")

#' @export
cross_validate <- function(x, ...) UseMethod("cross_validate")


#' @export
svd.wrapper <- function(XC, ncomp=min(dim(XC)), method=c("base", "fast", "irbla")) {
  res <- switch(method[1],
                base=svd(XC),
                fast.svd=corpcor:::fast.svd(XC),
                irlba=irlba:::irlba(XC, nu=min(ncomp, min(dim(XC)) -3), nv=min(ncomp, min(dim(XC)) -3)))
  
 
  res$d <- res$d[1:length(res$d)]
  res$u <- res$u[,1:length(res$d)]
  res$v <- res$v[,1:length(res$d)]
  res$ncomp <- length(res$d)
  res
}

#' @export
group_means <- function(Y, X) {
  Rs <- rowsum(X,Y)
  yt <- table(Y)
  sweep(Rs, 1, yt, "/")
}

#' @export
plscorr <- function(Y, X, ncomp=2, center=TRUE, scale=FALSE, svd.method="fast.svd") {
  if (is.factor(Y)) {
    plscorr_f(Y, X, ncomp, svd.method)
  }
}

#' @export
scores.pls_result_f <- function(x) {
  x$scores
}

nested_cv <- function(x, innerFolds, heldout, metric="AUC", min.comp=1) {
  res <- lapply(innerFolds, function(fidx) {
    exclude <- sort(c(fidx,heldout))
    prescv <- plscorr_f(x$Y[-exclude], x$X[-exclude,,drop=FALSE], ncomp=x$ncomp, 
                      center=x$center, scale=x$scale, svd.method=x$svd.method)
    
    pscores <- lapply(seq(min.comp, x$ncomp), function(n) {
      predict(prescv, x$X[fidx,,drop=FALSE], ncomp=n, type="crossprod")
    })
  })
  
  unlist(lapply(seq(min.comp, x$ncomp), function(n) {
    Pmat <- do.call(rbind, lapply(res, "[[", n))
    combinedAUC(Pmat, x$Y[-heldout])
  }))
  
}
  
#' @export
cross_validate.plscorr_result_f <- function(x, folds) {
  if (length(folds) == 1) {
    folds <- caret::createFolds(1:length(x$Y), folds)
  } else if (length(folds) == length(x$Y)) {
    folds <- split(1:length(x$Y), folds)
  }
  
  preds <- lapply(seq_along(folds), function(fnum) {
    print(fnum)
    res <- nested_cv(x, folds[-fnum], folds[[fnum]], min.comp=1)
    nc <- which.max(res)
    
    fidx <- folds[[fnum]]
    pmod <- plscorr_f(x$Y[-fidx], x$X[-fidx,,drop=FALSE], ncomp=nc, svd.method=x$svd.method) 
    pclass <- predict(pmod, x$X[folds[[fnum]],,drop=FALSE], ncomp=nc, type="class")
    prob <- predict(pmod, x$X[folds[[fnum]],,drop=FALSE], ncomp=nc, type="prob")
    list(class=pclass,prob=prob,K_perf=res, K_best=nc)
  })
  
  probmat <- do.call(rbind, lapply(preds, "[[", "prob"))
  class <- unlist(lapply(preds, "[[", "class"))
  ncomp <- unlist(lapply(preds, "[[", "K_best"))
  perfmat <- do.call(rbind, lapply(preds, "[[", "K_perf"))
  list(class=class, prob=probmat, ncomp=ncomp, ncomp_perf=perfmat)
}

#' @export
predict.plscorr_result_f <- function(x, newdata, type=c("class", "prob", "scores", "crossprod"), ncomp=NULL) {
  if (is.vector(newdata) && (length(newdata) == ncol(x$condMeans))) {
    x <- matrix(newdata, nrow=1, ncol=ncol(x$condMeans))
  }
  
  assert_that(is.matrix(newdata))
  assert_that(ncol(newdata) == ncol(x$condMeans))
  assert_that(type %in% c("class", "prob", "scores", "crossprod"))
  
  if (is.null(ncomp)) {
    ncomp <- x$ncomp
  }
  

  xc <- x$pre_process(newdata)
  fscores <- xc %*% x$v[,1:ncomp,drop=FALSE]
  
  if (type == "scores") {
    fscores
  } else if (type == "crossprod") {
    tcrossprod(fscores[,1:ncomp,drop=FALSE], x$scores[,1:ncomp,drop=FALSE])
  }else if (type =="class") {
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


#' @export
plscorr_f <- function(Y, X, ncomp=2, center=TRUE, scale=FALSE, svd.method="base") {
  assert_that(is.factor(Y))
  Xc <- scale(X, center=center, scale=scale)
  XB <- group_means(Y, Xc)
  svdres <- svd.wrapper(XB, ncomp, svd.method)

  scores <- svdres$u %*% diag(svdres$d, nrow=svdres$ncomp, ncol=svdres$ncomp)
  row.names(scores) <- levels(Y)

  ret <- list(Y=Y,X=X,ncomp=svdres$ncomp, condMeans=XB, center=center, scale=scale, pre_process=apply_scaling(Xc), 
              svd.method=svd.method, scores=scores, v=svdres$v, u=svdres$u, d=svdres$d)
  
  class(ret) <- "plscorr_result_f"
  ret
  
}

