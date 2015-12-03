
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
scores <- function(x,...) UseMethod("scores")

#' @export
loadings <- function(x,...) UseMethod("loadings")

#' @export
cross_validate <- function(x, ...) UseMethod("cross_validate")

#' @export
nested_cv <- function(x, ...) UseMethod("nested_cv")

#' @export
bootstrap <- function(x, niter, ...) UseMethod("bootstrap")

#' @export
jackstraw <- function(x, nsynth, niter, ...) UseMethod("jackstraw")

#' @export
permutation <- function(x, ...) UseMethod("permutation")

#' @export
optimal_components <- function(x, ...) UseMethod("optimal_components")

#' @export
reproducibility <- function(x, folds, metric, ...) UseMethod("reproducibility")

#' @export
project.rows <- function(xr, u, ncomp) {
  xr %*% u 
}

#' @export
project.cols <- function(xc, v, ncomp) {
  xc %*% v
}

#' @export
reconstruct <- function(x, ncomp) UseMethod("reconstruct")

#' @export
svd.wrapper <- function(XC, ncomp=min(dim(XC)), method=c("base", "fast", "irlba")) {
  assert_that(method %in% c("base", "fast", "irlba"))
  res <- switch(method[1],
                base=svd(XC),
                fast=corpcor:::fast.svd(XC),
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
  ret <- sweep(Rs, 1, yt, "/")
  row.names(ret) <- names(yt)
  ret
}

#' @export
plscorr <- function(Y, X, ncomp=2, center=TRUE, scale=FALSE, svd.method="fast.svd") {
  if (is.factor(Y)) {
    plscorr_da(Y, X, ncomp, svd.method)
  } else {
    stop("Y variable must be factor")
  }
}

#' @export
scores.pls_result_da <- function(x) {
  x$scores
}


#' @export
nested_cv.plscorr_result_da <- function(x, innerFolds, heldout, metric="AUC", min.comp=1) {
  res <- lapply(innerFolds, function(fidx) {
    exclude <- sort(c(fidx,heldout))
    prescv <- x$refit(x$Y[-exclude], x$X[-exclude,,drop=FALSE], ncomp=x$ncomp, 
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
reconstruct.plscorr_result_da <- function(x, ncomp=x$ncomp) {
  t(x$u[,1:ncomp,drop=FALSE] %*% t(x$scores[,1:ncomp,drop=FALSE]))
}

#' @export
reproducibility.pls_result_da <- function(x, folds, metric=c("norm-2", "norm-1", "avg_cor")) {
  if (length(folds) == 1) {
    folds <- caret::createFolds(1:length(x$Y), folds)
  } else if (length(folds) == length(x$Y)) {
    folds <- split(1:length(x$Y), folds)
  }
  
  metric <- metric[1]
  
  preds <- lapply(seq_along(folds), function(fnum) {
    message("plscorr_da: reproducibility iteration: ", fnum)
    fidx <- folds[[fnum]]
    pmod <- plscorr_da(x$Y[-fidx], x$X[-fidx,,drop=FALSE], ncomp=x$ncomp, svd.method=x$svd.method)
    
    xb <- group_means(x$Y[fidx], x$X[fidx,,drop=FALSE])
    xb <- x$pre_process(xb)
    
    res <- lapply(seq(1, x$ncomp), function(nc) {
      xrecon <- reconstruct(x, ncomp=nc)
      switch(metric,
            "norm-2"=sqrt(sum((xb - xrecon)^2)),
            "norm-1"=sum(abs(xb - xrecon)),
            "avg_cor"=mean(diag(cor(t(xb), t(xrecon))))
      )
    })
  })
    
           
}


#' @export
cross_validate.plscorr_result_da <- function(x, folds, metric="AUC") {
  if (length(folds) == 1) {
    folds <- caret::createFolds(1:length(x$Y), folds)
  } else if (length(folds) == length(x$Y)) {
    folds <- split(1:length(x$Y), folds)
  }
  
  preds <- lapply(seq_along(folds), function(fnum) {
    message("plscorr_da: cross validation iteration: ", fnum)
    res <- nested_cv(x, folds[-fnum], folds[[fnum]], min.comp=1)
    
    perc_of_best <- res/max(res)
    print(perc_of_best)
    nc <- which(perc_of_best > .95)[1]
    
    fidx <- folds[[fnum]]
    pmod <- x$refit(x$Y[-fidx], x$X[-fidx,,drop=FALSE], ncomp=nc, svd.method=x$svd.method) 
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
scorepred <- function(fscores, scores, type=c("class", "prob", "scores", "crossprod"), ncomp=2) {
  if (type == "scores") {
    fscores
  } else if (type == "crossprod") {
    tcrossprod(fscores[,1:ncomp,drop=FALSE], scores[,1:ncomp,drop=FALSE])
  }else if (type =="class") {
    D <- rdist(fscores[,1:ncomp,drop=FALSE], scores[,1:ncomp,drop=FALSE])
    D2 <- D^2
    min.d <- apply(D2, 1, which.min)
    row.names(scores)[min.d]
  } else if (type == "prob"){
    ## type is 'prob'
    probs <- tcrossprod(fscores[,1:ncomp,drop=FALSE], scores[,1:ncomp,drop=FALSE])
    maxid <- apply(probs,1,which.max)
    maxp <- probs[cbind(1:nrow(probs), maxid)]
    probs <- exp(probs - maxp)
    probs <- zapsmall(probs/rowSums(probs))
  } else {
    stop(paste("illegal 'type' argument: ", type))
  }
}

#' @export
predict.plscorr_result_da <- function(x, newdata, type=c("class", "prob", "scores", "crossprod"), ncomp=NULL) {
  if (is.vector(newdata) && (length(newdata) == ncol(x$condMeans))) {
    x <- matrix(newdata, nrow=1, ncol=ncol(x$condMeans))
  }
  
  assert_that(is.matrix(newdata))
  assert_that(ncol(newdata) == ncol(x$condMeans))
  assert_that(type %in% c("class", "prob", "scores", "crossprod"))
  
  if (is.null(ncomp)) {
    ncomp <- x$ncomp
  }
  
  xc <- x$pre_process(as.matrix(newdata))
  
  #fscores <- xc %*% x$v[,1:ncomp,drop=FALSE]
  fscores <- xc %*% x$u[,1:ncomp,drop=FALSE]
  scorepred(fscores, x$scores,type, ncomp)
  
}

#' @export
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
plscorr_sda <- function(Y, X, ncomp=2, center=TRUE, scale=FALSE, svd.method="base") {
  assert_that(is.factor(Y))
  
  if (any(table(Y) > 1)) {
    ## centroids
    XB <- sda(X, Y, diagonal=TRUE)
    XBc <- XB$beta
    global_mean <- colMeans(X)
  } else {
    stop("must have more than one observation per category")
  }
  
  do_scale <- function(X) {
    sweep(X, 2, global_mean)
  }
  
  #svdres <- svd.wrapper(XBc, ncomp, svd.method)
  svdres <- svd.wrapper(t(XBc), ncomp, svd.method)
  
  scores <- svdres$v %*% diag(svdres$d, nrow=svdres$ncomp, ncol=svdres$ncomp)
  row.names(scores) <- levels(Y)
  
  refit <- function(Y, X, ncomp, ...) { plscorr_sda(Y, X, ncomp,...) }
  
  ret <- list(Y=Y,X=X,ncomp=svdres$ncomp, condMeans=XBc, center=center, scale=scale, pre_process=do_scale, 
              svd.method=svd.method, scores=scores, v=svdres$v, u=svdres$u, d=svdres$d, refit=refit)
  
  class(ret) <- c("plscorr_result_da", "plscorr_result_sda")
  ret
  
}


#' @export
plscorr_da <- function(Y, X, ncomp=2, center=TRUE, scale=FALSE, svd.method="base") {
  assert_that(is.factor(Y))
  
  if (any(table(Y) > 1)) {
    ## compute barycenters
    XB <- group_means(Y, X)
    XBc <- scale(XB, center=center, scale=scale)
  } else {
    XB <- X
    XBc <- scale(XB, center=center, scale=scale)
  }
  
  #svdres <- svd.wrapper(XBc, ncomp, svd.method)
  svdres <- svd.wrapper(t(XBc), ncomp, svd.method)
  
  scores <- svdres$v %*% diag(svdres$d, nrow=svdres$ncomp, ncol=svdres$ncomp)
  row.names(scores) <- levels(Y)

  refit <- function(Y, X, ncomp, ...) { plscorr_da(Y, X, ncomp,...) }
  
  ret <- list(Y=Y,X=X,ncomp=svdres$ncomp, condMeans=XBc, center=center, scale=scale, pre_process=apply_scaling(XBc), 
              svd.method=svd.method, scores=scores, v=svdres$v, u=svdres$u, d=svdres$d, refit=refit)
  
  class(ret) <- "plscorr_result_da"
  ret
  
}

#' @export
scores.pls_result_da <- function(x) {
  x$scores
}

#' @export
optimal_components.pls_result_da <- function(x, method="smooth") {
  FactoMineR::estim_ncp(t(x$condMeans),method=method)
}


#' @export
permutation.pls_result_da <- function(x, nperms=100, threshold=.05, verbose=TRUE, seed=NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  dstat <- x$d[1:x$ncomp]^2/sum(x$d[1:x$ncomp]^2)
  dstat0 <- matrix(0, nperms, length(dstat))
  
  for (i in 1:nperms) {
    print(i)
    yperm <- x$Y[sample(1:length(x$Y))]
    pres <- plscorr_da(yperm, x$X, x$ncomp, x$center, x$scale, svd.method=x$svd.method) 
    dstat0[i, ] <- pres$d[1:length(dstat)]^2/sum(pres$d[1:length(dstat)]^2)
  }
  
  p <- rep(0, x$ncomp)
  for (i in 1:x$ncomp) {
    p[i] <- mean(dstat0[, i] >= dstat[i])
  }
  for (i in 2:x$ncomp) {
    p[i] <- max(p[(i - 1)], p[i])
  }
  r <- sum(p <= threshold)
  return(list(r = r, p = p))
  
}

#' @export
bootstrap.plscorr_result_da <- function(x, niter, strata=NULL) {
  do_boot <- function(indices) {
    XBoot <- x$X[indices,]
    YBoot <- x$Y[indices]
    
    #Xc <- x$pre_process(XBoot)
    
    ## barycenters of bootstrap samples
    XB <- group_means(YBoot, XBoot)
    
    ## apply pre-processing from full model
    XBc <- x$pre_process(XB) 
    
    br <- project.rows(XBc, x$u)
    bc <- project.cols(t(XBc), x$v)
    
    list(boot4R=br,
         boot4C=bc)
  }
  
  boot_ratio <- function(bootlist) {
    boot.mean <- Reduce("+", bootlist)/length(bootlist)
    boot.sd <- sqrt(Reduce("+", lapply(bootlist, function(mat) (mat - boot.mean)^2))/length(bootlist))
    boot.mean/boot.sd	
  }
  
  .resample <- function(x, ...)  { x[sample.int(length(x), ...)] }
  
  sample_indices <- function(Y, strata=NULL) {
    if (is.null(strata)) {
      sort(unlist(lapply(levels(x$Y), function(lev) {
        sample(which(x$Y == lev), replace=TRUE)
      })))
    } else {
      boot.strat <- sample(levels(strata), replace=TRUE)
      indices <- lapply(boot.strat, function(stratum) {
        unlist(lapply(levels(Y), function(lev) {
          .resample(which(Y == lev & strata==stratum), replace=TRUE)					
        }))
      })
    }
  }
  
  boot.res <- 
    lapply(1:niter, function(i) {
      message("bootstrap iteration: ", i)
      row_indices <- sample_indices(x$Y, strata)
      do_boot(sort(unlist(row_indices)))
    })
 
  
  ret <- list(boot.ratio.R=boot_ratio(lapply(boot.res, "[[", 1)),
              boot.ratio.C=boot_ratio(lapply(boot.res, "[[", 2)),
              boot.raw.R=lapply(boot.res, "[[", 1),
              boot.raw.C=lapply(boot.res, "[[", 2))
  
  class(res) <- c("list", "bootstrap_result")
  
}


