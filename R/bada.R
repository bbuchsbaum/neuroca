#' @export
#' @importFrom assertthat assert_that
bada <- function(Y, X, ncomp=length(levels(Y)), center=TRUE, scale=FALSE, svd.method="base", strata=NULL) {
  assert_that(is.factor(Y))
  assert_that(length(Y) == nrow(X)) 
  
  if (!is.null(strata)) {
    assert_that(length(strata) == length(Y))
    strata <- as.factor(strata)
  }
  
  reduce <- function(.Y, .X, .strata, .center, .scale) {
    if (any(table(.Y) > 1) && is.null(.strata)) {
      ## compute barycenters, no strata
      XB <- group_means(.Y, .X)
      list(XBc=scale(XB, center=.center, scale=.scale), XBlock=NULL)
    } else if (!is.null(strata)) {
      ## compute barycenters by strata
      XBlock <- lapply(levels(.strata), function(i) {
        idx <- which(.strata==i)
        scale(group_means(.Y[idx], .X[idx,]), center=.center, scale=.scale)
      })
      
      list(XBc=Reduce("+", XBlock)/length(XBlock), XBlock=XBlock)
    } else {
      list(XBc=scale(.X, center=.center, scale=.scale), XBlock=NULL)
    }
    
  }
  
  xb <- reduce(Y, X, strata, center, scale)
  pcres <- pca_core(t(xb$XBc), ncomp, center=FALSE, scale=FALSE, svd.method)
  
  scores <- pcres$scores
  row.names(scores) <- levels(Y)
  
  refit <- function(.Y, .X, .ncomp, .strata) { bada(.Y, .X, .ncomp, center, scale, svd.method, .strata) }
  
  permute_refit <- function(.ncomp=ncomp) {
    if (is.null(strata)) {
      idx <- sample(seq_along(Y))
      refit(Y[idx], X, .ncomp)
    } else {
      Xperm <- do.call(rbind, lapply(levels(strata), function(lev) {
        idx <- which(lev == strata)
        X[sample(idx),]
      }))
      refit(Y, Xperm, .ncomp)
    }
  }
  
  ret <- list(Y=Y,X=X,ncomp=ncomp, condMeans=xb$XBc, center=center, scale=scale, pre_process=apply_scaling(xb$XBc), 
              svd.method=svd.method, scores=scores, v=pcres$v, u=pcres$u, d=pcres$d, refit=refit, permute_refit=permute_refit, reduce=reduce,
              strata=strata, XBlock=xb$XBlock)
 
  class(ret) <- c("pca_result", "bada_result")
  ret
  
}


singular_values.bada_result <- function(x) {
  x$d
}

partial_scores.bada_result <- function(x) {
  if (is.null(x$strata)) {
    x$scores
  } else {
    levs <- levels(x$strata)
    lapply(1:length(levs), function(i) {
      x$XBlock[[i]] %*% x$u
    })
  }
}

#' @export
scores.bada_result <- function(x) {
  x$scores
}

loadings.bada_result <- function(x) {
  x$v
}

#' @export
#' @importFrom MatrixCorrelation PCAcv
optimal_components.bada_result <- function(x, ncomp=x$ncomp) {
  if (ncol(x$condMeans) > 5000) {
    ssize <- as.integer(min(c(.1 * ncol(x$condMeans)), 5000))
    iter <- min(round(2*ncol(x$condMeans)/ssize), 100)
    res <- lapply(1:iter, function(i) {
      idx <- sample(1:ncol(bres$condMeans), ssize)
      MatrixCorrelation::PCAcv(x$condMeans[,idx], ncomp)
    })
    
    colMeans(do.call(rbind, res))
  } else {
    MatrixCorrelation::PCAcv(x$condMeans, ncomp)
  }
}

#' @export
cross_validate.bada_result <- function(x, folds, metric="AUC") {
  if (length(folds) == 1) {
    folds <- caret::createFolds(1:length(x$Y), folds)
  } else if (length(folds) == length(x$Y)) {
    folds <- split(1:length(x$Y), folds)
  } else {
    stop("folds variable must be a single integer or a vector of indices equal to the length of x$Y")
  }
  
  preds <- lapply(seq_along(folds), function(fnum) {
    message("bada: cross validation iteration: ", fnum)
    res <- nested_cv(x, folds[-fnum], folds[[fnum]], "AUC", min.comp=1)
    
    perc_of_best <- res/max(res)
    print(perc_of_best)
    nc <- which(perc_of_best > .95)[1]
    
    fidx <- folds[[fnum]]
    pmod <- x$refit(x$Y[-fidx], x$X[-fidx,,drop=FALSE], .ncomp=nc, .strata=droplevels(x$strata[-fidx])) 
    pclass <- predict(pmod, x$X[folds[[fnum]],,drop=FALSE], ncomp=nc, type="class")
    d <- predict(pmod, x$X[folds[[fnum]],,drop=FALSE], ncomp=nc, type="distance")
    list(class=pclass,dist=d,K_perf=res, K_best=nc)
  })
  
  dmat <- do.call(rbind, lapply(preds, "[[", "dist"))
  class <- unlist(lapply(preds, "[[", "class"))
  ncomp <- unlist(lapply(preds, "[[", "K_best"))
  perfmat <- do.call(rbind, lapply(preds, "[[", "K_perf"))
  list(class=class, dist=dmat, ncomp=ncomp, ncomp_perf=perfmat)
}

#' @export
nested_cv.bada_result <- function(x, innerFolds, heldout, metric=c("ACC", "AUC"), min.comp=1) {
  metric <- match.arg(metric)
  res <- lapply(innerFolds, function(fidx) {
    
    exclude <- sort(c(fidx,heldout))
    prescv <- x$refit(x$Y[-exclude], x$X[-exclude,,drop=FALSE], .ncomp=x$ncomp, .strata=droplevels(x$strata[-exclude]))
                      
    distlist <- lapply(seq(min.comp, x$ncomp), function(n) {
      predict(prescv, x$X[fidx,,drop=FALSE], ncomp=n, type="distance")
    })
    
  })
  
  innerIdx <- unlist(innerFolds)
  obs <- x$Y[innerIdx]
  obs <- factor(obs, levels=rownames(x$scores))
  
  p <- unlist(lapply(seq(min.comp, x$ncomp), function(n) {
    
    Dn <- lapply(res, "[[", n)
    
    if (metric == "ACC") {
      preds <- unlist(lapply(1:length(Dn), function(i) {
        min.d <- apply(Dn[[i]], 1, which.min)
        row.names(scores)[min.d]
      }))
      
      sum(preds == obs)/length(innerIdx)
      
    } else if (metric == "AUC") {
      preds <- do.call(rbind, lapply(1:length(Dn), function(i) {
        1/(Dn[[i]]+1)
      }))
      combinedAUC(preds, obs)
    }
      
  }))
}


#' @export
predict.bada_result <- function(x, newdata, type=c("class", "prob", "scores", "cosine", "distance"), ncomp=x$ncomp, pre_process=TRUE) {
  assert_that(type %in% c("class", "prob", "scores", "crossprod", "distance"))
  if (is.vector(newdata) && (length(newdata) == ncol(x$condMeans))) {
    newdata <- matrix(newdata, nrow=1, ncol=ncol(x$condMeans))
  }
  assert_that(is.matrix(newdata))
  assert_that(ncol(newdata) == ncol(x$condMeans))
  type <- match.arg(type)
  
 
  if (is.null(ncomp)) {
    ncomp <- x$ncomp
  }
  
  xc <- if (pre_process) {
    x$pre_process(as.matrix(newdata))
  } else {
    newdata
  }
  
  scorepred(fscores, x$scores[, 1:ncomp,drop=FALSE],type, ncomp)
  
}

#' @export
cross_validate.bada_result <- function(x, folds, metric="AUC") {
  if (length(folds) == 1) {
    folds <- caret::createFolds(1:length(x$Y), folds)
  } else if (length(folds) == length(x$Y)) {
    folds <- split(1:length(x$Y), folds)
  }
  
  preds <- lapply(seq_along(folds), function(fnum) {
    message("cross validation iteration: ", fnum)
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
bootstrap.plscorr_result_da <- function(x, niter) {
  do_boot <- function(indices) {
    XBoot <- x$X[indices,]
    YBoot <- x$Y[indices]
    
    #Xc <- x$pre_process(XBoot)
    
    ## barycenters of bootstrap samples
    #XB <- x$reduce(XBoot, YBoot)
    XB <- group_means(YBoot, XBoot)
    
    ## apply pre-processing from full model
    XBc <- x$pre_process(XB) 
    
    br <- project.rows(XBc, x$u)
    bc <- project.cols(t(XBc), x$v)
    
    list(boot4R=br,
         boot4C=bc)
  }
  
  do_stratified_boot <- function(ids) {
    blocks <- x$XBlock[ids]
    XB <- Reduce("+", blocks)/length(blocks)
    br <- project.rows(XB, x$u)
    bc <- project.cols(t(XB), x$v)
    
    list(boot4R=br,
         boot4C=bc)
    
  }
  
  boot_ratio <- function(bootlist) {
    boot.mean <- Reduce("+", bootlist)/length(bootlist)
    boot.sd <- sqrt(Reduce("+", lapply(bootlist, function(mat) (mat - boot.mean)^2))/length(bootlist))
    boot.mean/boot.sd	
  }
  
  .resample <- function(x, ...)  { x[sample.int(length(x), ...)] }
  
  sample_indices <- function(Y) {
    sort(unlist(lapply(levels(x$Y), function(lev) {
      sample(which(x$Y == lev), replace=TRUE)
    })))
  }
  
  boot.res <- 
    lapply(1:niter, function(i) {
      message("bootstrap iteration: ", i)
      if (is.null(x$strata)) {
        row_indices <- sample_indices(x$Y)
        do_boot(sort(unlist(row_indices)))
      } else {
        levs <- levels(x$strata)
        slevs <- sample(levs, length(levs), replace=TRUE)
        ind <- sort(sapply(slevs, function(lev) which(lev == levs)))
        ret <- do_stratified_boot(ind)
      }
    })
  
  
  ret <- list(boot.ratio.R=boot_ratio(lapply(boot.res, "[[", 1)),
              boot.ratio.C=boot_ratio(lapply(boot.res, "[[", 2)),
              boot.raw.R=lapply(boot.res, "[[", 1),
              boot.raw.C=lapply(boot.res, "[[", 2))
  
  class(res) <- c("list", "bootstrap_result")
  
}






