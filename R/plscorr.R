

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
reduce <- function(x, Y, ...) UseMethod("reduce")

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
permute_refit <- function(x, ...) UseMethod("permute_refit")

#' @export
split_half_reliability <- function(x, ...) UseMethod("split_half_reliability")

#' @export
optimal_components <- function(x, ...) UseMethod("optimal_components")

#' @export
singular_values <- function(x, niter, ...) UseMethod("singular_values")

partial_scores <- function(x, ...) UseMethod("partial_scores")

contributions <- function(x, ...) UseMethod("contributions")

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



nested_cv.plscorr_result_contrast <- function(x, innerFolds, heldout, metric="AUC", min.comp=1) {
  res <- lapply(innerFolds, function(fidx) {
    print(fidx)
    exclude <- sort(c(fidx,heldout))
    prescv <- plscorr_contrast(x$X[-exclude,,drop=FALSE], x$G[-exclude,], x$strata[-exclude], ncomp=x$ncomp, 
                      center=x$center, scale=x$scale, svd.method=x$svd.method)
            
    pscores <- lapply(seq(min.comp, x$ncomp), function(n) {
      Xdat <- x$X[fidx,,drop=FALSE]
      Xb <- t(t(Xdat) %*% x$G[fidx,])
      list(pred=predict.plscorr_result_da(prescv, Xb, ncomp=n, type="prob"), obs=colnames(x$G))
    })
  })
  
  metric <- unlist(lapply(seq(min.comp, x$ncomp), function(n) {
    Plist <- lapply(res, "[[", n)
    obs <- unlist(lapply(Plist, "[[", "obs"))
    obs <- factor(obs, levels=colnames(x$G))
    pred <- do.call(rbind, lapply(Plist, "[[", "pred"))
    if (metric == "AUC") {
      combinedACC(pred, obs)
    } else if (metric == "ACC") {
      combinedAUC(pred, obs)
    }
  }))
}


#' @export
reconstruct.plscorr_result <- function(x, ncomp=x$ncomp) {
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


stratified_folds <- function(strata, nfolds=2) {
  blocks <- split(1:length(strata), strata)
  bfolds <- caret::createFolds(1:length(blocks), nfolds)
  lapply(bfolds, function(f) unlist(blocks[f]))
}

cross_validate.plscorr_result_contrast <- function(x, nfolds=2, nrepeats=10, metric=c("ACC", "distance"), nested=FALSE, max.comp=2) {
  
  M <- list()
  for (i in 1:nrepeats) {
    message("cross-validation iteration:", i)
    folds <- stratified_folds(x$strata, nfolds)
    
    ret <- lapply(1:length(folds), function(j) {
      exclude <- folds[[j]]
      prescv <- plscorr_contrast(x$X[-exclude,,drop=FALSE], x$G[-exclude,], x$strata[-exclude], ncomp=x$ncomp, 
                               center=x$center, scale=x$scale, svd.method=x$svd.method)
    
      pscores <- lapply(seq(1, max.comp), function(n) {
        Xdat <- x$X[exclude,,drop=FALSE]
        Xb <- blockwise_average(Xdat, x$G[exclude,], factor(strata[exclude]), center=FALSE, scale=FALSE)
        #Xb <- t(t(Xdat) %*% x$G[exclude,])
        list(pred=predict.plscorr_result_da(prescv, Xb, ncomp=n, type="class"), obs=colnames(x$G))
      })
    })
    
    M[[i]]  <- unlist(lapply(seq(1, max.comp), function(n) {
      Plist <- lapply(ret, "[[", n)
      obs <- unlist(lapply(Plist, "[[", "obs"))
      obs <- factor(obs, levels=colnames(x$G))
      pred <-unlist(lapply(Plist, "[[", "pred"))
      sum(obs == pred)/length(obs)
    }))
  }
  
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


plscorr_lm <- function(X, formula, design, random=NULL, ncomp=2, center=TRUE, scale=TRUE, svd.method="base") {
  tform <- terms(formula)
  facs <- attr(tform, "factors")
  
  termorder <- apply(facs,2, sum)
  orders <- seq(1, max(termorder))
  
  strata <- if (!is.null(random)) {
    design[[random]]
  }
  
  
  get_lower_form <- function(ord) {
    tnames <- colnames(facs)[termorder < ord]
    
    meat <- if (!is.null(random)) {
      paste0("(", paste(tnames, collapse= " + "), ")*", random)
    } else {
      paste(tnames, collapse= " + ")
    }
    
    as.formula(paste("X ~ ", meat))
  }
  
  fit_betas <- function(X, G) {
    if (!is.null(strata)) {
      dummy_matrix <- turner::factor_to_dummy(as.factor(strata))
      
      Xblocks <- lapply(1:ncol(dummy_matrix), function(i) {
        ind <- dummy_matrix[,i] == 1
        Xs <- X[ind,]
        lsfit(G[ind,], Xs, intercept=FALSE)$coefficients
      })
      
      Reduce("+", Xblocks)/length(Xblocks)
    
    } else {
      lsfit(G, X, intercept=FALSE)$coefficients
    }
    
  }
  
  res <- lapply(1:length(termorder), function(i) {
    ord <- termorder[i]
    if (ord == 1) {
      form <- as.formula(paste("X ~ ", colnames(facs)[i], "-1"))
      G <- model.matrix(form,data=design)
      betas <- fit_betas(X, G)
      plsres <- plscorr_contrast(betas, turner::factor_to_dummy(factor(colnames(G))), strata=NULL, 
                                 center=center, scale=scale, ncomp=ncomp, svd.method=svd.method)
      list(G=G, form=form, lower_form = ~ 1, plsres=plsres)
    } else {
      lower_form <- get_lower_form(ord)
      Glower <- model.matrix(lower_form, data=design)
      Xresid <- resid(lsfit(Glower, X, intercept=FALSE))
      form <- as.formula(paste("Xresid ~ ", colnames(facs)[i], "-1"))
      G <- model.matrix(form, data=design)
      betas <- fit_betas(Xresid, G)
      plsres <- plscorr_contrast(Xresid, turner::factor_to_dummy(factor(colnames(G))), center=center, scale=scale, ncomp=ncomp, svd.method=svd.method)
      list(G=G, form=form, lower_form = lower_form, plsres=plsres)
    }
  })
  
  names(res) <- colnames(facs)
  ret <- list(
    results=res,
    formula=formula,
    design=design,
    terms=colnames(facs)
  )
  
  
  class(ret) <- "plscorr_result_cpca"
  ret
  
}
 

#' plscorr_aov
#' 
#' @export
#' @param X the data matrix
#' @param formula a formula specifying the design
#' @param design a \code{data.frame} providing the variables provided in \code{formula} argument.
#' @param random a \code{character} string indicating the name of the random effects variable
#' @param ncomp of components to compute
plscorr_aov <- function(X, formula, design, random=NULL, ncomp=2, center=TRUE, scale=TRUE, svd.method="base") {
  tform <- terms(formula)
  facs <- attr(tform, "factors")
  
  termorder <- apply(facs,2, sum)
  orders <- seq(1, max(termorder))
  
  #if (!is.null(random)) {
  #  random2 <- update.formula(random,  ~ . -1)
  #  Grandom <- model.matrix(random, data=design)
  #  random_list <- turner::dummy_to_list(Grandom)
  #}
  
  strata <- if (!is.null(random)) {
    as.factor(design[[random]])
  } else {
    NULL
  }
  
  
  get_lower_form <- function(ord) {
    tnames <- colnames(facs)[termorder < ord]
    
    meat <- if (!is.null(random)) {
      paste0("(", paste(tnames, collapse= " + "), ")*", random)
    } else {
      paste(tnames, collapse= " + ")
    }
    
    
    as.formula(paste("X ~ ", meat))
  }
  
  res <- lapply(1:length(termorder), function(i) {
    ord <- termorder[i]
    if (ord == 1) {
      form <- as.formula(paste("X ~ ", colnames(facs)[i], "-1"))
      G <- model.matrix(form,data=design)
      cnames <- colnames(G)
      Y <- factor(cnames[apply(G, 1, function(x) which(x==1))], levels=cnames)
      plsres <- plscorr_contrast(X, G, strata, center=center, scale=scale, ncomp=ncomp, svd.method=svd.method)
      list(G=G, Y=Y, Glower=NULL, form=form, lower_form = ~ 1, plsres=plsres)
    } else {
      lower_form <- get_lower_form(ord)
      Glower <- model.matrix(lower_form, data=design)
      Xresid <- resid(lsfit(Glower, X, intercept=FALSE))
      form <- as.formula(paste("Xresid ~ ", colnames(facs)[i], "-1"))
      G <- model.matrix(form, data=design)
      cnames <- colnames(G)
      Y <- factor(cnames[apply(G, 1, function(x) which(x==1))])
      plsres <- plscorr_contrast(Xresid, G, strata=strata, center=center, scale=scale, ncomp=ncomp, svd.method=svd.method)
      list(G=G, Glower, form=form, lower_form = lower_form, plsres=plsres)
    }
  })
  
  permute <- function(obj) {
    idx <- sample(1:nrow(obj$X))
    list(X=obj$X, Y=obj$Y[idx,], idx=idx, group=group[idx])
  }
  
  names(res) <- colnames(facs)
  ret <- list(
    results=res,
    formula=formula,
    strata=strata,
    design=design,
    terms=colnames(facs)
  )
  
  
  class(ret) <- "plscorr_result_aov"
  ret
  
}

# 
#' @export
#' @import turner
plscorr_behav <- function(Y, X, group=NULL, random=NULL, ncomp=2, svd.method="base") {
    if (is.vector(Y)) {
      Y <- as.matrix(Y)
    }
    
    if (is.null(colnames(Y))) {
      colnames(Y) <- paste0("V", 1:ncol(Y))
    }
    
    assert_that(nrow(Y) == nrow(X))
    assert_that(length(group) == nrow(X))
    
    reduce <- function(obj) {
      blockids <- split(seq(1,nrow(obj$X)), obj$group)
      
      Xs <- turner::matrix_to_blocks(obj$X, blockids)
      Ys <- turner::matrix_to_blocks(obj$Y, blockids)
   
      Xs <- lapply(Xs, scale)
      Ys <- lapply(Ys, scale)
    
      R <- lapply(1:length(Xs), function(i) {
        t(cor(Xs[[i]], Ys[[i]]))
      })
    
      do.call(cbind, R)
    }
    
    bootstrap_sample <- function() {
      idx <- unlist(lapply(blockids, function(ids) sort(sample(ids, replace=TRUE))))
      YBoot <- Y[idx,]
      XBoot <- X[idx,]
      list(XBoot=XBoot, YBoot=YBoot, idx=idx, group=group[idx])
    }
    
    permute <- function(obj) {
      idx <- sample(1:nrow(obj$X))
      list(X=obj$X, Y=obj$Y[idx,], idx=idx, group=group[idx])
    }
    
    Yvars <- rep(colnames(Y), length(levels(group)))
    Gvars <- rep(levels(group), each=ncol(Y))
    
    Xred <- reduce(list(X=X, Y=Y, group=group))
    
    svdres <- svd.wrapper(t(Xred), ncomp, svd.method)
    scores <- svdres$v %*% diag(svdres$d, nrow=svdres$ncomp, ncol=svdres$ncomp)
    
    refit <- function(Y, X, group, ncomp, ...) { plscorr_behav(Y, X, group, ncomp, ...) }
      
    
    ret <- list(X=X, Y=Y, Xred=Xred, design=data.frame(Y=Yvars, group=Gvars), ncomp=svdres$ncomp, 
                svd.method=svd.method, scores=scores, v=svdres$v, u=svdres$u, d=svdres$d, 
                refit=refit, bootstrap=bootstrap_sample, reduce=reduce, permute=permute)
    
    class(ret) <- c("plscorr_result", "plscorr_result_behav")
    ret
}

blockwise_average <- function(X, G, strata, center=TRUE, scale=FALSE) {
  dummy_matrix <- turner::factor_to_dummy(as.factor(strata))
  
  Xblocks <- lapply(1:ncol(dummy_matrix), function(i) {
    ind <- dummy_matrix[,i] == 1
    scale(t(t(X[ind,]) %*% G[ind,]), center=center, scale=scale)
  })
  
  X0c <- Reduce("+", Xblocks)/length(Xblocks)
  
  if (center) {
    xc <- colMeans(do.call(rbind, lapply(Xblocks, function(x) attr(x, "scaled:center"))))
    attr(X0c, "scaled:center") <- xc
  }
  if (scale) {
    xs <- colMeans(do.call(rbind, lapply(Xblocks, function(x) attr(x, "scaled:scale"))))
    attr(X0c, "scaled:scale") <- xs
  }
  
  X0c
    
}


#' @export
#' @import turner
plscorr_contrast <- function(X, G, strata=NULL, ncomp=2, center=TRUE, scale=FALSE, svd.method="base") {
  assert_that(is.matrix(G))
  assert_that(nrow(G) == nrow(X))
  
  if (!is.null(strata)) {
    X0c <- blockwise_average(X, G, strata, center,scale)
  } else {
    X0 <- t(t(X) %*% G)
    X0c <- scale(X0, center=center, scale=scale)
  }
  
  svdres <- svd.wrapper(t(X0c), ncomp, svd.method)
  
  scores <- svdres$v %*% diag(svdres$d, nrow=svdres$ncomp, ncol=svdres$ncomp)
  row.names(scores) <- colnames(G)
  
  refit <- function(X, G, ncomp, ...) { plscorr_contrast(X, G, strata, ncomp,...) }
  
  ret <- list(X=X, G=G, ncomp=svdres$ncomp, condMeans=X0c, center=center, scale=scale, pre_process=apply_scaling(X0c), 
              svd.method=svd.method, scores=scores, v=svdres$v, u=svdres$u, d=svdres$d, refit=refit, strata=strata)
  
  class(ret) <- c("plscorr_result", "plscorr_result_contrast")
  ret
  
}



get_perm <- function(G, strata) {
  if (!is.null(x$strata)) {
    Gperm <- do.call(rbind, lapply(levels(x$strata), function(lev) {
      Gs <- G[strata==lev,]
      Gs <- Gs[sample(1:nrow(Gs)),]
    }))
  } else {
    Gperm <- x$G[sample(1:nrow(x$G)),]
  }
}

permutation_ <- function(x, nperms=100, threshold=.05, ncomp=1, verbose=TRUE, seed=NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  d <- singular_values(x)[1:ncomp]^2
  
  if (ncomp > x$ncomp) {
    ncomp <- x$ncomp
  }
  
  dstat0 <- matrix(0, nperms, ncomp)
  
  for (i in 1:nperms) {
    print(i)
    message("permutation ", i)
    fit <- x$permute_refit()
    dp <- singular_values(fit)[1:ncomp]^2
    dstat0[i, ] <- dp
  }
  
  p <- rep(0, ncomp)
  for (i in 1:ncomp) {
    p[i] <- mean(dstat0[, i] >= dstat[i])
  }
  
  for (i in 2:ncomp) {
    p[i] <- max(p[(i - 1)], p[i])
  }
  r <- sum(p <= threshold)
  return(list(perm_mat=dstat0, pval=p))

}


permutation.musubada_result <- function(x, nperms=100, threshold=.05, ncomp=1, verbose=TRUE, seed=NULL) {
  permutation_(x,nperms,threshold, ncomp, verbose, seed)
}

permutation.bada_result <- function(x, nperms=100, threshold=.05, ncomp=1, verbose=TRUE, seed=NULL) {
  permutation_(x,nperms,threshold, ncomp, verbose, seed)
}

#' @export
permutation.pls_result_contrast <- function(x, nperms=100, threshold=.05, ncomp=2, verbose=TRUE, seed=NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  dstat <- x$d[1:x$ncomp]^2/sum(x$d[1:x$ncomp]^2)
  
  dummy_matrix <- turner::factor_to_dummy(strata)
  
  subtract_recon <- function(recon) {
    ret <- lapply(1:length(Xblocks), function(i) {
      Xblocks[[i]] - recon
    })
  }
  
  permute_blocks <- function(xb, permset) {
    lapply(1:length(xb), function(i) {
      xb[[i]][permset[[i]],]
    })
  }
  
  Xblocks <- lapply(1:ncol(dummy_matrix), function(i) {
    ind <- dummy_matrix[,i] == 1
    scale(t(t(x$X[ind,]) %*% x$G[ind,]), center=x$center, scale=x$scale) 
  })
  
  if (ncomp > x$ncomp) {
    ncomp <- x$ncomp
  }
  
  dstat0 <- matrix(0, nperms, ncomp)
  
  for (i in 1:nperms) {
    
    message("permutation", i)
    for (j in 1:ncomp) {
      permset <- lapply(Xblocks, function(x) sample(1:nrow(x)))
      message("factor: ", j)
      if (j > 1) {
        message("recon")
        recon <- t(x$u[,1:(j-1),drop=FALSE] %*% t(x$scores[,1:(j-1),drop=FALSE]))
        message("subtract recon")
        Xresid <- subtract_recon(recon)
        message("permute blocks")
        Xperm <- permute_blocks(Xresid, permset)
        message("compute avg")
        Xpermavg <- Reduce("+", Xperm)/length(Xperm)
        #Xresidavg <- Reduce("+", Xresid)/length(Xresid)
        message("svd")
        dstat0[i,j] <- svd(Xpermavg)$d[1]
      } else {
        
        Xperm <- permute_blocks(Xblocks, permset)
        Xpermavg <- Reduce("+", Xperm)/length(Xperm)
        dstat0[i,j] <- svd(Xpermavg)$d[1]
      }
    }
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




