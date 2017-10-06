


## take a new design and reduce original data?

#' @importFrom assertthat assert_that
reduce_rows <- function(Xlist, Ylist) {
  
  assert_that(length(Xlist) == length(Ylist))
 
  ## compute barycenters for each table, averaging over replications
  XB <- lapply(seq_along(Xlist), function(i) group_means(Ylist[[i]], Xlist[[i]]))
  Xr <- block_matrix(XB)
}




#' musu_bada
#' 
#' @importFrom assertthat assert_that 
#' @param Y dependent \code{factor} variable. If All X matrices have same number of rows, Y can be a single factor.
#'        If there are a different number of rows (e.g. different numbers of replications per subject), Y can be a list of factors.
#' @param Xlist a list of X matrices, one per subject. 
#' @param ncomp number of common components to estimate
#' @param center center the variables
#' @param scale scale the variable by 1/sd
#' @param normalization the type of normalization
#' @param rank_k use reduce data to k components per block
#' @export
musu_bada <- function(Y, Xlist, ncomp=2, center=TRUE, scale=FALSE,  
                     normalization=c("MFA", "RV", "None", "RV-MFA"), rank_k=NULL) {


  normalization <- match.arg(normalization)

  assertthat::assert_that(all(sapply(Xlist, is.matrix)))
  
  if (is.null(names(Xlist))) {
    names(Xlist) <- paste0(1:length(Xlist))
  }
   
  table_names <- names(Xlist)
 
  ## work out the vategory structure of the X matrices
  if (is.factor(Y) && length(Y) == sum(sapply(Xlist, nrow))) {
    ## Y is a single factor with length equal to the sum of rows of each block
    Yl <- split(Y, rep(1:length(Xlist), sapply(Xlist, nrow)))
  } else if (is.factor(Y)) {
    ## Y is a single factor, therefore all matrices must have same number of rows.
    assert_that(all(sapply(Xlist, nrow) == nrow(Xlist[[1]])))
    assert_that(length(Y) == nrow(Xlist[[1]]))
    Yl <- replicate(length(Xlist), Y, simplify=FALSE)
  } else if (is.list(Y)) {
    ## Y is a list of factors
    assert_that(all(sapply(Y, is.factor)))
    Yl <- Y
    ## check that all Ys have same levels
    assert_that(length(unlist(Y)) == sum(sapply(Xlist, nrow)))
  }
  
  ## determine number of repetitions per category/block
  has_reps <- any(sapply(Yl, function(y) any(table(y) > 1)))
  
  Y_reps <- if (has_reps) {
    lapply(Yl, function(f) {
      m <- outer(f, unique(f), "==")
      apply(m * apply(m,2,cumsum), 1, sum)
    })
  } else {
    lapply(Yl, function(y) rep(1, length(y)))
  }
  
  YIndices <- rep(1:length(Xlist), sapply(Xlist, nrow))
  
  ## create block variable
  blockInd <- blockIndices(Xlist)
  
  refit <- function(.Y, .Xlist, .ncomp=ncomp) { 
    musu_bada(.Y, .Xlist, .ncomp, center, scale, normalization, rank_k) 
  }
  
  permute_refit <- function(.ncomp=ncomp) {
    nreps <- sapply(lapply(Yl, table), min)
    if (any(nreps == 1)) {
      ## permute data
      .Xlist <- lapply(Xlist, function(x) x[sample(1:nrow(x)),])
      refit(Y, .Xlist, .ncomp)
    } else {
      ## permute labels
      YPerm <- lapply(Yl, function(y) y[sample(1:length(y))])
      refit(YPerm, Xlist, .ncomp)
    }
    
  }
    


  conditions <- factor(levels(Yl[[1]]), levels=levels(Yl[[1]]))
  
  ## average categories within block to get barycenters
  Xr <- reduce_rows(Xlist, Yl)
  
  mfa_fit <- mfa(Xr, ncomp=ncomp, center=center, scale=scale, normalization=normalization, rank_k=rank_k)
  ncomp <- length(mfa_fit$ncomp)
  
 
  sc <- scores(mfa_fit)
  row.names(sc) <- levels(Yl[[1]])
  

  result <- list(
    Xlist=Xlist,
    Y=Yl,
    Y_reps=Y_reps,
    conditions=conditions,
    XB=Xr,
    scores=sc,
   
    ntables=length(Xlist),
    ncond=nrow(Xr),
    mfa_fit=mfa_fit,
    center=center,
    scale=scale,
    ncomp=mfa_fit$ncomp,
    block_indices=mfa_fit$block_indices,
    alpha=mfa_fit$alpha,
    normalization=mfa_fit$normalization,
    refit=refit,
    table_names=table_names,
    reprocess=mfa_fit$reprocess,
    rank_k=rank_k,
    permute_refit=permute_refit
  )
  
  class(result) <- c("musu_bada", "list")
  result
}

#' @export
scores.musu_bada <- function(x) {
  x$scores
}


#' @export
contributions.musu_bada <- function(x, type=c("table", "column", "row")) {
  contributions(x$mfa_fit, type)
}


#' @export
permute_refit.musu_bada <- function(x) {
  x$permute_refit()
}


#' @export
singular_values.musu_bada <- function(x) {
  x$mfa_fit$pca_fit$d
}



# project_copy.musu_bada <- function(x, Ylist, Xlist) {
#   ## reduce but no centering or scaling
#   XB <- lapply(1:length(Xlist), function(i) group_means(Ylist[[i]], Xlist[[i]]))
#   
#   XBlocks <- lapply(1:length(XB), function(i) x$reprocess(XB[[i]], i))
#   browser()
#   partial_scores = 
#     lapply(1:length(XBlocks), function(i) {
#       ind <- x$blockInd[i,1]:x$blockInd[i,2]
#       length(XBlocks) * XBlocks[[i]] %*% x$pca_fit$v[ind,]
#     })
#   
#   lds <- do.call(cbind, XBlocks %*% x$pca_fit$u)
#   
#   scores <- Reduce("+", partial_fscores)/length(partial_fscores)
#   list(scores=scores, partial_scores=partial_scores, loadings=lds)
#   
# }

#' @importFrom assertthat assert_that 
#' @export
# project_table.musu_bada <- function(x, supY, supX, ncomp=x$ncomp) { #table_index=NULL, table_name = NULL) {
#   assert_that(ncomp >= 1)
#   assert_that(is.factor(supY))
#   assert_that(length(levels(supY)) == x$ncond)
#   
#   ## pre-process new table
#   suptab <- block_reduce(list(supX), list(supY), normalization=x$normalization, scale=x$scale, center=x$center)
#   
#   ## compute supplementary Q
#   Qsup <- t(suptab$Xr) %*% (x$pca_fit$u[, 1:ncomp, drop=FALSE] %*% diag(1/x$pca_fit$d[1:ncomp]))
#   
#   ##Qsup <- t(suptab$Xr) %*% (x$pca_fit$u[, 1:ncomp, drop=FALSE] )
#   
#   ## compute supplementary factor scores
#   Fsup <- x$ntables * (suptab$Xr %*% Qsup) 
#   list(scores=Fsup, loadings=Qsup)
# }

#' @export
supplementary_loadings.musu_bada <- function(x, suptab, ncomp=x$ncomp) {
  suptab <- x$reprocess(suptab)
  Qsup <- t(suptab) %*% (x$pca_fit$u[,1:ncomp,drop=FALSE]) %*% diag(1/x$pca_fit$d[1:ncomp])
}


#' @export
loadings.musu_bada <- function(x, table_index=1:nrow(x$blockIndices), comp=1:ncol(x$pca_fit$v)) {
  do.call(cbind, lapply(table_index, function(i) {
    ind <- x$block_indices[[i]]
    1/x$alpha[i] * x$pca_fit$v[ind, comp,drop=FALSE]
  }))
}

#' @export
correlations.musu_bada <- function(x, table_index=1:nrow(x$blockIndices), comp=1:ncol(x$pca_fit$u)) {
  do.call(cbind, lapply(table_index, function(i) {
    ind <- x$blockIndices[i,1]:x$blockIndices[i,2]
    cor(x$XB[,ind,drop=FALSE], x$pca_fit$u[, comp])
  }))
}


# rename "embed_table"?

#' @export
#' @importFrom assertthat assert_that
supplementary_predictor.musu_bada <- function(x, supX, supY, type=c("class", "prob", "scores", "crossprod", "distance", "cosine"), 
                                             ncomp=x$ncomp) {
  assert_that(length(supY) == nrow(supX))
  type <- match.arg(type)
  ## pre-process new table
  suptab <- block_reduce(list(supX), list(supY), normalization=x$normalization, scale=x$scale, center=x$center)
  
  ## compute supplementary Q
  Qsup <- t(suptab$Xr) %*% (x$pca_fit$u[, 1:ncomp, drop=FALSE] %*% diag(1/x$pca_fit$d[1:ncomp]))
  
  ##Qsup <- t(suptab$Xr) %*% (x$pca_fit$u[, 1:ncomp, drop=FALSE])
  
  predfun <- function(newdat) {
    if (!is.null(suptab$centroids)) {
      newdat <- sweep(newdat, 2, suptab$centroids[[1]], "-")
    }
    if (!is.null(suptab$scales[[1]])) {
      sweep(newdat, 2, suptab$scales[[1]], "/")
    }
    
    newdat <- newdat * suptab$alpha
    fscores <- x$ntables * newdat %*% Qsup
    scorepred(fscores, x$scores, type, ncomp)
  }
}

 
#' @importFrom abind abind
project.musu_bada <- function(x, newdata, ncomp=x$ncomp, table_index=1:x$ntables) {
 
  # if (is.null(newX)) {
  #   ## project each replication onto the subject's components
  #   lapply(table_index, function(i) {
  #     yr <- x$Y_reps[[i]]
  #     ysplit <- split(1:length(yr), yr)
  #     newdata <- x$reprocess(x$Xlist[[i]],i)
  #     #browser()
  #     lv <- lapply(ysplit, function(idx) {
  #       xnewdat <- newdata[idx,]
  #       #t(t(xnewdat) %*% x$pca_fit$u[, 1:ncomp])
  #       t(t(xnewdat) %*% (x$pca_fit$u[, 1:ncomp] %*% diag(1/x$pca_fit$d[1:ncomp])))
  #     })
  #     
  #     abind::abind(lv, along=3)
  #     
  #   })
  #   
  # } else {
  
  project(x$mfa_fit, newdata, ncomp=ncomp, table_index=table_index)

}
  


## project from existing table
#' @export
predict.musu_bada <- function(x, newdata, type=c("class", "prob", "scores", "crossprod", "distance", "cosine"), 
                                    ncomp=x$ncomp, table_index=1:x$ntables, pre_process=TRUE) {
  type <- match.arg(type)

  fscores <- predict(x$mfa_fit, newdata, ncomp=ncomp, table_index=table_index, pre_process=pre_process)
  
  if (type == "scores") {
    fscores
  } else {
    scorepred(as.matrix(fscores), scores(x), type=type, ncomp=ncomp)
  }
  
}


#' @export
reconstruct.musu_bada <- function(x, ncomp=x$ncomp) {
  reconstruct(x$mfa_fit, ncomp=ncomp)
}

reproducibility.musu_bada <- function(x, blocks, nrepeats=5, metric=c("norm-2", "norm-1", "avg_cor")) {
  
  metric <- match.arg(metric)
  
  recon_error <- lapply(1:nrepeats, function(fnum) {
    message("musu_bada: reproducibility fold: ", fnum)
    fidx <- lapply(blocks, function(bind) {
      buniq <- unique(bind)
      bsam <- sample(buniq, length(buniq)/2)
      which(bind %in% bsam)
    })
    
    xsub <- musu_subset(x, lapply(fidx, "*", -1))
    rfit <- x$refit(xsub$y,xsub$x, x$ncomp) 
    xsubout <- musu_subset(x, fidx)
    
    xb <- lapply(1:length(xsubout$x), function(i) {
      yi <- xsubout$y[[i]]
      xi <- xsubout$x[[i]]
      gm <- group_means(yi, xi)
    })
    
    xb <- do.call(cbind, xb)
    
    res <- lapply(seq(0, x$ncomp), function(nc) {
      if (nc ==0) {
        switch(metric,
               "norm-2"=sqrt(sum((xb)^2)),
               "norm-1"=sum(abs(xb)),
               "avg_cor"=0)
      } else {
        xrecon <- reconstruct(rfit, ncomp=nc)
        switch(metric,
             "norm-2"=sqrt(sum((xb - xrecon)^2)),
             "norm-1"=sum(abs(xb - xrecon)),
             "avg_cor"=mean(diag(cor(t(xb), t(xrecon))))
        )
      }
    })
    
    unlist(res)
    
  })
  
  recon_error <- do.call(rbind, recon_error)
  out <- as.data.frame(recon_error)
  names(out) <- paste("NC_", seq(0, x$ncomp))
  row.names(out) <- paste0("Fold_", 1:nrow(recon_error))
  out
  
}

musu_subset <- function(x, fidx) {
  .Xlist <- lapply(1:length(x$Xlist), function(i) {
    xi <- x$Xlist[[i]]
    xi[fidx[[i]],,drop=FALSE]
  })
  
  .Y <- lapply(1:length(x$Xlist), function(i) {
    yi <- x$Y[[i]]
    yi[fidx[[i]]]
  })
  
  list(x=.Xlist, y=.Y)
  
}

performance.musu_bada <- function(x, ncomp=x$ncomp, folds=10, metric=c("ACC", "AUC")) {
  if (length(folds) == 1) {
    folds <- lapply(1:length(x$Y), function(i) caret::createFolds(x$Y[[i]], folds))
  } else {
    ## folds must be a list of blocking variables
    folds <- lapply(folds, function(bind) split(1:length(bind), bind))
  }
  
  ncomp <- min(x$ncomp, ncomp)
  
  yobs <- x$Y
  
  res <- lapply(seq_along(folds[[1]]), function(fnum) {
    message("musu_bada: performance fold: ", fnum)
    fidx <- lapply(folds, "[[", fnum)
    xsub <- musu_subset(x, lapply(fidx, "*", -1))
    
    rfit <- x$refit(xsub$y,xsub$x, ncomp) 
    
    xsubout <- musu_subset(x, fidx)
    
    preds <- lapply(1:rfit$ntables, function(i) {
      predict.musu_bada(rfit, xsubout$x[[i]], type="class", ncomp=ncomp, table_index=i)
    })
    
  })
  
  perf <- lapply(1:x$ntables, function(tind) {
    p <- unlist(lapply(res, "[[", tind))
    p == yobs[[tind]][unlist(folds[[tind]])]
  })
  
  
}

#' @importFrom vegan procrustes 
procrusteanize.musu_bada <- function(x, ncomp=2) {
  F <- scores(x)[,1:ncomp]
  
  res <- lapply(1:length(x$Xlist), function(i) {
    xp <- x$Xlist[[i]] * x$alpha[i]
    pres <- my_procrustes(F, xp)
    list(H=pres$rotation,
         scalef=pres$scale)
  })
  
  ret <- list(rot_matrices=res, ncomp=ncomp, musufit=x)
  class(ret) <- c("procrusteanized_musu_bada", "list")
  ret
}

predict.procrusteanized_musu_bada <- function(x, newdata, type=c("class", "prob", "scores", "crossprod", "distance", "cosine"), 
                              ncomp=x$ncomp, table_index=1:x$ntables, pre_process=TRUE) {
  
  
  

  type <- match.arg(type)
  
  mf <- x$musufit
  assert_that(is.matrix(newdata))
  assert_that(length(table_index) == 1 || length(table_index) == mf$ntables)
  
  fscores <- if (length(table_index) == mf$ntables) {
    assert_that(ncol(newdata) == sum(sapply(mf$Xlist, ncol)))
    Reduce("+", lapply(table_index, function(i) {
      ind <- mf$blockIndices[i,]
      
      Xp <- if (pre_process) {
        mf$reprocess(newdata[, ind[1]:ind[2]], i)
      } else {
        newdata[, ind[1]:ind[2]]
      }
      
     
      fscores <- (Xp * mf$alpha[i] * x$rot_matrices[[i]]$scalef) %*% x$rot_matrices[[i]]$H
    
      }))
  } else if (length(table_index) == 1) {
    
    ind <- mf$block_indices[[table_index]]
    
    Xp <- if (pre_process) {
      mf$reprocess(newdata, table_index)
    } else {
      newdata[, ind]
    }

    (Xp * x$rot_matrices[[table_index]]$scalef) %*% x$rot_matrices[[table_index]]$H
    
    
  }
  
  if (type == "scores") {
    fscores[,1:ncomp]
  } else {
    scorepred(fscores[, 1:ncomp], mf$scores, type=type, ncomp=ncomp)
  }
  
}


my_procrustes <- function (X, Y, scale = TRUE, symmetric = FALSE, scores = "sites", 
          ...) 
{

  if (nrow(X) != nrow(Y)) 
    stop("Matrices have different number of rows: ", nrow(X), 
         " and ", nrow(Y))
  if (ncol(X) < ncol(Y)) {
    warning("X has fewer axes than Y: X adjusted to comform Y\n")
    addcols <- ncol(Y) - ncol(X)
    for (i in 1:addcols) X <- cbind(X, 0)
  }
  ctrace <- function(MAT) sum(MAT^2)
  c <- 1
  if (symmetric) {
    X <- scale(X, scale = FALSE)
    Y <- scale(Y, scale = FALSE)
    X <- X/sqrt(ctrace(X))
    Y <- Y/sqrt(ctrace(Y))
  }
  xmean <- apply(X, 2, mean)
  ymean <- apply(Y, 2, mean)
  if (!symmetric) {
    X <- scale(X, scale = FALSE)
    Y <- scale(Y, scale = FALSE)
  }
  XY <- crossprod(X, Y)
  sol <- svd(XY)
  A <- sol$v %*% t(sol$u)
  if (scale) {
    c <- sum(sol$d)/ctrace(Y)
  }
  Yrot <- c * Y %*% A
  b <- xmean - c * ymean %*% A
  R2 <- ctrace(X) + c * c * ctrace(Y) - 2 * c * sum(sol$d)
  reslt <- list(Yrot = Yrot, X = X, ss = R2, rotation = A, 
                translation = b, scale = c, xmean = xmean, symmetric = symmetric, 
                call = match.call())
  reslt$svd <- sol
  class(reslt) <- "procrustes"
  reslt
}
  


