#' @export
supplementary_loadings <- function(x,...) UseMethod("supplementary_loadings")


colCors = function(x, y) { 
  sqr = function(x) x*x
  if (!is.matrix(x)||!is.matrix(y)||any(dim(x)!=dim(y))) {
    stop("Please supply two matrices of equal size.")
  }
  
  x   <- sweep(x, 2, colMeans(x))
  y   <-  sweep(y, 2, colMeans(y))
  cor <- colSums(x*y) /  sqrt(colSums(sqr(x))*colSums(sqr(y)))
  return(cor)
}

corMat <- function(Xl) {
  pairs <- combn(length(Xl),2)
  M <- matrix(0, length(Xl), length(Xl))
  for (i in 1:ncol(pairs)) {
    p1 <- pairs[1,i]
    p2 <- pairs[2,i]
    rmean <- mean(colCors(t(Xl[[p1]]), t(Xl[[p2]])))
    M[p1,p2] <- rmean
    M[p2,p1] <- rmean
  }
  M
  
}

SimMat <- function(Xlist, FUN, ...) {
  pairs <- combn(length(Xlist),2)
  M <- matrix(0, length(Xlist), length(Xlist))
  for (i in 1:ncol(pairs)) {
    p1 <- pairs[1,i]
    p2 <- pairs[2,i]
    rv <- FUN(Xlist[[p1]], Xlist[[p2]], ...)
    M[p1,p2] <- rv
    M[p2,p1] <- rv
  }
  M
}

RVMat <- function(Xlist) {
  pairs <- combn(length(Xlist),2)
  M <- matrix(0, length(Xlist), length(Xlist))
  for (i in 1:ncol(pairs)) {
    p1 <- pairs[1,i]
    p2 <- pairs[2,i]
    rv <- coefficientRV(Xlist[[p1]], Xlist[[p2]], center=FALSE, scale=FALSE)
    M[p1,p2] <- rv
    M[p2,p1] <- rv
  }
  M
  
}
coefficientRV <- function(X, Y, center=TRUE, scale=FALSE) {
  if (dim(X)[[1]] != dim(Y)[[1]]) 
    stop("no the same dimension for X and Y")
  if (dim(X)[[1]] == 1) {
    print("1 configuration RV is  NA")
    rv = NA
  }
  else {
    if (center || scale) {
      Y <- scale(Y, center=center, scale=scale)
      X <- scale(X, center=center, scale=scale)
    }
    W1 <- X %*% t(X)
    W2 <- Y %*% t(Y)
    rv <- sum(diag(W1 %*% W2))/(sum(diag(W1 %*% W1)) * 
                                  sum(diag(W2 %*% W2)))^0.5
  }
  return(rv)
}

normalization_factors <- function(Xl, type=c("MFA", "RV", "DCor", "None")) {
  type <- match.arg(type)
  
  alpha <- if (type == "MFA") {
    sapply(Xl, function(X) 1/(svd.wrapper(X, ncomp=1, method="propack")$d[1]^2))
  } else if (type == "RV" && length(Xl) > 2) {
    message("rv normalization")
    smat <- SimMat(Xl, function(x1,x2) coefficientRV(x1,x2, center=FALSE, scale=FALSE))
    diag(smat) <- 1
    wts <- abs(svd.wrapper(smat, ncomp=1, method="propack")$u[,1])
    wts/sum(wts)
    
    ## make alpha average 1
    (wts - mean(wts)) + 1
    
  } else if (type == "DCor" && length(Xl) > 2) {
    message("dcor normalization")
    smat <- SimMat(Xl, function(x1,x2) energy::dcor.ttest(x1,x2)$estimate)
    diag(smat) <- 1
    wts <- abs(svd.wrapper(smat, ncomp=1, method="propack")$u[,1])
    wts/sum(wts)
    (wts - mean(wts)) + 1
    ## make alpha average 1
  } else {
    rep(1, length(Xl))
  }
  
  alpha <- (alpha - mean(alpha)) + 1
  
  
}


## take a new design and reduce original data
#' @importFrom assertthat assert_that
block_reduce <- function(Xlist, Ylist, normalization=c("MFA", "RV", "DCor", "None", "Pre"), center=TRUE, scale=FALSE, alpha=rep(1, length(Xlist))) {
  assert_that(length(Xlist) == length(Ylist))
  normalization <- match.arg(normalization)
  
  ## compute barycenters for each table
  XB <- lapply(1:length(Xlist), function(i) group_means(Ylist[[i]], Xlist[[i]]))
  
  # center/scale barcycenters of each table
  XBc <- lapply(XB, pre_process, center, scale)
  
  centroids <- lapply(XBc, attr, "scaled:center")
  scales <- lapply(XBc, attr, "scaled:scale")
  
  if (normalization != "Pre") {
    # compute normalization factor
    alpha <- sqrt(normalization_factors(XBc, normalization))
  }
  
  # normalize
  normXBc <- lapply(1:length(alpha), function(i) XBc[[i]] * alpha[i])
  
  # construct matrix of barycenters
  Xr <- do.call(cbind, normXBc)
  
  list(Xr=Xr, XB=XB, alpha=alpha, centroids=centroids, scales=scales)
}



### implement soft-thresholding that spans datasets...? similar to spls?
### musubada is really a bada with block structure -- bada can be engine
#' @importFrom assertthat assert_that 
#' @importFrom energy dcor.ttest
#' @param Y dependent \code{factor} variable. If All X matrices have same number of rows, Y can be a single factor.
#'        If there are a different nume rof rows (e..g different numbers of replications per subject), Y can be a list of factors.
#' @param Xlist a list of X matrices, one per subject. 
#' @importFrom assertthat assert_that 
#' @export
musubada <- function(Y, Xlist, ncomp=2, center=TRUE, scale=FALSE, svd.method="fast", 
                     normalization=c("MFA", "RV", "None","DCor")) {
  normalization <- normalization[1]

  assert_that(all(sapply(Xlist, is.matrix)))
  
  if (is.null(names(Xlist))) {
    names(Xlist) <- paste0(1:length(Xlist))
  }
   
  table_names <- names(Xlist)
 
  if (is.factor(Y)) {
    ## Y is a single factor, therefore all matrices must have same number of rows.
    assert_that(all(sapply(Xlist, nrow) == nrow(Xlist[[1]])))
    Yl <- replicate(length(Xlist), Y, simplify=FALSE)
  } else if (is.list(Y)) {
    ## Y is a list of factors
    assert_that(all(sapply(Y, is.factor)))
    Yl <- Y
    ## check that all Ys have same levels
    assert_that(length(unlist(Y)) == sum(sapply(Xlist, nrow)))
  }
  
  Y_reps <- 
    lapply(Yl, function(f) {
      m <- outer(f, unique(f), "==")
      apply(m* apply(m,2,cumsum), 1, sum)
    })
  
  YIndices <- rep(1:length(Xlist), sapply(Xlist, nrow))
  
  ## create block variable
  blockInd <- blockIndices(Xlist)
  
  refit <- function(.Y, .Xlist, .ncomp) { 
    musubada(.Y, .Xlist, .ncomp, center, scale, svd.method, normalization) 
  }
  
  permute_refit <- function(.ncomp=ncomp) {
    nreps <- sapply(lapply(Yl, table), min)
    if (any(nreps == 1)) {
      ## permute data
      .Xlist <- lapply(Xlist, function(x) x[sample(1:nrow(x)),])
      refit(Y, .Xlist, .ncomp)
    } else {
      YPerm <- lapply(Yl, function(y) y[sample(1:length(y))])
      refit(YPerm, Xlist, .ncomp)
    }
    
  }
    
  Xreduced <- block_reduce(Xlist, Yl, normalization, center, scale)
  
  reprocess <- function(newdat, table_index) {
    ## given a new observation(s), pre-process it in the same way the original observations were processed
    if (!is.null(Xreduced$centroids)) {
      ## recenter
      newdat <- sweep(newdat, 2, Xreduced$centroids[[table_index]], "-")
    }
    if (!is.null(Xreduced$scales[[table_index]])) {
      ## rescale
      sweep(newdat, 2, Xreduced$scales[[table_index]], "/")
    }
    
    newdat * Xreduced$alpha[table_index]
    
  }
  
  
  YB <- factor(levels(Yl[[1]]), levels=levels(Yl[[1]]))
  
  pca_fit <- pca_core(Xreduced$Xr, ncomp=ncomp, 
                      center=FALSE, 
                      scale=FALSE, 
                      svd.method=svd.method)
  
  ncomp <- length(pca_fit$d)
  
 
  partial_fscores = 
    lapply(1:length(Xlist), function(i) {
      ind <- blockInd[i,1]:blockInd[i,2]
      length(Xlist) * Xreduced$Xr[, ind] %*% pca_fit$v[ind,]
    })
  
  sc <- pca_fit$scores
  row.names(sc) <- levels(Yl[[1]])
  
  result <- list(
    Xlist=Xlist,
    Y=Yl,
    Y_reps=Y_reps,
    conditions=YB,
    XB=Xreduced$Xr,
    scores=sc,
    partial_scores=partial_fscores,
    table_contr = do.call(cbind, lapply(1:ncomp, function(i) {
      sapply(1:length(Xlist), function(j) {
        ind <- blockInd[j,1]:blockInd[j,2]
        sum(pca_fit$v[ind,i]^2)
      })
    })),
    
    ntables=length(Xlist),
    ncond=nrow(Xreduced$Xr),
    pca_fit=pca_fit,
    center=center,
    scale=scale,
    ncomp=ncomp,
    blockIndices=blockInd,
    alpha=Xreduced$alpha,
    normalization=normalization,
    refit=refit,
    table_names=table_names,
    reprocess=reprocess,
    permute_refit=permute_refit
  )
  
  class(result) <- c("musubada")
  result
}

contributions.musubada <- function(x, fac, component=1) {
  out <- lapply(1:x$ntables, function(j) {
    ind <- x$blockInd[j,1]:x$blockInd[j,2]
    contr <- pca_fit$v[ind,component]^2 * x$alpha[j]
    tapply(contr,fac, sum)
  })
  
  do.call(cbind, out)
  
}

permute_refit.musubada <- function(x) {
  x$permute_refit()
}

singular_values.musubada <- function(x) {
  x$pca_core$d
}


#' @export
project_table <- function(x, supY, supX, ncomp, ...) UseMethod("project_table")

project_copy <- function(x, ...) UseMethod("project_copy")

project_copy.musubada <- function(x, Ylist, Xlist) {
  ## reduce but no centering or scaling
  XB <- lapply(1:length(Xlist), function(i) group_means(Ylist[[i]], Xlist[[i]]))
  
  XBlocks <- lapply(1:length(XB), function(i) x$reprocess(XB[[i]], i))
  browser()
  partial_scores = 
    lapply(1:length(XBlocks), function(i) {
      ind <- x$blockInd[i,1]:x$blockInd[i,2]
      length(XBlocks) * XBlocks[[i]] %*% x$pca_fit$v[ind,]
    })
  
  lds <- do.call(cbind, XBlocks %*% x$pca_fit$u)
  
  scores <- Reduce("+", partial_fscores)/length(partial_fscores)
  list(scores=scores, partial_scores=partial_scores, loadings=lds)
  
}

#' @importFrom assertthat assert_that 
#' @export
project_table.musubada <- function(x, supY, supX, ncomp=x$ncomp) { #table_index=NULL, table_name = NULL) {
  assert_that(ncomp >= 1)
  assert_that(is.factor(supY))
  assert_that(length(levels(supY)) == x$ncond)
  
  ## pre-process new table
  suptab <- block_reduce(list(supX), list(supY), normalization=x$normalization, scale=x$scale, center=x$center)
  
  ## compute supplementary Q
  Qsup <- t(suptab$Xr) %*% (x$pca_fit$u[, 1:ncomp, drop=FALSE] %*% diag(1/x$pca_fit$d[1:ncomp]))
  
  ##Qsup <- t(suptab$Xr) %*% (x$pca_fit$u[, 1:ncomp, drop=FALSE] )
  
  ## compute supplementary factor scores
  Fsup <- x$ntables * (suptab$Xr %*% Qsup) 
  list(scores=Fsup, loadings=Qsup)
}

#' @export
supplementary_loadings.musubada <- function(x, suptab, ncomp=x$ncomp) {
  suptab <- x$pre_process(suptab)
  Qsup <- t(suptab) %*% (x$plsfit$u[,1:ncomp,drop=FALSE]) %*% diag(1/x$pca_fit$d[1:ncomp])
}


#' @export
loadings.musubada <- function(x, table_index=1:nrow(x$blockIndices), comp=1:ncol(x$pca_fit$v)) {
  do.call(cbind, lapply(table_index, function(i) {
    ind <- x$blockIndices[i,1]:x$blockIndices[i,2]
    1/x$alpha[i] * x$pca_fit$v[ind, comp,drop=FALSE]
  }))
}

#' @export
correlations.musubada <- function(x, table_index=1:nrow(x$blockIndices), comp=1:ncol(x$pca_fit$u)) {
  do.call(cbind, lapply(table_index, function(i) {
    ind <- x$blockIndices[i,1]:x$blockIndices[i,2]
    cor(x$XB[,ind,drop=FALSE], x$pca_fit$u[, comp])
  }))
}


# rename "embed_table"?

#' @export
#' @importFrom assertthat assert_that
supplementary_predictor.musubada <- function(x, supX, supY, type=c("class", "prob", "scores", "crossprod", "distance", "cosine"), 
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

 
#' @import From abind abind
project_cols.musubada <- function(x, newdata=NULL, ncomp=x$ncomp, table_index=1:x$ntables) {
  assert_that(length(table_index) == 1 || length(table_index) == x$ntables)
  
  if (is.null(newdata)) {
    ## project each replication onto the subject's components
    lapply(table_index, function(i) {
      yr <- x$Y_reps[[i]]
      ysplit <- split(1:length(yr), yr)
      newdata <- x$reprocess(x$Xlist[[i]],i)
      #browser()
      lv <- lapply(ysplit, function(idx) {
        xnewdat <- newdata[idx,]
        #t(t(xnewdat) %*% x$pca_fit$u[, 1:ncomp])
        t(t(xnewdat) %*% (x$pca_fit$u[, 1:ncomp] %*% diag(1/x$pca_fit$d[1:ncomp]))
      })
      
      abind::abind(lv, along=3)
      
    })
    
  } else {
    stop()
    ## project new data-point from existing table.
    lapply(1:length(table_index), function(i) {
      browser()
      tbind <- table_index[i]
      ind <- x$blockIndices[tbind,]
      xnewdat <- x$reprocess(newdata[[i]], i)
      t(t(xnewdat) %*% x$pca_fit$u[ind[1]:ind[2], 1:ncomp])
    })
  }
}
  


## project from existing table
#' @export
predict.musubada <- function(x, newdata, type=c("class", "prob", "scores", "crossprod", "distance"), 
                                    ncomp=x$ncomp, table_index=1:x$ntables) {
  type <- match.arg(type)
  assert_that(is.matrix(newdata))
  assert_that(length(table_index) == 1 || length(table_index) == x$ntables)
  
  fscores <- if (length(table_index) == x$ntables) {
    assert_that(ncol(newdata) == sum(sapply(x$Xlist, ncol)))
    Reduce("+", lapply(table_index, function(i) {
      ind <- x$blockIndices[i,]
      Xp <- x$reprocess(newdata[, ind[1]:ind[2]], i)
      fscores <- Xp %*% x$pca_fit$v[ind[1]:ind[2],,drop=FALSE]
    }))
  } else if (length(table_index) == 1) {
    ind <- x$blockIndices[table_index,]
    Xp <- x$reprocess(newdata, table_index)
    Xp %*% x$pca_fit$v[ind[1]:ind[2], , drop=FALSE]
  }
  
  
    
  scorepred(fscores, x$scores, type=type, ncomp=ncomp)
  
}
