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
  
 
  if (is.factor(Y)) {
    ## Y is a single factor, therefore all matrices must have same number of rows.
    assert_that(all(sapply(Xlist, nrow) == nrow(Xlist[[1]])))
    Yl <- replicate(length(Xlist), Y, simplify=FALSE)
  } else if (is.list(Y)) {
    ## Y is a list of factors
    assert_that(all(sapply(Y, is.factor)))
    Yl <- Y
    assert_that(length(unlist(Y)) == sum(sapply(Xlist, nrow)))
  }
  
  YIndices <- rep(1:length(Xlist), sapply(Xlist, nrow))
  
  ## create block variable
  blockInd <- blockIndices(Xlist)
  
  normalization <- function(X) {
    if (normalization == "MFA") {
      alpha <- 1/svd.wrapper(X, ncomp=1, method="propack")$d[1]
    } else {
      1
    }
  }
 
  pre_process <- function(X) {
    Xc <- scale(X, center=center, scale=scale)
    
    applyFun <- function(X) {
      if (center) {
        X <- sweep(X, 2, attr(Xc, "scaled:center"), "-")
      }
      if (scale) {
        X <- sweep(X, 2, attr(Xc, "scaled:scale"), "/")
      }
      
    }
    
    attr(Xc, "applyFun") <- applyFun
    Xc
  }
  
  ## take a new design and reduce original data
  reduce <- function(Xlist) {
    
    ## compute barycenters for each table
    XB <- lapply(1:length(Xlist), function(i) group_means(Yl[[i]], Xlist[[i]]))
    
    # center/scale barcycenters of each table
    normXBc <- lapply(XB, pre_process)
    
    # compute normalization factor
    alpha <- sapply(XB, normalization)
    
    if (normalization == "MFA") {
      normXBc <- lapply(1:length(alpha), function(i) normXBc[[i]] * alpha[i])
    } else if (normalization == "RV") {
      message("rv normalization")
      smat <- SimMat(normXBc, function(x1,x2) coefficientRV(x1,x2, center=FALSE, scale=FALSE))
      diag(smat) <- 1
      wts <- abs(svd.wrapper(smat, ncomp=1, method="propack")$u[,1])
      alpha <- wts/sum(wts)
      normXBc <- lapply(1:length(normXBc), function(i) normXBc[[i]] * alpha[i])
      attr(normXBc, "alpha") <- alpha
      normXBC
    } else if (normalization == "DCor") {
      message("dcor normalization")
      smat <- SimMat(normXBc, function(x1,x2) energy::dcor.ttest(x1,x2)$estimate)
      diag(smat) <- 1
      wts <- abs(svd.wrapper(smat, ncomp=1, method="propack")$u[,1])
      alpha <- wts/sum(wts)
      normXBc <- lapply(1:length(normXBc), function(i) normXBc[[i]] * alpha[i])
      attr(normXBc, "alpha") <- alpha
      normXBc
    }

    Xr <- do.call(cbind, normXBc)
    list(Xred=Xr, XB=XB, alpha=sapply(normXBc, function(x) attr(x, "alpha")))
  }
  
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
    
  Xred <- reduce(Xlist)
  YB <- factor(row.names(normXBc[[1]]), levels=row.names(normXBc[[1]]))
  
  pca_fit <- pca_core(Xred$Xred, ncomp=ncomp, 
                      center=FALSE, 
                      scale=FALSE, 
                      svd.method=svd.method)
 
  result <- list(
    Y=Y,
    scores=pca_fit$scores,
    partial_fscores = 
      lapply(1:length(Xlist), function(i) {
        ind <- blockInd[i,1]:blockInd[i,2]
        length(Xlist) * normXBc[[i]] %*% pca_fit$v[ind,]
      }),
    
    table_contr = do.call(cbind, lapply(1:ncomp, function(i) {
      sapply(1:length(Xlist), function(j) {
        ind <- blockInd[j,1]:blockInd[j,2]
        sum(pca_fit$v[ind,i]^2)
      })
    })),
    
    ntables=length(Xlist),
    ncond=nrow(XB[[1]]),
    pca_fit=pca_fit,
    center=center,
    scale=scale,
    ncomp=ncomp,
    blockIndices=blockInd,
    alpha=Xred$alpha,
    pre_process=pre_process,
    refit=refit,
    permute_refit=permute_refit
  )
  
  class(result) <- c("musubada_result")
  result
}

contributions.musubada_result <- function(x, fac, component=1) {
  out <- lapply(1:x$ntables, function(j) {
    ind <- x$blockInd[j,1]:x$blockInd[j,2]
    contr <- pca_fit$v[ind,component]^2 * x$alpha[j]
    tapply(contr,fac, sum)
  })
  
  do.call(cbind, out)
  
}

permute_refit.musubada_result <- function(x) {
  x$permute_refit()
}

singular_values.musubada_result <- function(x) {
  x$pca_core$d
}


#' @export
project_table <- function(x, newdata, ncomp, ...) UseMethod("project_table")


#' @importFrom assertthat assert_that 
#' @export
project_table.musubada_result <- function(x, suptab, ncomp=x$ncomp, table_index=NULL) {
  assert_that(is.matrix(suptab))
  assert_that(nrow(suptab) == x$ncond )
  assert_that(ncomp >= 1)
  
  ## pre-process new table
  suptab <- x$pre_process(suptab)
  #Qsup <- t(suptab) %*% (x$plsfit$v[,1:ncomp,drop=FALSE] %*% diag(1/x$plsfit$d[1:ncomp]))
  #Qsup <- t(suptab) %*% (x$plsfit$v[,1:ncomp,drop=FALSE] %*% matrix(1/x$plsfit$d[1:ncomp], ncomp, nrow(Xr), byrow=TRUE)
  #Qsup <- t(suptab) %*% (x$scores)    
  Qsup <- t(suptab) %*% (x$plsfit$u)  
  Fsup <- x$ntables * (suptab %*% Qsup) 
}

#' @export
supplementary_loadings.musubada_result <- function(x, suptab, ncomp=x$ncomp) {
  suptab <- x$pre_process(suptab)
  Qsup <- t(suptab) %*% (x$plsfit$u[,1:ncomp,drop=FALSE])
}


#' @export
loadings.musubada_result <- function(x, table_index=1) {
 ind <- x$blockIndices[table_index,1]:x$blockIndices[table_index,2]
 x$pca_fit$v[ind, ,drop=FALSE]
}


#' @export
supplementary_predictor.musubada_result <- function(x, suptab, type=c("class", "prob", "scores", "crossprod"), ncomp=x$ncomp) {
  suptab <- x$pre_process(suptab)
  Qsup <- t(suptab) %*% (x$pca_fit$u[,1:ncomp,drop=FALSE])
  
  predfun <- function(newdat) {
    pfun <- attr(suptab, "applyFun")
    newdat <- pfun(newdat)
    fscores <- x$ntables * newdat %*% Qsup
    scorepred(fscores, x$scores, type, ncomp)
  }
}


#' @export
predict.musubada_result <- function(x, newdata, type=c("class", "prob", "scores", "crossprod"), ncomp=x$ncomp, table_index=NULL) {
  #if (is.vector(newdata) && (length(newdata) == ncol(x$condMeans))) {
  #  x <- matrix(newdata, nrow=1, ncol=ncol(x$condMeans))
  #
  #}
  
  assert_that(is.matrix(newdata))
  #assert_that(nrow(newdata) == x$ncond)
  assert_that(type[1] %in% c("class", "prob", "scores", "crossprod"))
  
  Xp <- x$pre_process(newdata)
  
  fscores <- if (!is.null(table_index)) {
    ind <- x$blockIndices[[table_index]]
    loadings <- x$pca_fit$u[ind,1:ncomp]
    Xp %*% loadings
  } else {
    Xp %*% x$pca_fit$u[,1:ncomp]
  }
  
  scorepred(fscores, x$scores, ncomp)
  
}
