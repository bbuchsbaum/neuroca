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

blockIndices <- function(Xlist) {
  ncols <- sapply(Xlist, ncol)
  csum <- cumsum(ncols)
  csum1 <- c(0, csum[-length(csum)])
  m <- cbind(csum1+1, csum)
}
  
#' @importFrom assertthat assert_that 
### implement soft-thresholding that spans datasets...? similar to spls?
#' @param Y dependent \code{factor} variable. If All X matrices have same number of rows, Y can be a single factor.
#'        If there are a different nume rof rows (e..g different numbers of replications per subject), Y can be a list of factors.
#' @param Xlist a list of X matrices, one per subject. 
#' @importFrom assertthat assert_that 
#' @export
musubada <- function(Y, Xlist, ncomp=2, center=TRUE, scale=FALSE, svd.method="fast", normalization=c("MFA", "RV", "None")) {
  normalization <- normalization[1]
  
  assert_that(all(sapply(Xlist, is.matrix)))
  
  message("svd.method", svd.method)
  
  print("expanding Y")
  if (is.factor(Y)) {
    ## Y is a single factor, therefore all matrices must have same number of rows.
    assert_that(all(sapply(Xlist, nrow) == nrow(Xlist[[1]])))
    Yl <- replicate(length(Xlist), Y, simplify=FALSE)
  } else if (is.list(Y)) {
    assert_that(all(sapply(Y, is.factor)))
    Yl <- Y
    assert_that(length(unlist(Y)) == sum(sapply(Xlist, nrow)))
  }
  print("done")
  print("Y indices")
  YIndices <- rep(1:length(Xlist), sapply(Xlist, nrow))
  
  
  print("done")
  print("blockVar")
  ## create block variable
  blockInd <- blockIndices(Xlist)
  print("done")
  
  pre_process <- function(X) {
    Xc <- scale(X, center=center, scale=scale)
    
    if (normalization == "MFA") {
      alpha <- 1/svd.wrapper(X, ncomp=1, method="propack")$d[1]
      Xc <- Xc * alpha
    } else if (normalization == "None") {
      alpha <- 1
    } else if (normalization == "RV" || "VoxelCor") {
      alpha <- 1
    } else {
      stop(paste0("normalization of type: '", normalization, "' is not supported."))
    }
    
    attr(Xc, "alpha") <- alpha
    
    applyFun <- function(X) {
      if (center) {
        X <- sweep(X, 2, attr(Xc, "scaled:center"), "-")
      }
      if (scale) {
        X <- sweep(X, 2, attr(Xc, "scaled:scale"), "/")
      }
      
      X * alpha
    }
    
    attr(Xc, "applyFun") <- applyFun
   
    Xc
  }
  
  ## take a new design and reduce original data
  reduce <- function(Xlist) {
    ## compute barycenters for each table
    print("group_means")
    XB <- lapply(1:length(Xlist), function(i) group_means(Yl[[i]], Xlist[[i]]))
    print("done")
    # center/scale barcycenters of each table
    
    print("pre_process")
    normXBc <- lapply(XB, pre_process)
    
    if (normalization == "RV") {
      message("rv normalization")
      rvmat <- RVMat(normXBc)
      diag(rvmat) <- 1
      wts <- abs(svd.wrapper(rvmat, ncomp=1, method="propack")$u[,1])
      alpha <- wts/sum(wts)
      normXBc <- lapply(1:length(normXBc), function(i) normXBc[[i]] * alpha[i])
    } else if (normalization == "VoxelCor") {
      rmat <- corMat(normXBc)
      rmat <- rmat + min(rmat)
      diag(rmat) <- 1
      wts <- abs(svd.wrapper(rmat, ncomp=1, method="propack")$u[,1])
      alpha <- wts/sum(wts)
      normXBc <- lapply(1:length(normXBc), function(i) normXBc[[i]] * alpha[i])
    }
    

    print("done")
    
    print("cbind")
    Xr <- do.call(cbind, normXBc)
    print("done")
    list(Xred=Xr, XB=XB, alpha=sapply(normXBc, function(x) attr(x, "alpha")))
  }
  
  refit <- function(.Y, .Xlist, .ncomp) { musubada(.Y, .Xlist, .ncomp, center, scale, svd.method, normalization) }
  
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
  pls_res <- plscorr_da(YB, Xred$Xred, ncomp=ncomp, center=FALSE, scale=FALSE, svd.method=svd.method)
 
  
  result <- list(
    Y=Y,
    scores=pls_res$scores,
    partial_fscores = 
      lapply(1:length(Xlist), function(i) {
        ind <- blockInd[i,1]:blockInd[i,2]
        length(Xlist) * normXBc[[i]] %*% pls_res$u[ind,]
      }),
    
    table_contr = do.call(cbind, lapply(1:ncomp, function(i) {
      sapply(1:length(Xlist), function(j) {
        ind <- blockInd[j,1]:blockInd[j,2]
        sum(pls_res$u[ind,i]^2)
      })
    })),
    ntables=length(Xlist),
    ncond=nrow(XB[[1]]),
    plsfit=pls_res,
    center=center,
    scale=scale,
    ncomp=ncomp,
    blockIndices=blockInd,
    blockVar=blockVar,
    alpha=Xred$alpha,
    pre_process=pre_process,
    refit=refit,
    permute_refit=permute_refit
  )
  
  class(result) <- c("pca_result", "musubada_result")
  result
}


permute_refit.musubada_result <- function(x) {
  x$permute_refit()
}

singular_values.musubada_result <- function(x) {
  x$plsfit$d
}


#' @export
project_table <- function(x, newdata, ncomp, ...) UseMethod("project_table")


#' @importFrom assertthat assert_that 
#' @export
project_table.musubada_result <- function(x, suptab, ncomp=x$ncomp, table_index=NULL) {
  assert_that(is.matrix(suptab))
  assert_that(nrow(suptab) == x$ncond )
  assert_that(ncomp >= 1)
  
  suptab <- x$pre_process(suptab)
  #Qsup <- t(suptab) %*% (x$plsfit$v[,1:ncomp,drop=FALSE] %*% diag(1/x$plsfit$d[1:ncomp]))
  #Qsup <- t(suptab) %*% (x$plsfit$v[,1:ncomp,drop=FALSE] %*% matrix(1/x$plsfit$d[1:ncomp], ncomp, nrow(Xr), byrow=TRUE)
  #Qsup <- t(suptab) %*% (x$scores)    
  Qsup <- t(suptab) %*% (x$plsfit$v)  
  Fsup <- x$ntables * (suptab %*% Qsup) 
}

supplementary_loadings.musubada_result <- function(x, suptab, ncomp=x$ncomp) {
  suptab <- x$pre_process(suptab)
  Qsup <- t(suptab) %*% (x$plsfit$v[,1:ncomp,drop=FALSE])
}


#' @export
loadings.musubada_result <- function(x, table_index) {
 ind <- x$blockIndices[table_index,1]:x$blockIndices[table_index,2]
 x$plsfit$u[ind, ,drop=FALSE]
}


#' @export
supplementary_predictor.musubada_result <- function(x, suptab, type=c("class", "prob", "scores", "crossprod"), ncomp=x$ncomp) {
  suptab <- x$pre_process(suptab)
  Qsup <- t(suptab) %*% (x$plsfit$v[,1:ncomp,drop=FALSE])
  
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
    loadings <- x$plsfit$u[ind,1:ncomp]
    Xp %*% loadings
  } else {
    Xp %*% x$plsfit$u[,1:ncomp]
  }
  
  scorepred(fscores, x$scores, ncomp)
  
}
