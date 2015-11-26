#' @export
supplementary_loadings <- function(x,...) UseMethod("supplementary_loadings")


coefficientRV <- function(X, Y) {
  if (dim(X)[[1]] != dim(Y)[[1]]) 
    stop("no the same dimension for X and Y")
  if (dim(X)[[1]] == 1) {
    print("1 configuration RV is  NA")
    rv = NA
  }
  else {
    Y <- scale(Y, scale = FALSE)
    X <- scale(X, scale = FALSE)
    W1 <- X %*% t(X)
    W2 <- Y %*% t(Y)
    rv <- sum(diag(W1 %*% W2))/(sum(diag(W1 %*% W1)) * 
                                  sum(diag(W2 %*% W2)))^0.5
  }
  return(rv)
}



#' @importFrom assertthat assert_that 
### implement soft-thresholding that spans datasets...? similar to spls?
#' @importFrom assertthat assert_that 
#' @export
musubada <- function(Y, Xlist, ncomp=2, center=TRUE, scale=FALSE, svd.method="fast.svd", normalization=c("MFA", "RV", "None")) {
  normalization <- normalization[1]
  
  assert_that(all(sapply(Xlist, is.matrix)))
  #assert_that(all(sapply(Xlist, nrow) == nrow(Xlist[[1]])))
  
  assert_that(is.factor(Y))
  
  YIndices <- rep(1:length(Xlist), sapply(Xlist, nrow))
  
  ## create block variable
  blockVar <- factor(unlist(lapply(1:length(Xlist), function(i) rep(i, ncol(Xlist[[i]])))))
  
  ## get column indices associated with each block
  blockIndices <- lapply(1:length(Xlist), function(i) which(blockVar == i))
  
  pre_process <- function(X) {
    Xc <- scale(X, center=center, scale=scale)
    
    if (normalization == "MFA") {
      alpha <- 1/svd.wrapper(X, ncomp=1, method="irlba")$d[1]
      Xc <- Xc * alpha
    } else if (normalization == "None") {
      alpha <- 1
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
  
  ## compute barycenters for each table
  XB <- lapply(1:length(Xlist), function(i) group_means(Y[YIndices==i], Xlist[[i]]))
  
  # center/scale barcycenters of each table
  normXBc <- lapply(XB, pre_process)
  
  
  
  ## compute scaling factor based on reciprocal of first squared singular value
  #alpha <- lapply(XBc, function(X) 1/svd.wrapper(X, ncomp=1, method="irlba")$d[1])
  
  ## scale each table by 'alpha'
  #normXBc <- lapply(1:length(XBc), function(i) XBc[[i]] * alpha[[i]])
  
  YB <- factor(row.names(normXBc[[1]]))
  pls_res <- plscorr_da(YB, do.call(cbind, normXBc), ncomp=ncomp, center=FALSE, scale=FALSE)
  
  
  result <- list(
    scores=pls_res$scores,
    partial_fscores = 
      lapply(1:length(Xlist), function(i) {
        ind <- blockIndices[[i]]
        length(Xlist) * normXBc[[i]] %*% pls_res$u[ind,]
      }),
    
    table_contr = do.call(cbind, lapply(1:ncomp, function(i) {
      sapply(1:length(blockIndices), function(j) {
        sum(pls_res$u[blockIndices[[j]],i]^2)
      })
    })),
    ntables=length(Xlist),
    ncond=nrow(XB[[1]]),
    plsfit=pls_res,
    center=center,
    scale=scale,
    ncomp=ncomp,
    blockIndices=blockIndices,
    blockVar=blockVar,
    alpha=sapply(normXBc, function(x) attr(x, "alpha")),
    pre_process=pre_process
  )
  
  class(result) <- "musubada_result"
  result
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
 x$plsfit$u[x$blockIndices[[table_index]], ,drop=FALSE]
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
